#include <zmq.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/thread/thread.hpp>

#include <chrono>
// #include <complex>
#include <csignal>
#include <fstream>
#include <iostream>
#include <thread>

// VRT
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include <vrt/vrt_read.h>
#include <vrt/vrt_string.h>
#include <vrt/vrt_types.h>
#include <vrt/vrt_util.h>

// #include <complex.h>
#include <complex>

#include "vrt-tools.h"
#include "tracker-extended-context.h"

const double pi = std::acos(-1.0);
const std::complex<double> complexi(0.0, 1.0);

namespace po = boost::program_options;

static bool stop_signal_called = false;
void sig_int_handler(int)
{
    stop_signal_called = true;
}

template <typename samp_type> inline float get_abs_val(samp_type t)
{
    return std::fabs(t);
}

inline float get_abs_val(std::complex<int16_t> t)
{
    return std::fabs(t.real());
}

inline float get_abs_val(std::complex<int8_t> t)
{
    return std::fabs(t.real());
}

int main(int argc, char* argv[])
{

    // variables to be set by po
    std::string file, type, zmq_address;
    uint16_t pub_instance, instance, main_port, port, pub_port;
    uint32_t channel;
    int hwm;
    float freq_offset, bandwidth, doppler_rate;
    double frequency;
    size_t num_requested_samples;
    double total_time;

    uint32_t decimation;
    uint32_t taps_per_decimation;
    uint32_t num_taps;
    double *taps;
    float **poly_taps;

    std::complex<double> f0;
    std::complex<double> alpha;
    std::complex<double> alpha_dop;
    std::complex<float>  alpha2;
    std::complex<double> step;
    std::complex<double> step_dop;
    float polyfir_channel;

    std::complex<float>*x;
    std::complex<float>*y;
    std::complex<float>*tmp_acc;

    // setup the program options
    po::options_description desc("Allowed options");
    // clang-format off

    desc.add_options()
        ("help", "help message")
        ("nsamps", po::value<size_t>(&num_requested_samples)->default_value(0), "total number of samples to receive")
        ("duration", po::value<double>(&total_time)->default_value(0), "total number of seconds to receive")
        ("channel", po::value<uint32_t>(&channel)->default_value(0), "VRT channel")
        ("progress", "periodically display short-term bandwidth")
        ("int-second", "align start of reception to integer second")
        ("null", "run without writing to file")
        ("continue", "don't abort on a bad packet")
        ("channel-mode", "use frequency/offset to select channel")
        ("tracking", "use VRT tracking data")
        ("decimation", po::value<uint32_t>(&decimation)->default_value(2), "decimation factor")
        ("taps-per-decimation", po::value<uint32_t>(&taps_per_decimation)->default_value(20), "taps per decimation")
        ("bandwidth", po::value<float>(&bandwidth)->default_value(0), "bandwidth")
        ("doppler", po::value<float>(&doppler_rate)->default_value(0), "doppler rate in Hz/s")
        ("freq-offset", po::value<float>(&freq_offset)->default_value(0), "frequency offset")
        ("frequency", po::value<double>(&frequency)->default_value(0), "center frequency")
        ("address", po::value<std::string>(&zmq_address)->default_value("localhost"), "VRT ZMQ address")
        ("zmq-split", "create a ZeroMQ stream per VRT channel, increasing port number for additional streams")
        ("instance", po::value<uint16_t>(&instance)->default_value(0), "VRT ZMQ instance")
        ("port", po::value<uint16_t>(&port), "VRT ZMQ port")
        ("pub-port", po::value<uint16_t>(&pub_port), "VRT ZMQ PUB port")
        ("pub-instance", po::value<uint16_t>(&pub_instance)->default_value(1), "VRT ZMQ instance")
        ("hwm", po::value<int>(&hwm)->default_value(10000), "VRT ZMQ HWM")
    ;
    // clang-format on
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // print the help message
    if (vm.count("help")) {
        std::cout << boost::format("VRT channelizer. %s") % desc << std::endl;
        std::cout << std::endl
                  << "This application channelizer a VRT stream.\n"
                  << std::endl;
        return ~0;
    }

    bool progress               = vm.count("progress") > 0;
    bool stats                  = vm.count("stats") > 0;
    bool null                   = vm.count("null") > 0;
    bool continue_on_bad_packet = vm.count("continue") > 0;
    bool int_second             = (bool)vm.count("int-second");
    bool zmq_split              = vm.count("zmq-split") > 0;
    bool channel_mode           = vm.count("channel-mode") > 0;
    bool tracking               = vm.count("tracking") > 0;

    context_type vrt_context;
    init_context(&vrt_context);
    tracker_ext_context_type tracker_ext_context;

    packet_type vrt_packet;

    if (vm.count("port") > 0) {
        main_port = port;
    } else {
        main_port = DEFAULT_MAIN_PORT + MAX_CHANNELS*instance;
    }

    if (!(vm.count("pub-port") > 0)) {
        pub_port = DEFAULT_MAIN_PORT + MAX_CHANNELS*pub_instance;
    }

    if (zmq_split) {
        main_port += channel;
        vrt_packet.channel_filt = 1;
    } else {
        vrt_packet.channel_filt = 1<<channel;
    }

    // ZMQ
    void *context = zmq_ctx_new();
    void *subscriber = zmq_socket(context, ZMQ_SUB);
    int rc = zmq_setsockopt (subscriber, ZMQ_RCVHWM, &hwm, sizeof hwm);
    std::string connect_string = "tcp://" + zmq_address + ":" + std::to_string(main_port);
    rc = zmq_connect(subscriber, connect_string.c_str());
    assert(rc == 0);
    zmq_setsockopt(subscriber, ZMQ_SUBSCRIBE, "", 0);

    void *responder = zmq_socket(context, ZMQ_PUB);
    rc = zmq_setsockopt (responder, ZMQ_SNDHWM, &hwm, sizeof hwm);
    assert(rc == 0);
    connect_string = "tcp://*:" + std::to_string(pub_port);
    rc = zmq_bind(responder, connect_string.c_str());
    assert (rc == 0);

    // time keeping
    auto start_time = std::chrono::steady_clock::now();
    auto stop_time = start_time + std::chrono::milliseconds(int64_t(1000 * total_time));

    uint32_t rx_buffer[ZMQ_BUFFER_SIZE];
    uint32_t tx_buffer[ZMQ_BUFFER_SIZE];

    unsigned long long num_total_samps = 0;

    // Track time and samps between updating the BW summary
    auto last_update                     = start_time;
    unsigned long long last_update_samps = 0;

    bool first_frame = true;
    uint64_t last_fractional_seconds_timestamp = 0;

    // set to true to process data before context
    bool start_rx = false;

    uint32_t signal_pointer = 0;

    /* VRT init */
    struct vrt_packet p;
    vrt_init_packet(&p);
    vrt_init_data_packet(&p);
    p.fields.stream_id = 1;

    std::complex<int16_t> iq_buff[VRT_SAMPLES_PER_PACKET];
    uint32_t iq_counter = 0;
    uint32_t fir_pointer = 0;
    uint32_t frame_count = 0;
    uint32_t t_samp = 0;

    std::complex<double> phasor = 1;

    double total_phase = 0;

    while (not stop_signal_called
           and (num_requested_samples > num_total_samps or num_requested_samples == 0)
           and (total_time == 0.0 or std::chrono::steady_clock::now() <= stop_time)) {

        int len = zmq_recv(subscriber, rx_buffer, ZMQ_BUFFER_SIZE, 0);

        const auto now = std::chrono::steady_clock::now();

        if (not vrt_process(rx_buffer, sizeof(rx_buffer), &vrt_context, &vrt_packet)) {
            printf("Not a Vita49 packet?\n");
            continue;
        }

        if (not start_rx and vrt_packet.context
            and not (tracking and not tracker_ext_context.tracker_ext_context_received)) {
            vrt_print_context(&vrt_context);
            start_rx = true;

            if (tracking) {
                frequency = tracker_ext_context.frequency+tracker_ext_context.doppler;
                doppler_rate = tracker_ext_context.doppler_rate;
                printf("# Setting freq. to %f Hz with %f Hz/s dopppler rate for \"%s\" (source \"%s\")\n",
                    frequency, tracker_ext_context.doppler_rate, tracker_ext_context.object_name, tracker_ext_context.tracking_source);
            }

            if (frequency > 0) {
                freq_offset = frequency-(double)vrt_context.rf_freq;
            }

            if (bandwidth > 0) {
                decimation = vrt_context.sample_rate/bandwidth;
            } else {
                bandwidth = vrt_context.sample_rate/decimation;
            }

            // check for valid frequency offset
            if ( fabs(freq_offset) > vrt_context.sample_rate/2) {
                printf("Selected frequency outside of stream\n");
                exit(1);
            }

            // check for valid decimation
            if (VRT_SAMPLES_PER_PACKET % decimation != 0) {
                printf("decimation needs to be a divisor of %u.\n", VRT_SAMPLES_PER_PACKET);
                exit(1);
            }

            if ((uint64_t)vrt_context.sample_rate % decimation != 0) {
                printf("decimation needs to be a divisor of the sample rate (%u).\n", vrt_context.sample_rate);
                exit(1);
            }

            // create FIR filter
            uint32_t fir_order = taps_per_decimation*decimation-1;
            num_taps = fir_order+1;
            taps = (double*)malloc(sizeof(double)*num_taps);

            double K = 0.97*(fir_order/decimation);

            // Blackman window
            // double a0 = 0.42;
            // double a1 = 0.50;
            // double a2 = 0.08;
            // double a3 = 0.00;

            // Blackman-Harris window
            double a0 = 0.35875;
            double a1 = 0.48829;
            double a2 = 0.14128;
            double a3 = 0.01168;

            for (int i=0;i<fir_order;i++) {
                int j = -(i - fir_order/2);
                double blackman_window = a0 - a1*cos(2*pi*(double)i/((double)fir_order-1)) +
                                            a2*cos(4*pi*(double)i/((double)fir_order-1)) +
                                            a3*cos(6*pi*(double)i/((double)fir_order-1));
                if (j==0) {
                    taps[i] = blackman_window*((double)K/(double)fir_order);
                } else {
                    taps[i] = blackman_window*(1.0/(double)fir_order)*sin(pi*(double)j*(double)K/(double)fir_order)/sin(pi*(double)j/(double)fir_order);
                }
            }
            taps[fir_order] = 0;

            // Create polyphase partitions of filter
            // float poly_taps[decimation][taps_per_decimation];

            poly_taps = (float **)malloc(sizeof(float *)*decimation);
            for (size_t dec=0; dec < decimation; dec++)
                poly_taps[dec] = (float*)malloc(sizeof(float)*taps_per_decimation);

            for (uint32_t i=0; i<num_taps; i++) {
                poly_taps[i%decimation][i/decimation] = (float)taps[i];
            }

            if (channel_mode) {
                polyfir_channel = round(freq_offset/bandwidth);
                printf("# Selected channel: %.0f\n", polyfir_channel);
                freq_offset = polyfir_channel*bandwidth;
                printf("# New offset: %.0f Hz\n", freq_offset);
            } else {
                polyfir_channel = 0;
            }

            f0 = -freq_offset;
            alpha = complexi*2.0*pi*f0;
            alpha_dop = complexi*2.0*pi*(double)-doppler_rate;
            step = std::exp(alpha/(double)vrt_context.sample_rate);
            step_dop = std::exp(alpha_dop/pow((double)vrt_context.sample_rate,2));
            alpha2 = (std::complex<float>)complexi*polyfir_channel*2.0f*(float)pi/(float(decimation));

            int M = decimation;
            int L = vrt_packet.num_rx_samps;

            x = (std::complex<float>*)malloc(sizeof(std::complex<float>)*(M+L+num_taps));
            y = (std::complex<float>*)malloc(sizeof(std::complex<float>)*(L/M));
            tmp_acc = (std::complex<float>*)malloc(sizeof(std::complex<float>)*(L/M));

            for (uint32_t i = 0; i < M+L+num_taps; i++)
                x[i] = std::complex<float>(0,0);
        }

        if (start_rx and vrt_packet.context) {
            // construct new context
            struct vrt_packet pc;
            vrt_init_packet(&pc);
            vrt_init_context_packet(&pc);

            pc.fields.stream_id = vrt_context.stream_id;
            pc.fields.integer_seconds_timestamp = vrt_context.integer_seconds_timestamp;
            pc.fields.fractional_seconds_timestamp = vrt_context.fractional_seconds_timestamp;

            double doppler_offset = total_phase/(double)vrt_context.sample_rate;
            pc.if_context.rf_reference_frequency = (double)vrt_context.rf_freq+(double)freq_offset-doppler_offset;

            pc.if_context.bandwidth = vrt_context.bandwidth;
            pc.if_context.sample_rate = vrt_context.sample_rate/decimation;
            pc.if_context.rf_reference_frequency_offset = 0;
            pc.if_context.if_reference_frequency = 0;
            pc.if_context.if_band_offset = 0;
            pc.if_context.gain.stage1 = vrt_context.gain;
            pc.if_context.gain.stage2 = 0;

            pc.if_context.state_and_event_indicators.has.reference_lock = true;
            pc.if_context.state_and_event_indicators.reference_lock = vrt_context.reflock;
            pc.if_context.state_and_event_indicators.has.calibrated_time = true;
            pc.if_context.state_and_event_indicators.calibrated_time = vrt_context.time_cal;

             // TODO: check if present
            pc.if_context.has.temperature  = true;
            pc.if_context.temperature = vrt_context.temperature;
            pc.if_context.has.timestamp_calibration_time = true;
            pc.if_context.timestamp_calibration_time = vrt_context.timestamp_calibration_time;


            int32_t rv = vrt_write_packet(&pc, tx_buffer, VRT_DATA_PACKET_SIZE, true);
            if (rv < 0) {
                fprintf(stderr, "Failed to write packet: %s\n", vrt_string_error(rv));
            }

            // ZMQ
            zmq_send (responder, tx_buffer, rv*4, 0);

        }

        if (start_rx and vrt_packet.data) {

            if (vrt_packet.lost_frame)
               if (not continue_on_bad_packet)
                    break;

            if (int_second) {
                // check if fractional second has wrapped
                if (vrt_packet.fractional_seconds_timestamp > last_fractional_seconds_timestamp ) {
                        last_fractional_seconds_timestamp = vrt_packet.fractional_seconds_timestamp;
                        continue;
                } else {
                    int_second = false;
                    last_update = now;
                    start_time = now;
                    stop_time = start_time + std::chrono::milliseconds(int64_t(1000 * total_time));
                }
            }

            // Assumes ci16_le

            int M = decimation;
            int L = vrt_packet.num_rx_samps;

            for (uint32_t i = 0; i < L/M; i++)
                y[i] = std::complex<float>(0,0);

            for (uint32_t i = 0; i < vrt_packet.num_rx_samps; i++) {
                int16_t re;
                memcpy(&re, (char*)&rx_buffer[vrt_packet.offset+i], 2);
                int16_t img;
                memcpy(&img, (char*)&rx_buffer[vrt_packet.offset+i]+2, 2);

                x[M+i+num_taps] = std::complex<float>(re,img);
            }

            // nomalize phasor and step (for doppler)
            phasor = phasor/std::abs(phasor);
            step = step/std::abs(step);

            if (!channel_mode && freq_offset!=0 && doppler_rate!=0) {
                for (uint32_t i = 0; i < vrt_packet.num_rx_samps; i++) {
                    total_phase -= doppler_rate;
                    step = step * step_dop;
                    phasor = phasor * step;
                    x[M+i+num_taps] *= (std::complex<float>)phasor;
                }
            } else if (!channel_mode && freq_offset!=0 && doppler_rate==0) {
                for (uint32_t i = 0; i < vrt_packet.num_rx_samps; i++) {
                    phasor = phasor * step;
                    x[M+i+num_taps] *= (std::complex<float>)phasor;
                }
            }

            for (uint32_t i = 0; i < M; i++) {

                for (uint32_t k = 0; k < L/M; k++) {

                    tmp_acc[k] = std::complex<float>(0,0);

                    for (uint32_t j = 0; j < taps_per_decimation; j++) {
                        tmp_acc[k] += poly_taps[i][taps_per_decimation-j-1] * x[M-i+j*M+k*M];
                    }
                }

                if (channel_mode && polyfir_channel !=0) {
                    std::complex<float> phase_rotation = std::exp(alpha2*(float)i );
                    for (uint32_t k = 0; k < L/M; k++)
                        y[k] += phase_rotation*tmp_acc[k];
                } else {
                    for (uint32_t k = 0; k < L/M; k++)
                        y[k] += tmp_acc[k];
                }

            }

            // overlap between blocks
            for (uint32_t i = 0; i < num_taps; i++) {
                x[M+i] = x[M+i+L];
            }

            for (uint32_t k = 0; k < L/M; k++) {
                iq_buff[iq_counter] = y[k];
                iq_counter++;
            }

            if (iq_counter==VRT_SAMPLES_PER_PACKET) {

                iq_counter = 0;
                t_samp = 0;

                p.fields.integer_seconds_timestamp = vrt_packet.integer_seconds_timestamp;
                p.fields.fractional_seconds_timestamp = vrt_packet.fractional_seconds_timestamp;
                p.header.packet_count = (uint8_t)frame_count%16;
                frame_count++;

                p.body = (char*)iq_buff;
                p.fields.stream_id = 1;

                zmq_msg_t msg;
                int rc = zmq_msg_init_size (&msg, VRT_DATA_PACKET_SIZE*4);
                int32_t rv = vrt_write_packet(&p, zmq_msg_data(&msg), VRT_DATA_PACKET_SIZE, true);

                zmq_msg_send(&msg, responder, 0);
                zmq_msg_close(&msg);

                num_total_samps += vrt_packet.num_rx_samps;
            }

            if (start_rx and first_frame) {
                std::cout << boost::format(
                                 "# First frame: %u samples, %u full secs, %.09f frac secs")
                                 % vrt_packet.num_rx_samps
                                 % vrt_packet.integer_seconds_timestamp
                                 % ((double)vrt_packet.fractional_seconds_timestamp/1e12)
                          << std::endl;
                first_frame = false;
            }
        }

        if (vrt_packet.extended_context) {
            if (tracking) {
                tracker_process(rx_buffer, sizeof(rx_buffer), &vrt_packet, &tracker_ext_context);
                if (!isnan(tracker_ext_context.doppler_rate)) {
                    doppler_rate = tracker_ext_context.doppler_rate;
                    alpha_dop = complexi*2.0*pi*(double)-doppler_rate;
                    step_dop = std::exp(alpha_dop/pow((double)vrt_context.sample_rate,2));
                    // printf("# Doppler rate update (%s): %f\n", tracker_ext_context.object_name, doppler_rate);
                }
            }
            zmq_msg_t msg;
            zmq_msg_init_size (&msg, len);
            memcpy (zmq_msg_data(&msg), rx_buffer, len);
            zmq_msg_send(&msg, responder, 0);
            zmq_msg_close(&msg);
        }

        if (progress) {
            if (vrt_packet.data)
                last_update_samps += vrt_packet.num_rx_samps;
            const auto time_since_last_update = now - last_update;
            if (time_since_last_update > std::chrono::seconds(1)) {
                const double time_since_last_update_s =
                    std::chrono::duration<double>(time_since_last_update).count();
                const double rate = double(last_update_samps) / time_since_last_update_s;
                std::cout << "\t" << (rate / 1e6) << " Msps, ";

                last_update_samps = 0;
                last_update       = now;

                float sum_i = 0;
                uint32_t clip_i = 0;

                double datatype_max = 32768.;

                for (int i=0; i<vrt_packet.num_rx_samps; i++ ) {
                    auto sample_i = get_abs_val((std::complex<int16_t>)rx_buffer[vrt_packet.offset+i]);
                    sum_i += sample_i;
                    if (sample_i > datatype_max*0.99)
                        clip_i++;
                }
                sum_i = sum_i/vrt_packet.num_rx_samps;
                std::cout << boost::format("%.0f") % (100.0*log2(sum_i)/log2(datatype_max)) << "% I (";
                std::cout << boost::format("%.0f") % ceil(log2(sum_i)+1) << " of ";
                std::cout << (int)ceil(log2(datatype_max)+1) << " bits), ";
                std::cout << "" << boost::format("%.0f") % (100.0*clip_i/vrt_packet.num_rx_samps) << "% I clip, ";
                std::cout << std::endl;

            }
        }
    }

    zmq_close(subscriber);
    zmq_close(responder);
    zmq_ctx_destroy(context);

    return 0;

}
