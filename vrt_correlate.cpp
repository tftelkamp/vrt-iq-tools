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
#include <complex>
#include <csignal>
#include <fstream>
#include <iostream>
#include <thread>

// VRT
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include <string.h>

#include <vrt/vrt_read.h>
#include <vrt/vrt_string.h>
#include <vrt/vrt_types.h>
#include <vrt/vrt_util.h>

#include <complex.h>
#include <fftw3.h>

#include "vrt-tools.h"
#include "dt-extended-context.h"
#include "tracker-extended-context.h"

namespace po = boost::program_options;

#define REAL 0
#define IMAG 1

const double pi = std::acos(-1.0);
const std::complex<double> complexi(0.0, 1.0);

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

    // FFTW
    fftw_plan fft_plan[2];
    fftw_plan ifft_plan;

    uint32_t num_bins;

    // variables to be set by po
    std::string file, type, zmq_address, channel_list;
    size_t num_requested_samples;
    uint32_t bins;
    int gain;
    double total_time;
    uint16_t port;
    uint16_t buffer_depth;
    uint32_t channel;
    uint32_t integrations;
    int hwm;
    float amplitude;
    float bin_size, integration_time = 0.0;

    double delta_range = 0;
    double delta_range_dot = 0;
    double cable_delay = 0;
    double clock_offset = 0;

    double range_u = 0;
    double range_v = 0;

    std::complex<double> clock_phase = 1;

    double current_delta_range = 0;

    const double c = 299792458.0;

    double t0 = 0;
    double t_ephem = 0;

    double current_delay;
    double fractional_delay;
    double current_delay_samples;
    int32_t current_sample_delay;

    std::complex<int16_t> **iq_samples;
    std::complex<double> **signal;
    std::complex<double> **fft_result;
    std::complex<double> *xcorr;
    std::complex<double> *xcorr_integrated;
    std::complex<double> *xcorr_time;
    double *signal_mag;


    // setup the program options
    po::options_description desc("Allowed options");
    // clang-format off

    desc.add_options()
        ("help", "help message")
        ("nsamps", po::value<size_t>(&num_requested_samples)->default_value(0), "total number of samples to receive")
        ("duration", po::value<double>(&total_time)->default_value(0), "total number of seconds to receive")
        // ("progress", "periodically display short-term bandwidth")
        ("channel", po::value<std::string>(&channel_list)->default_value("0,1"), "which VRT channel(s) to use (specify \"0\", \"1\", \"0,1\", etc)")
        ("int-second", "align start of reception to integer second")
        ("num-bins", po::value<uint32_t>(&num_bins)->default_value(1000), "number of bins")
        ("integration-time", po::value<float>(&integration_time)->default_value(1.0), "integration time (seconds)")
        ("amplitude", po::value<float>(&amplitude)->default_value(1), "amplitude correction of second channel")
        ("delta-range", po::value<double>(&delta_range)->default_value(0), "delta range (m)")
        ("delta-range-dot", po::value<double>(&delta_range_dot)->default_value(0), "delate range dot (m/s)")
        ("cable-delay", po::value<double>(&cable_delay)->default_value(0), "delay offset (s)")
        ("clock-offset", po::value<double>(&clock_offset)->default_value(0), "clock offset")
        ("buffer-depth", po::value<uint16_t>(&buffer_depth)->default_value(10), "Correlation buffer depth in VRT frames")
        ("lags", "output lags instead of cross-spectrum")
        ("null", "run without writing to file")
        ("ecsv", "output in ECSV format (Astropy)")
        ("continue", "don't abort on a bad packet")
        ("address", po::value<std::string>(&zmq_address)->default_value("127.0.0.1"), "VRT ZMQ address")
        ("port", po::value<uint16_t>(&port)->default_value(50100), "VRT ZMQ port")
        ("hwm", po::value<int>(&hwm)->default_value(10000), "VRT ZMQ HWM")

    ;
    // clang-format on
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // print the help message
    if (vm.count("help")) {
        std::cout << boost::format("VRT correlator. %s") % desc << std::endl;
        std::cout << std::endl
                  << "This application correlates data from "
                     "a dual channel VRT stream.\n"
                  << std::endl;
        return ~0;
    }

    // bool progress               = vm.count("progress") > 0;
    bool stats                  = vm.count("stats") > 0;
    bool null                   = vm.count("null") > 0;
    bool continue_on_bad_packet = vm.count("continue") > 0;
    bool int_second             = (bool)vm.count("int-second");
    bool ecsv                   = vm.count("ecsv") > 0;
    // bool dt_trace               = vm.count("dt-trace") > 0;
    bool lags                   = vm.count("lags") > 0;

    context_type vrt_context;
    init_context(&vrt_context);

    dt_ext_context_type dt_ext_context;
    tracker_ext_context_type tracker_ext_context;

    packet_type vrt_packet;

    // detect which channels to use
    std::vector<std::string> channel_strings;
    std::vector<size_t> channel_nums;
    boost::split(channel_strings, channel_list, boost::is_any_of("\"',"));
    for (size_t ch = 0; ch < channel_strings.size(); ch++) {
        size_t chan = std::stoi(channel_strings[ch]);
        channel_nums.push_back(std::stoi(channel_strings[ch]));
        vrt_packet.channel_filt |= 1<<std::stoi(channel_strings[ch]);
    }

    if (channel_nums.size() != 2) {
        printf("Only 2 channels supported.\n");
        exit(1);
    }

    // ZMQ

    void *context = zmq_ctx_new();

    void *subscriber = zmq_socket(context, ZMQ_SUB);
    int rc = zmq_setsockopt (subscriber, ZMQ_RCVHWM, &hwm, sizeof hwm);
    std::string connect_string = "tcp://" + zmq_address + ":" + std::to_string(port);
    rc = zmq_connect(subscriber, connect_string.c_str());
    assert(rc == 0);
    zmq_setsockopt(subscriber, ZMQ_SUBSCRIBE, "", 0);

    std::string data_address = "localhost";
    int data_port = 70001;

    void *zmq_client = zmq_socket(context, ZMQ_REQ);
    rc = zmq_setsockopt (zmq_client, ZMQ_RCVHWM, &hwm, sizeof hwm);
    connect_string = "tcp://" + data_address + ":" + std::to_string(data_port);
    rc = zmq_connect(zmq_client, connect_string.c_str());
    assert(rc == 0);

    bool first_frame = true;

    // time keeping
    auto start_time = std::chrono::steady_clock::now();
    auto stop_time =
        start_time + std::chrono::milliseconds(int64_t(1000 * total_time));

    uint32_t buffer[ZMQ_BUFFER_SIZE];
    char fringe_stop_buffer[2000];

    unsigned long long num_total_samps = 0;

    // Track time and samps between updating the BW summary
    auto last_update                     = start_time;
    auto last_req                        = start_time;
    unsigned long long last_update_samps = 0;

    // trigger context update
    last_req -= std::chrono::seconds(1);

    bool start_rx = false;
    uint64_t last_fractional_seconds_timestamp = 0;

    uint32_t signal_pointer = 0;

    bool first_block = true;
    uint32_t integration_counter = 0;


    while (not stop_signal_called
           and (num_requested_samples > num_total_samps or num_requested_samples == 0) ) {

        int fringe_stop_len = zmq_recv(zmq_client, fringe_stop_buffer, 100000, ZMQ_NOBLOCK);

        if (fringe_stop_len > 0) {
            fringe_stop_buffer[fringe_stop_len] = 0;
            double values[10];
            char *pt;
            pt = strtok (fringe_stop_buffer,",");
            uint32_t item = 0;
            while (pt != NULL) {
                values[item] = atof(pt);
                pt = strtok (NULL, ",");
                item++;
            }
            if (values[0] != 0)
                t_ephem = values[0]+values[1]/1e12;
            delta_range = values[2];
            delta_range_dot = values[3];
            cable_delay = values[4];
            clock_offset = values[5];
            range_u = values[6];
            range_v = values[7];

            // printf("# DBG: fringe stop, %f %f %f %f %e %e %f %f\n", values[0], values[1], delta_range, delta_range_dot, cable_delay, clock_offset, range_u, range_v);

            // for next FFT
            current_delay_samples = (double)vrt_context.sample_rate*(current_delta_range/c+cable_delay);
            current_sample_delay = (int32_t)floor(current_delay_samples+0.5);
            fractional_delay = current_delay_samples - (double)current_sample_delay;
        }

        int len = zmq_recv(subscriber, buffer, ZMQ_BUFFER_SIZE, 0);

        const auto now = std::chrono::steady_clock::now();

        if (not vrt_process(buffer, sizeof(buffer), &vrt_context, &vrt_packet)) {
            printf("Not a Vita49 packet?\n");
            continue;
        }

        if (not vrt_packet.context and not vrt_packet.data)
            continue;

        uint32_t ch = 0;
        for(ch = 0; ch<channel_nums.size(); ch++)
            if (vrt_packet.stream_id & (1 << channel_nums[ch]) )
                break;

        uint32_t channel = channel_nums[ch];

        if (vrt_packet.extended_context) {
            dt_process(buffer, sizeof(buffer), &vrt_packet, &dt_ext_context);
        }

        if (not start_rx and vrt_packet.context) {

            if (!ecsv)
                vrt_print_context(&vrt_context);
                start_rx = true;

                integrations = (uint32_t)round((double)integration_time/((double)num_bins/(double)vrt_context.sample_rate));

                if (total_time > 0)
                    num_requested_samples = 2 * total_time * vrt_context.sample_rate; // 2 channels

            iq_samples = (std::complex<int16_t> **) malloc(sizeof(std::complex<int16_t>*)*channel_nums.size());
            signal = (std::complex<double> **) malloc(sizeof(std::complex<double>*)*channel_nums.size());
            fft_result = (std::complex<double> **) malloc(sizeof(std::complex<double>*)*channel_nums.size());

            for (size_t ch=0; ch < channel_nums.size(); ch++)
                signal[ch] = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * num_bins);

            for (size_t ch=0; ch < channel_nums.size(); ch++)
                iq_samples[ch] = (std::complex<int16_t>*) calloc(VRT_SAMPLES_PER_PACKET*(buffer_depth+1), sizeof(std::complex<int16_t>));

            for (size_t ch=0; ch < channel_nums.size(); ch++)
                fft_result[ch] = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * num_bins);

            xcorr_time = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * num_bins);

            xcorr = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * num_bins);
            xcorr_integrated = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * num_bins);
            signal_mag = (double*)calloc(num_bins, sizeof(double));

            for (uint32_t i = 0; i < num_bins; i++) {
                xcorr_integrated[i] = 0;
                signal_mag[i] = 0;
            }

            for (size_t ch=0; ch < channel_nums.size(); ch++)
                fft_plan[ch] = fftw_plan_dft_1d(
                    num_bins,
                    reinterpret_cast<fftw_complex*>(signal[ch]),
                    reinterpret_cast<fftw_complex*>(fft_result[ch]),
                    FFTW_FORWARD,
                    FFTW_ESTIMATE
                );

            ifft_plan = fftw_plan_dft_1d(
                num_bins,
                reinterpret_cast<fftw_complex*>(xcorr_integrated),
                reinterpret_cast<fftw_complex*>(xcorr_time),
                FFTW_BACKWARD,
                FFTW_ESTIMATE
            );

            bin_size = (double)vrt_context.sample_rate/(double)num_bins;

            if (!ecsv) {
                printf("# Correlation parameters:\n");
                printf("#    Mode: %s\n", lags ? "lags" : "cross-spectrum");
                printf("#    Bins: %u\n", num_bins);
                printf("#    Bin size [Hz]: %.2f\n", ((double)vrt_context.sample_rate)/((double)num_bins));
                printf("#    Integrations: %u\n", integrations);
                printf("#    Integration Time [sec]: %.2f\n", (double)integrations*(double)num_bins/(double)vrt_context.sample_rate);
            } else {
                uint32_t first_col = 4;
                printf("# %%ECSV 1.0\n");
                printf("# ---\n");

                uint32_t ch=0;
                while(not (vrt_context.stream_id & (1 << ch) ) )
                    ch++;

                printf("# delimiter: \',\'\n");
                printf("# meta: !!omap\n");
                printf("# - vrt: !!omap\n");
                printf("#   - {stream_id: %u}\n", vrt_context.stream_id);
                printf("#   - {channel: %u}\n", ch);
                printf("#   - {sample_rate: %.1f}\n", (float)vrt_context.sample_rate);
                printf("#   - {frequency: %.1f}\n", (double)vrt_context.rf_freq);
                printf("#   - {bandwidth: %.1f}\n", (float)vrt_context.bandwidth);
                printf("#   - {rx_gain: %.1f}\n", (float)vrt_context.gain);
                printf("#   - {reference: %s}\n", vrt_context.reflock == 1 ? "external" : "internal");
                printf("#   - {time_source: %s}\n", vrt_context.time_cal == 1? "pps" : "internal");
                printf("# - correlation: !!omap\n");
                printf("#   - {mode: %s}\n", lags ? "lags" : "cross-spectrum");
                printf("#   - {bins: %u}\n", num_bins);
                printf("#   - {col_first_bin: %u}\n", first_col);
                printf("#   - {bin_size: %.2f}\n", ((double)vrt_context.sample_rate)/((double)num_bins));
                printf("#   - {integrations: %u}\n", integrations);
                printf("#   - {integration_time: %.2f}\n", (double)integrations*(double)num_bins/(double)vrt_context.sample_rate);

                printf("# datatype:\n");
                printf("# - {name: timestamp, datatype: float64}\n");
                printf("# - {name: u, unit: m, datatype: float64}\n");
                printf("# - {name: v, unit: m, datatype: float64}\n");
                printf("# - {name: w, unit: m, datatype: float64}\n");

                if (lags) {
                    for (int32_t i = 0; i < num_bins; ++i) {
                        printf("# - {name: \'%.4e\', datatype: complex128}\n", (double)(i - (double)num_bins / 2)/(double)vrt_context.sample_rate);
                    }
                } else {
                    for (int32_t i = 0; i < num_bins; ++i) {
                        printf("# - {name: \'%.0f\', datatype: complex128}\n", (double)((double)vrt_context.rf_freq + (i*bin_size - vrt_context.sample_rate/2)));
                    }
                }
                printf("# schema: astropy-2.0\n");
            }
            // Header
            printf("timestamp, u ,v, w");

            if (lags) {
                for (int32_t i = 0; i < num_bins; i++) {
                        printf(", %.4e", (double)(i - (double)num_bins / 2)/(double)vrt_context.sample_rate);
                }
            } else {
                for (int32_t i = 0; i < num_bins; ++i) {
                            printf(", %.0f", (double)((double)vrt_context.rf_freq + (i*bin_size - vrt_context.sample_rate/2)));
                }
            }
            printf("\n");

            current_delay = (double)vrt_context.sample_rate*(delta_range/c+cable_delay);
            current_sample_delay = (int32_t)floor(current_delay+0.5);
            fractional_delay = current_delay - (double)current_sample_delay;

        }

        const auto time_since_last_req = now - last_req;

        if (vrt_packet.data and time_since_last_req > std::chrono::milliseconds(100)) {

            last_req = now;

            char message[512] = "";
            zmq_msg_t msg;
            snprintf(message, 512, "%llu %llu",vrt_packet.integer_seconds_timestamp,vrt_packet.fractional_seconds_timestamp);
            zmq_msg_init_size(&msg, strlen(message));
            memcpy(zmq_msg_data(&msg), message, strlen(message));
            zmq_msg_send(&msg, zmq_client, 0);
            zmq_msg_close(&msg);
        }

        if (t0 == 0 and vrt_packet.data) {
            uint64_t seconds = vrt_packet.integer_seconds_timestamp;
            uint64_t frac_seconds = vrt_packet.fractional_seconds_timestamp;
            t0 = (double)seconds + (double)frac_seconds/1e12;
            t_ephem = t0;
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
                }
            }

            if (first_block and (ch!=0) ) {
                continue;
            } else {
                first_block = false;
            }

            for (int32_t b = 0; b < buffer_depth; b++) { // TODO: only copy the channel that is actually shifted
                memcpy((char*)&iq_samples[ch][b*VRT_SAMPLES_PER_PACKET], (char*)&iq_samples[ch][(b+1)*VRT_SAMPLES_PER_PACKET], VRT_SAMPLES_PER_PACKET*sizeof(std::complex<int16_t>)); // sizeof
            }

            for (int32_t i = 0; i < vrt_packet.num_rx_samps; i++) {

                int16_t re;
                memcpy(&re, (char*)&buffer[vrt_packet.offset+i], 2);
                int16_t img;
                memcpy(&img, (char*)&buffer[vrt_packet.offset+i]+2, 2);

                iq_samples[ch][i+VRT_SAMPLES_PER_PACKET*buffer_depth] = std::complex<int16_t>(re,img);

            }

            if (ch==1) {
                // both channels received (we assume they are in order)

                for (int32_t k = 0; k < vrt_packet.num_rx_samps; k++) {

                    int32_t ch0_shift = 0;
                    int32_t ch1_shift = 0;

                    if (current_sample_delay < 0) {
                        // shift ch1
                        ch1_shift = abs(current_sample_delay);
                    } else {
                        // chift ch0
                        ch0_shift = abs(current_sample_delay);
                    }

                    signal[0][signal_pointer] = std::complex<double>(
                        (double)iq_samples[0][buffer_depth*VRT_SAMPLES_PER_PACKET-ch0_shift+k].real() / 32768.0,
                        (double)iq_samples[0][buffer_depth*VRT_SAMPLES_PER_PACKET-ch0_shift+k].imag() / 32768.0
                    );

                    signal[1][signal_pointer] = std::complex<double>(
                        (double)iq_samples[1][buffer_depth*VRT_SAMPLES_PER_PACKET-ch1_shift+k].real() / 32768.0,
                        (double)iq_samples[1][buffer_depth*VRT_SAMPLES_PER_PACKET-ch1_shift+k].imag() / 32768.0
                    );

                    signal_pointer++;

                    if (signal_pointer == num_bins) {

                        int64_t seconds = vrt_packet.integer_seconds_timestamp;
                        int64_t frac_seconds = vrt_packet.fractional_seconds_timestamp;
                        frac_seconds += (k-num_bins/2)*1e12/vrt_context.sample_rate;
                        if (frac_seconds > 1e12) {
                            frac_seconds -= 1e12;
                            seconds++;
                        } else if (frac_seconds < 0) {
                            frac_seconds += 1e12;
                            seconds--;
                        }

                        double t = (double)seconds + (double)frac_seconds/1e12;

                        signal_pointer = 0;

                        fftw_execute(fft_plan[0]);
                        fftw_execute(fft_plan[1]);

                        integration_counter++;

                        current_delta_range = delta_range + delta_range_dot * (t-t_ephem);
                        current_delay = current_delta_range/c + cable_delay;

                        std::complex<double> phase_correction = std::exp(complexi*-2.0*pi*(double)vrt_context.rf_freq*current_delay);

                        double df = clock_offset * vrt_context.rf_freq;

                        clock_phase =  std::exp(-2.0 * complexi * pi * df * (t-t0));

                        // correlate and integrate
                        for (int32_t i = 0; i < num_bins; i++) {
                            std::complex<double> frac_corr = std::exp(-1.0 * complexi * pi * fractional_delay *
                                (double)(2.0 * (((int)i + (int)num_bins / 2) % (int)num_bins) / (double)num_bins - 1.0)
                            );
                            xcorr[i] = fft_result[0][i] * conj(fft_result[1][i]);
                            xcorr_integrated[i] +=  phase_correction * clock_phase * frac_corr * xcorr[i];
                            signal_mag[i] += std::abs(fft_result[0][i]) * std::abs(fft_result[1][i]);
                        }

                        if (integration_counter == integrations) {

                            for (int32_t i = 0; i < num_bins; i++) {
                                if (signal_mag[i]>0)
                                    xcorr_integrated[i] = xcorr_integrated[i] / signal_mag[i];
                            }

                            printf("%llu.%09lli", seconds, (int64_t)(frac_seconds/1e3));

                            printf(", %.12e, %.12e, %.12e", range_u, range_v, delta_range);

                            if (lags) {
                                // inverse FFT
                                fftw_execute(ifft_plan);

                                for (uint32_t i = 0; i < num_bins; i++) {
                                    printf(", (%.6e%s%.6ej)", xcorr_time[i].real(), (xcorr_time[i].imag() > 0) ? "+" : "-", abs(xcorr_time[i].imag()) );
                                }
                            } else {
                                for (uint32_t i = 0; i < num_bins; i++) {
                                    printf(", (%.6e%s%.6ej)", xcorr_integrated[i].real(), (xcorr_integrated[i].imag() > 0) ? "+" : "-", abs(xcorr_integrated[i].imag()) );
                                }
                            }

                            printf("\n");
                            fflush(stdout);

                            integration_counter = 0;

                            for (uint32_t i = 0; i < num_bins; i++) {
                                xcorr_integrated[i] = 0;
                                signal_mag[i] = 0;
                            }
                        }

                        // set for next FFT
                        current_delay_samples = (double)vrt_context.sample_rate*(current_delta_range/c+cable_delay);
                        current_sample_delay = (int32_t)floor(current_delay_samples+0.5);
                        fractional_delay = current_delay_samples - (double)current_sample_delay;
                    }
                }
            }

            // end
            num_total_samps += vrt_packet.num_rx_samps;
        }
    }

    zmq_close(subscriber);
    zmq_ctx_destroy(context);

    return 0;
}
