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

#include <complex.h>
// #include <complex>

// CUDA FFT
#include <cufft.h>

#include "vrt-tools.h"
#include "tracker-extended-context.h"

#define REAL 0
#define IMAG 1

const double pi = std::acos(-1.0);

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

// CUDA
__global__
void pfb(cuFloatComplex* __restrict__ data, cuFloatComplex* __restrict__ shift_reg, float* __restrict__ hh2, int bins, int samples, int pfb_filter_taps, int osr)
{
    int k = blockIdx.x*blockDim.x + threadIdx.x;
    int sample = blockIdx.y;
    int t, b;

    uint32_t index = samples - ( sample + 1) * ( bins / osr );
    uint32_t shift = ( ( sample + 1 ) * ( bins / osr ) ) % bins;

    b = ( bins + k - shift ) % bins;

    cuFloatComplex tmp = make_cuFloatComplex(0, 0);

    for (t=0; t < pfb_filter_taps; t++) {
        tmp.x += shift_reg[ index + ((k+t*bins)%(bins*pfb_filter_taps)) ].x * hh2[ k + t * bins ];
        tmp.y += shift_reg[ index + ((k+t*bins)%(bins*pfb_filter_taps)) ].y * hh2[ k + t * bins ];
    }

    data[ sample * bins + b ] = tmp;
}

int main(int argc, char* argv[])
{

    // variables to be set by po
    std::string file, type, zmq_address;
    uint16_t pub_instance, instance, main_port, port, pub_port;
    uint32_t channel;
    int hwm, io_threads;
    float freq_offset, rate, doppler_rate;
    double frequency;
    size_t num_requested_samples;
    double total_time;
    float channel_bw;
    uint32_t gpu_frames;

    uint32_t decimation, osr, taps_per_decimation, num_taps;

    double *taps;
    float *hh2;

    // C
    std::complex<int16_t> **iq_buff;
    std::complex<float> *shift_reg;
    std::complex<float> *ifft_out;

    // CUDA
    cuFloatComplex *cuda_shift_reg;
    cuFloatComplex *cuda_poly_filter_out;

    cufftHandle cuda_plan;

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
        ("decimation", po::value<uint32_t>(&decimation)->default_value(2), "decimation factor")
        ("osr", po::value<uint32_t>(&osr)->default_value(1), "channel oversampling rate")
        ("taps-per-decimation", po::value<uint32_t>(&taps_per_decimation)->default_value(20), "taps per decimation")
        ("rate", po::value<float>(&rate)->default_value(0), "channel rate")
        ("channel-bw", po::value<float>(&channel_bw)->default_value(0.97), "channel bandwidth as fraction of rate")
        ("gpu-frames", po::value<uint32_t>(&gpu_frames)->default_value(1), "output frames per GPU call")
        ("address", po::value<std::string>(&zmq_address)->default_value("localhost"), "VRT ZMQ address")
        ("zmq-split", "use a ZeroMQ stream per VRT channel, increasing port number for additional streams")
        ("instance", po::value<uint16_t>(&instance)->default_value(0), "VRT ZMQ instance")
        ("port", po::value<uint16_t>(&port), "VRT ZMQ port")
        ("pub-zmq-split", "create a ZeroMQ stream per VRT channel, increasing port number for additional streams")
        ("pub-port", po::value<uint16_t>(&pub_port), "VRT ZMQ PUB port")
        ("pub-instance", po::value<uint16_t>(&pub_instance)->default_value(1), "VRT ZMQ instance")
        ("io-threads", po::value<int>(&io_threads)->default_value(1), "ZMQ IO threads")
        ("hwm", po::value<int>(&hwm)->default_value(10000), "VRT ZMQ HWM")
    ;
    // clang-format on
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // print the help message
    if (vm.count("help")) {
        std::cout << boost::format("VRT pfb channelizer. %s") % desc << std::endl;
        std::cout << std::endl
                  << "This application channelizes a VRT stream.\n"
                  << std::endl;
        return ~0;
    }

    bool progress               = vm.count("progress") > 0;
    bool stats                  = vm.count("stats") > 0;
    bool null                   = vm.count("null") > 0;
    bool continue_on_bad_packet = vm.count("continue") > 0;
    bool int_second             = (bool)vm.count("int-second");
    bool zmq_split              = vm.count("zmq-split") > 0;
    bool pub_zmq_split          = vm.count("pub-zmq-split") > 0;

    // CUDA
    cudaFree(0);

    context_type vrt_context;
    init_context(&vrt_context);
    tracker_ext_context_type tracker_ext_context;

    packet_type vrt_packet;

    if (rate > 0 && pub_zmq_split) {
        printf("specify --decimation instead of --rate when using --pub-zmq-split.\n");
        exit(1);
    }

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

    void *zmq_server[100]; // TODO use MAX_CHANNELS
    void *context = zmq_ctx_new();
    void *responder;
    int rc;

    zmq_ctx_set (context, ZMQ_IO_THREADS, io_threads);

    // ZMQ
    void *subscriber = zmq_socket(context, ZMQ_SUB);
    rc = zmq_setsockopt (subscriber, ZMQ_RCVHWM, &hwm, sizeof hwm);
    std::string connect_string = "tcp://" + zmq_address + ":" + std::to_string(main_port);
    rc = zmq_connect(subscriber, connect_string.c_str());
    assert(rc == 0);
    zmq_setsockopt(subscriber, ZMQ_SUBSCRIBE, "", 0);


    if (pub_zmq_split) {
        for (size_t ch = 0; ch < decimation; ch++) {
            responder = zmq_socket(context, ZMQ_PUB);
            rc = zmq_setsockopt (responder, ZMQ_SNDHWM, &hwm, sizeof hwm);
            assert(rc == 0);
            std::string connect_string = "tcp://*:" + std::to_string(pub_port+ch);
            rc = zmq_bind(responder, connect_string.c_str());
            assert (rc == 0);
            zmq_server[ch] = responder;
        }
    } else {
        responder = zmq_socket(context, ZMQ_PUB);
        rc = zmq_setsockopt (responder, ZMQ_SNDHWM, &hwm, sizeof hwm);
        assert(rc == 0);
        connect_string = "tcp://*:" + std::to_string(pub_port);
        rc = zmq_bind(responder, connect_string.c_str());
        assert (rc == 0);
        zmq_server[0] = responder;
    }

    // time keeping
    auto start_time = std::chrono::steady_clock::now();
    auto stop_time = start_time + std::chrono::milliseconds(int64_t(1000 * total_time));

    uint32_t rx_buffer[ZMQ_BUFFER_SIZE];
    uint32_t tx_buffer[ZMQ_BUFFER_SIZE];

    uint64_t num_total_samps = 0;

    // Track time and samps between updating the BW summary
    auto last_update                     = start_time;
    uint64_t last_update_samps = 0;

    bool first_frame = true;
    uint64_t last_fractional_seconds_timestamp = 0;

    // set to true to process data before context
    bool start_rx = false;

    /* VRT init */
    struct vrt_packet p;
    vrt_init_packet(&p);
    vrt_init_data_packet(&p);
    p.fields.stream_id = 1;

    uint32_t iq_counter = 0;
    uint32_t frame_count = 0;

    bool first_context = true;

    uint64_t next_integer_seconds_timestamp = 0;
    uint64_t next_fractional_seconds_timestamp = 0;

    uint32_t buffer_frames = 0;
    uint32_t frame_counter = 0;

    while (not stop_signal_called
           and (num_requested_samples > num_total_samps or num_requested_samples == 0)
           and (total_time == 0.0 or std::chrono::steady_clock::now() <= stop_time)) {

        int len = zmq_recv(subscriber, rx_buffer, ZMQ_BUFFER_SIZE, 0);

        const auto now = std::chrono::steady_clock::now();

        if (not vrt_process(rx_buffer, sizeof(rx_buffer), &vrt_context, &vrt_packet)) {
            printf("Not a Vita49 packet?\n");
            continue;
        }

        if (not start_rx and vrt_packet.context) {
            vrt_print_context(&vrt_context);
            start_rx = true;

            if (frequency > 0) {
                freq_offset = frequency-(double)vrt_context.rf_freq;
            }

            if (rate > 0) {
                decimation = vrt_context.sample_rate/rate;
            } else {
                rate = vrt_context.sample_rate/decimation;
            }

            // check for valid bandwidth
            if ( fmod((float)vrt_context.sample_rate,rate) != 0) {
                printf("channel rate needs to be a divisor of the sample rate (%u).\n", vrt_context.sample_rate);
                exit(1);
            }

            // check for valid frequency offset
            if ( fabs(freq_offset) > vrt_context.sample_rate/2) {
                printf("Selected frequency outside of stream\n");
                exit(1);
            }

            // // check for valid decimation
            // if (VRT_SAMPLES_PER_PACKET % decimation != 0) {
            //     printf("decimation needs to be a divisor of %u.\n", VRT_SAMPLES_PER_PACKET);
            //     exit(1);
            // }

            // // check for valid decimation
            // if (VRT_SAMPLES_PER_PACKET % (osr*VRT_SAMPLES_PER_PACKET/decimation) != 0) {
            //     printf("osr*10000/decimation needs to be a divisor of %u.\n", VRT_SAMPLES_PER_PACKET);
            //     exit(1);
            // }
            
            // check for valid decimation
            if ((uint64_t)vrt_context.sample_rate % decimation != 0) {
                printf("decimation needs to be a divisor of the sample rate (%u).\n", vrt_context.sample_rate);
                exit(1);
            }

            // check for valid decimation
            if (decimation > 32 && !pub_zmq_split) {
                printf("maximum decimation is 32 when not using --pub-zmq-split.\n");
                exit(1);
            }

            buffer_frames = gpu_frames * decimation / osr;

            // create FIR filter
            uint32_t fir_order = taps_per_decimation*decimation-1;
            num_taps = fir_order+1;
            taps = (double*)malloc(sizeof(double)*num_taps);
            // hh2 = (float*)malloc(sizeof(float)*num_taps);

            cudaMallocManaged((void**) &hh2, num_taps * sizeof(float));
            cudaStreamAttachMemAsync(NULL, hh2, 0, cudaMemAttachHost);
            cudaStreamSynchronize(NULL);

            double K = channel_bw*(fir_order/decimation);

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

            double norm_sum = 0;

            for (int i=0;i<fir_order;i++) {
                int j = -(i - fir_order/2);
                double blackman_window = a0 - a1*cos(2*pi*(double)i/((double)fir_order-1)) +
                                            a2*cos(4*pi*(double)i/((double)fir_order-1)) +
                                            a3*cos(6*pi*(double)i/((double)fir_order-1));
                if (j==0) {
                    taps[i] = blackman_window*((double)K/(double)fir_order);
                } else {
                    taps[i] = blackman_window*(1.0/(double)fir_order)*sin(pi*(double)j*(double)K/(double)fir_order)/(pi*(double)j/(double)fir_order);
                }

                norm_sum += taps[i];
            }

            taps[fir_order] = 0;

            for (uint32_t i=0; i<num_taps; i++) {
                hh2[i] = taps[i]/norm_sum;
            }

            cudaStreamAttachMemAsync(NULL, hh2, 0, cudaMemAttachGlobal);

            iq_buff = (std::complex<int16_t> **)malloc(sizeof(std::complex<int16_t> *) * decimation);

            for (size_t dec=0; dec < decimation; dec++)
                iq_buff[dec] = (std::complex<int16_t> *)malloc(sizeof(std::complex<int16_t>) * gpu_frames * VRT_SAMPLES_PER_PACKET);

            int rank = 1;
            int n[] = {(int)decimation};
            int howmany = buffer_frames*osr*VRT_SAMPLES_PER_PACKET/decimation;
            int idist = decimation, odist = decimation;
            int istride = 1, ostride = 1; // array is contiguous in memory
            int *inembed = n, *onembed = n;

            cufftPlanMany(&cuda_plan, rank, n, inembed, istride, idist,
                            onembed, ostride, odist, CUFFT_C2C, howmany);

            // CUDA
            cudaMalloc((void**) &cuda_shift_reg, (buffer_frames+1)*VRT_SAMPLES_PER_PACKET * sizeof(cuFloatComplex));
            cudaMalloc((void**) &cuda_poly_filter_out, sizeof(cufftComplex)*buffer_frames*VRT_SAMPLES_PER_PACKET*osr);
            cudaMallocHost((void **)&shift_reg,  (buffer_frames+1)*VRT_SAMPLES_PER_PACKET * sizeof(std::complex<float>));
            cudaMallocHost((void **)&ifft_out,  buffer_frames*VRT_SAMPLES_PER_PACKET*osr*sizeof(std::complex<float>));
            cudaStreamSynchronize(NULL);

        }

        if (start_rx and vrt_packet.context) {
            // construct new context
            struct vrt_packet pc;
            vrt_init_packet(&pc);
            vrt_init_context_packet(&pc);

            pc.fields.integer_seconds_timestamp = vrt_context.integer_seconds_timestamp;
            pc.fields.fractional_seconds_timestamp = vrt_context.fractional_seconds_timestamp;

            if (doppler_rate!=0 || first_context) {
                pc.if_context.context_field_change_indicator = true;
                first_context = false;
            }
            else
                pc.if_context.context_field_change_indicator = false;

            pc.if_context.bandwidth = vrt_context.bandwidth;
            pc.if_context.sample_rate = osr*vrt_context.sample_rate/decimation;
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

            for (int dec = 0; dec < decimation; dec++) {
     
                if (pub_zmq_split)
                    pc.fields.stream_id = 1;
                else 
                    pc.fields.stream_id = 1<<dec;

                pc.if_context.rf_reference_frequency = (double)vrt_context.rf_freq + ((double)dec-(double)decimation/2)*(double)vrt_context.sample_rate/(double)decimation;

                int32_t rv = vrt_write_packet(&pc, tx_buffer, VRT_DATA_PACKET_SIZE, true);
                if (rv < 0) {
                    fprintf(stderr, "Failed to write packet: %s\n", vrt_string_error(rv));
                }

                // ZMQ
                if (pub_zmq_split)
                    zmq_send (zmq_server[dec], tx_buffer, rv*4, 0);
                else
                    zmq_send (zmq_server[0], tx_buffer, rv*4, 0);
            }

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
            for (uint32_t i = 0; i < vrt_packet.num_rx_samps; i++) {
                int16_t re;
                memcpy(&re, (char*)&rx_buffer[vrt_packet.offset+i], 2);
                int16_t img;
                memcpy(&img, (char*)&rx_buffer[vrt_packet.offset+i]+2, 2);

                shift_reg[(buffer_frames-frame_counter)*VRT_SAMPLES_PER_PACKET+vrt_packet.num_rx_samps-i-1] = std::complex<float>(re,img);
            }

            frame_counter += 1;

            if (frame_counter == buffer_frames) {

                int samples_per_channel_out = buffer_frames*osr*VRT_SAMPLES_PER_PACKET/decimation;

                int block_samples_in = buffer_frames*VRT_SAMPLES_PER_PACKET;

                int bins = decimation;

                int pfb_filter_taps = taps_per_decimation;

                cudaMemcpy(
                    cuda_shift_reg,
                    shift_reg,
                    (buffer_frames*VRT_SAMPLES_PER_PACKET + bins * pfb_filter_taps / 2) * sizeof(std::complex<float>),
                    cudaMemcpyHostToDevice
                );

                // PFB
                dim3 dimBlock(bins,samples_per_channel_out);

                pfb<<<dimBlock,1>>>(
                    cuda_poly_filter_out,
                    cuda_shift_reg,
                    hh2,
                    bins,
                    block_samples_in,
                    pfb_filter_taps,
                    osr
                );

                // iFFT
                if (cufftExecC2C(cuda_plan, cuda_poly_filter_out, cuda_poly_filter_out, CUFFT_INVERSE) != CUFFT_SUCCESS){
                    printf("CUFFT error: ExecC2C Forward failed\n");
                    exit(1);
                }

                cudaDeviceSynchronize();

                cudaMemcpy(
                    ifft_out,
                    cuda_poly_filter_out,
                    sizeof(cufftComplex)*VRT_SAMPLES_PER_PACKET*osr*buffer_frames,
                    cudaMemcpyDeviceToHost
                );

                // reorder
                for (int loops = 0; loops < samples_per_channel_out; loops++) {
                    for (int k = 0; k < bins; k++) {
                        int offset = k + ((k < bins/2) * 2 - 1) * bins/2;

                        iq_buff[k][iq_counter+loops] = ifft_out[loops*bins+offset];

                    }
                }

                int noverlap = bins * pfb_filter_taps / 2;
                memcpy(&shift_reg[block_samples_in], shift_reg, noverlap * sizeof(std::complex<float>));

                iq_counter += samples_per_channel_out;

                if (next_integer_seconds_timestamp == 0) {
                    next_integer_seconds_timestamp = vrt_packet.integer_seconds_timestamp;
                    next_fractional_seconds_timestamp = vrt_packet.fractional_seconds_timestamp;
                }

                frame_counter = 0;

                while (iq_counter>=VRT_SAMPLES_PER_PACKET) {

                    iq_counter -= VRT_SAMPLES_PER_PACKET;

                    p.fields.integer_seconds_timestamp = next_integer_seconds_timestamp;
                    p.fields.fractional_seconds_timestamp = next_fractional_seconds_timestamp;
                    p.header.packet_count = (uint8_t)frame_count%16;
                    frame_count++;

                    next_integer_seconds_timestamp = 0;
                    next_fractional_seconds_timestamp = 0;

                    for (int dec = 0; dec < decimation; dec++) {
                        p.body = (char*)iq_buff[dec];
                        if (pub_zmq_split)
                            p.fields.stream_id = 1;
                        else 
                            p.fields.stream_id = 1<<dec;

                        zmq_msg_t msg;
                        int rc = zmq_msg_init_size (&msg, VRT_DATA_PACKET_SIZE*4);
                        int32_t rv = vrt_write_packet(&p, zmq_msg_data(&msg), VRT_DATA_PACKET_SIZE, true);

                        if (pub_zmq_split)
                            zmq_msg_send(&msg, zmq_server[dec], 0);
                        else
                            zmq_msg_send(&msg, zmq_server[0], 0);
                        zmq_msg_close(&msg);
                    }

                    num_total_samps += vrt_packet.num_rx_samps;

                }
            
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
            zmq_msg_t msg;
            zmq_msg_init_size (&msg, len);
            memcpy (zmq_msg_data(&msg), rx_buffer, len);
            zmq_msg_send(&msg, zmq_server[0], 0); // TODO send to all channels
            zmq_msg_close(&msg);
        }

        if (progress && vrt_packet.data)
            show_progress_stats(
                now,
                &last_update,
                &last_update_samps,
                &rx_buffer[vrt_packet.offset],
                vrt_packet.num_rx_samps, channel
            );           
    }

    zmq_close(subscriber);
    if (pub_zmq_split)
        for (size_t ch = 0; ch < decimation; ch++)
            zmq_close(zmq_server[ch]);
    else
        zmq_close(zmq_server[0]);
    zmq_ctx_destroy(context);

    return 0;

}
