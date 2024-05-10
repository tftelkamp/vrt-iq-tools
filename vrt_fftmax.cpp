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
#include <fftw3.h>

#include "vrt-tools.h"

namespace po = boost::program_options;

#define REAL 0
#define IMAG 1

#define SCALE_MAX 32768

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
    fftw_complex *signal, *result;
    fftw_plan plan;
    uint32_t num_points = 0;
    uint32_t fft_len = 1;

    int32_t min_bin, max_bin;

    // variables to be set by po
    std::string file, type, zmq_address;
    uint16_t port;
    uint32_t channel;
    int hwm;
    size_t num_requested_samples;
    double total_time, min_offset, max_offset;

    // setup the program options
    po::options_description desc("Allowed options");
    // clang-format off

    desc.add_options()
        ("help", "help message")
        ("nsamps", po::value<size_t>(&num_requested_samples)->default_value(0), "total number of samples to receive")
        ("duration", po::value<double>(&total_time)->default_value(0), "total number of seconds to receive")
        ("min-offset", po::value<double>(&min_offset), "min. freq. offset to track")
        ("max-offset", po::value<double>(&max_offset), "max. freq. offset to track")
        ("fft-duration", po::value<uint32_t>(&fft_len), "number of seconds to integrate")
        ("channel", po::value<uint32_t>(&channel)->default_value(0), "VRT channel")
        ("progress", "periodically display short-term bandwidth")
        // ("stats", "show average bandwidth on exit")
        ("int-second", "align start of reception to integer second")
        ("null", "run without writing to file")
        ("continue", "don't abort on a bad packet")
        ("ignore-dc", "Ignore  DC bin")
        ("address", po::value<std::string>(&zmq_address)->default_value("localhost"), "VRT ZMQ address")
        ("port", po::value<uint16_t>(&port)->default_value(50100), "VRT ZMQ port")
        ("hwm", po::value<int>(&hwm)->default_value(10000), "VRT ZMQ HWM")
    ;
    // clang-format on
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // print the help message
    if (vm.count("help")) {
        std::cout << boost::format("VRT samples to fftmax. %s") % desc << std::endl;
        std::cout << std::endl
                  << "This application streams data from a VRT stream "
                     "to fftmax.\n"
                  << std::endl;
        return ~0;
    }

    bool progress               = vm.count("progress") > 0;
    bool null                   = vm.count("null") > 0;
    bool continue_on_bad_packet = vm.count("continue") > 0;
    bool int_second             = (bool)vm.count("int-second");
    bool ignore_dc              = (bool)vm.count("ignore-dc");

    context_type vrt_context;
    init_context(&vrt_context);

    packet_type vrt_packet;

    vrt_packet.channel_filt = 1<<channel;

    // ZMQ
    void *context = zmq_ctx_new();
    void *subscriber = zmq_socket(context, ZMQ_SUB);
    int rc = zmq_setsockopt (subscriber, ZMQ_RCVHWM, &hwm, sizeof hwm);
    std::string connect_string = "tcp://" + zmq_address + ":" + std::to_string(port);
    rc = zmq_connect(subscriber, connect_string.c_str());
    assert(rc == 0);
    zmq_setsockopt(subscriber, ZMQ_SUBSCRIBE, "", 0);

    // time keeping
    auto start_time = std::chrono::steady_clock::now();
    auto stop_time = start_time + std::chrono::milliseconds(int64_t(1000 * total_time));

    uint32_t buffer[ZMQ_BUFFER_SIZE];

    unsigned long long num_total_samps = 0;

    // Track time and samps between updating the BW summary
    auto last_update                     = start_time;
    unsigned long long last_update_samps = 0;

    bool first_frame = true;
    bool start_rx = false;
    uint64_t last_fractional_seconds_timestamp = 0;

    uint32_t signal_pointer = 0;

    while (not stop_signal_called
           and (num_requested_samples > num_total_samps or num_requested_samples == 0)
           and (total_time == 0.0 or std::chrono::steady_clock::now() <= stop_time)) {

        int len = zmq_recv(subscriber, buffer, ZMQ_BUFFER_SIZE, 0);

        const auto now = std::chrono::steady_clock::now();

        if (not vrt_process(buffer, sizeof(buffer), &vrt_context, &vrt_packet)) {
            printf("Not a Vita49 packet?\n");
            continue;
        }

        // if (vrt_context.context_changed) {
        //     printf("Context changed\n");
        //     break;
        // }

        if (not start_rx and vrt_packet.context) {
            vrt_print_context(&vrt_context);
            start_rx = true;
            num_points = fft_len*vrt_context.sample_rate;

            min_bin = 0;
            max_bin = num_points;

            if (vm.count("min-offset")) {
                min_bin = min_offset+num_points/2;
                min_bin = min_bin < 0 ? 0 : min_bin;
                min_bin = min_bin > num_points ? num_points : min_bin;
            }

            if (vm.count("max-offset")) {
                max_bin = max_offset+num_points/2;
                max_bin = max_bin < 0 ? 0 : max_bin;
                max_bin = max_bin > num_points ? num_points : max_bin;
            }

            signal = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_points);
            result = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_points);
            plan = fftw_plan_dft_1d(num_points, signal, result, FFTW_FORWARD, FFTW_ESTIMATE);
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

            int mult = 1;
            for (uint32_t i = 0; i < vrt_packet.num_rx_samps; i++) {
                int16_t re;
                memcpy(&re, (char*)&buffer[vrt_packet.offset+i], 2);
                int16_t img;
                memcpy(&img, (char*)&buffer[vrt_packet.offset+i]+2, 2);
                signal[signal_pointer][REAL] = mult*re;
                signal[signal_pointer][IMAG] = mult*img;
                mult *= -1;

                signal_pointer++;

                if (signal_pointer >= num_points) {

                    signal_pointer = 0;

                    fftw_execute(plan);

                    double max = 0;
                    int32_t max_i = -1;

                    uint32_t dc = num_points/2;

                    for (uint32_t i = 0; i < num_points; ++i) {
                        double mag = sqrt(result[i][REAL] * result[i][REAL] +
                                  result[i][IMAG] * result[i][IMAG]);
                        if ( (mag > max) and (i >= min_bin) and (i <= max_bin) and not (ignore_dc && i==dc)) {
                            max = mag;
                            max_i = i;
                        }
                    }

                    uint64_t seconds = vrt_packet.integer_seconds_timestamp;
                    uint64_t frac_seconds = vrt_packet.fractional_seconds_timestamp;
                    frac_seconds += (i+1)*1e12/vrt_context.sample_rate;
                    if (frac_seconds > 1e12) {
                        frac_seconds -= 1e12;
                        seconds++;
                    }

                    double peak_hz = vrt_context.rf_freq + (double)max_i/(double)fft_len - vrt_context.sample_rate/2;
                    printf("%lu.%09li, %.2f, %.3f\n", static_cast<unsigned long>(seconds), static_cast<long>(frac_seconds/1e3), peak_hz, 20*log10(max/(double)num_points));
                    fflush(stdout);
                }
            }


            num_total_samps += vrt_packet.num_rx_samps;

            if (start_rx and first_frame) {
                std::cout << boost::format(
                                 "# First frame: %u samples, %u full secs, %.09f frac secs")
                                 % vrt_packet.num_rx_samps
                                 % vrt_packet.integer_seconds_timestamp
                                 % ((double)vrt_packet.fractional_seconds_timestamp/1e12)
                          << std::endl;
                first_frame = false;
                // Header
                printf("timestamp, frequency, power\n");
            }
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
                    auto sample_i = get_abs_val((std::complex<int16_t>)buffer[vrt_packet.offset+i]);
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
    zmq_ctx_destroy(context);

    return 0;

}
