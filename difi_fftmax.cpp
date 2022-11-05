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

#include "difi-tools.h"

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

    // variables to be set by po
    std::string file, type, zmq_address;
    uint16_t port;
    int hwm;
    size_t num_requested_samples;
    double total_time;

    // setup the program options
    po::options_description desc("Allowed options");
    // clang-format off

    desc.add_options()
        ("help", "help message")
        // ("file", po::value<std::string>(&file)->default_value("usrp_samples.dat"), "name of the file to write binary samples to")
        // ("type", po::value<std::string>(&type)->default_value("short"), "sample type: double, float, or short")
        ("nsamps", po::value<size_t>(&num_requested_samples)->default_value(0), "total number of samples to receive")
        ("duration", po::value<double>(&total_time)->default_value(0), "total number of seconds to receive")
        ("progress", "periodically display short-term bandwidth")
        // ("stats", "show average bandwidth on exit")
        ("int-second", "align start of reception to integer second")
        ("null", "run without writing to file")
        ("continue", "don't abort on a bad packet")
        ("ignore-dc", "Ignore 10 perc. of bins around DC")
        ("address", po::value<std::string>(&zmq_address)->default_value("localhost"), "DIFI ZMQ address")
        ("port", po::value<uint16_t>(&port)->default_value(50100), "DIFI ZMQ port")
        ("hwm", po::value<int>(&hwm)->default_value(10000), "DIFI ZMQ HWM")
    ;
    // clang-format on
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // print the help message
    if (vm.count("help")) {
        std::cout << boost::format("DIFI samples to fftmax %s") % desc << std::endl;
        std::cout << std::endl
                  << "This application streams data from a DIFI stream "
                     "to fftmax.\n"
                  << std::endl;
        return ~0;
    }

    bool progress               = vm.count("progress") > 0;
    bool stats                  = vm.count("stats") > 0;
    bool null                   = vm.count("null") > 0;
    bool continue_on_bad_packet = vm.count("continue") > 0;
    bool int_second             = (bool)vm.count("int-second");
    bool ignore_dc              = (bool)vm.count("ignore-dc");

    context_type difi_context;
    init_context(&difi_context);

    difi_packet_type difi_packet;

    // std::vector<std::shared_ptr<std::ofstream>> outfiles;
    std::vector<size_t> channel_nums = {0}; // single channel (0)

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

        if (not difi_process(buffer, sizeof(buffer), &difi_context, &difi_packet)) {
            printf("Not a Vita49 packet?\n");
            continue;
        }

        if (not start_rx and difi_packet.context) {
            difi_print_context(&difi_context);
            start_rx = true;
            num_points = difi_context.sample_rate;
            signal = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_points);
            result = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_points);
            plan = fftw_plan_dft_1d(num_points, signal, result, FFTW_FORWARD, FFTW_ESTIMATE);
        }
        
        if (start_rx and difi_packet.data) {

            if (difi_packet.lost_frame)
               if (not continue_on_bad_packet)
                    break;

            if (int_second) {
                // check if fractional second has wrapped
                if (difi_packet.fractional_seconds_timestamp > last_fractional_seconds_timestamp ) {
                        last_fractional_seconds_timestamp = difi_packet.fractional_seconds_timestamp;
                        continue;
                } else {
                    int_second = false;
                    last_update = now; 
                    start_time = now;
                    stop_time = start_time + std::chrono::milliseconds(int64_t(1000 * total_time));
                }
            }

            int mult = 1;
            for (uint32_t i = 0; i < difi_packet.num_rx_samps; i++) {
                signal_pointer++;
                if (signal_pointer >= num_points)  
                    break;
                int16_t re;
                memcpy(&re, (char*)&buffer[difi_packet.offset+i], 2);
                int16_t img;
                memcpy(&img, (char*)&buffer[difi_packet.offset+i]+2, 2);
                signal[signal_pointer][REAL] = mult*re;
                signal[signal_pointer][IMAG] = mult*img;
                mult *= -1;
            }

            if (signal_pointer >= num_points) {

                signal_pointer = 0;

                fftw_execute(plan);

                double max = 0;
                int32_t max_i = -1;

               for (uint32_t i = 0; i < num_points; ++i) {
                    double mag = sqrt(result[i][REAL] * result[i][REAL] +
                              result[i][IMAG] * result[i][IMAG]);
                    // ignore 10% of bins around DC (exp.)
                    if ( (mag > max) and (not ignore_dc or (abs((int32_t)i-(int32_t)num_points/2) ) > num_points/10)  ) {
                        max = mag;
                        max_i = i;
                    }
                }

                uint64_t seconds = difi_packet.integer_seconds_timestamp;
                uint64_t frac_seconds = difi_packet.fractional_seconds_timestamp;
                frac_seconds += 10000*1e12/difi_context.sample_rate;
                if (frac_seconds > 1e12) {
                    frac_seconds -= 1e12;
                    seconds++;
                }

                int64_t peak_hz = difi_context.rf_freq + max_i - difi_context.sample_rate/2;
                printf("%lu.%09li, %li, %.3f\n", seconds, (int64_t)(frac_seconds/1e3), peak_hz, 20*log10(max/(double)num_points));
                fflush(stdout);
            }

            num_total_samps += difi_packet.num_rx_samps;

            if (start_rx and first_frame) {
                std::cout << boost::format(
                                 "# First frame: %u samples, %u full secs, %.09f frac secs")
                                 % difi_packet.num_rx_samps
                                 % difi_packet.integer_seconds_timestamp
                                 % ((double)difi_packet.fractional_seconds_timestamp/1e12)
                          << std::endl;
                first_frame = false;
                // Header
                printf("timestamp, frequency, power\n");
            }
        }

        if (progress) {
            if (difi_packet.data)
                last_update_samps += difi_packet.num_rx_samps;
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
                // if (cpu_format == "sc8" || cpu_format == "s8")
                //     datatype_max = 128.;

                for (int i=0; i<difi_packet.num_rx_samps; i++ ) {
                    auto sample_i = get_abs_val((std::complex<int16_t>)buffer[difi_packet.offset+i]);
                    sum_i += sample_i;
                    if (sample_i > datatype_max*0.99)
                        clip_i++;
                }
                sum_i = sum_i/difi_packet.num_rx_samps;
                std::cout << boost::format("%.0f") % (100.0*log2(sum_i)/log2(datatype_max)) << "% I (";
                std::cout << boost::format("%.0f") % ceil(log2(sum_i)+1) << " of ";
                std::cout << (int)ceil(log2(datatype_max)+1) << " bits), ";
                std::cout << "" << boost::format("%.0f") % (100.0*clip_i/difi_packet.num_rx_samps) << "% I clip, ";
                std::cout << std::endl;

            }
        }
    }

    zmq_close(subscriber);
    zmq_ctx_destroy(context);

    return 0;

}  
