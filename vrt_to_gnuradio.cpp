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

#include <vrt/vrt_read.h>
#include <vrt/vrt_string.h>
#include <vrt/vrt_types.h>
#include <vrt/vrt_util.h>

#include "vrt-tools.h"

// gnuradio pmt
#include <pmt/pmt.h>

namespace po = boost::program_options;

static bool stop_signal_called = false;
void sig_int_handler(int)
{
    stop_signal_called = true;
}

//! Change to filename, e.g. from usrp_samples.dat to usrp_samples.chan0.dat,
//  but only if multiple names are to be generated.
std::string generate_out_filename(
    const std::string& base_fn, size_t n_names, size_t this_name)
{
    if (n_names == 1) {
        return base_fn;
    }

    boost::filesystem::path base_fn_fp(base_fn);
    base_fn_fp.replace_extension(boost::filesystem::path(
        str(boost::format("chan%d%s") % this_name % base_fn_fp.extension().string())));
    return base_fn_fp.string();
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
    size_t num_requested_samples;
    double total_time;
    uint16_t port;
    uint32_t channel;
    int hwm;

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
        ("channel", po::value<uint32_t>(&channel)->default_value(0), "VRT channel")
        // ("stats", "show average bandwidth on exit")
        // ("int-second", "align start of reception to integer second")
        ("null", "run without writing to file")
        ("continue", "don't abort on a bad packet")
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
        std::cout << "Split VRT stream into metadata and data for GNURadio\n"
                  << "Output streams: port + 10 (default 50110): data\n"
                  << "                port + 20 (default 50120): frequency\n"
                  << "                port + 30 (default 50130): bandwidth\n"
                  << desc << std::endl;
        std::cout << std::endl
                  << "This application streams data from a VRT stream "
                     "to a ZeroMQ socket to be used in GNURadio.\n"
                  << std::endl;
        return ~0;
    }

    bool progress               = vm.count("progress") > 0;
    bool stats                  = vm.count("stats") > 0;
    bool null                   = vm.count("null") > 0;
    bool continue_on_bad_packet = vm.count("continue") > 0;

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

    void *zmq_gr_data;
    void *zmq_gr_freq;
    void *zmq_gr_rate;

    void *responder = zmq_socket(context, ZMQ_PUB);
    rc = zmq_bind(responder, ("tcp://*:" + std::to_string(port + 10)).c_str());
    assert (rc == 0);
    zmq_gr_data = responder;

    responder = zmq_socket(context, ZMQ_PUB);
    rc = zmq_bind(responder, ("tcp://*:" + std::to_string(port + 11)).c_str());
    assert (rc == 0);
    zmq_gr_freq = responder;

    responder = zmq_socket(context, ZMQ_PUB);
    rc = zmq_bind(responder, ("tcp://*:" + std::to_string(port + 12)).c_str());
    assert (rc == 0);
    zmq_gr_rate = responder;

    bool first_frame = true;

    // time keeping
    auto start_time = std::chrono::steady_clock::now();
    auto stop_time =
        start_time + std::chrono::milliseconds(int64_t(1000 * total_time));

    uint32_t buffer[ZMQ_BUFFER_SIZE];
    
    unsigned long long num_total_samps = 0;

    // Track time and samps between updating the BW summary
    auto last_update                     = start_time;
    auto last_gnuradio_forward           = start_time;
    unsigned long long last_update_samps = 0;

    bool start_rx = false;

    bool alternate = true;

    while (not stop_signal_called
           and (num_requested_samples > num_total_samps or num_requested_samples == 0)
           and (total_time == 0.0 or std::chrono::steady_clock::now() <= stop_time)) {

        int len = zmq_recv(subscriber, buffer, ZMQ_BUFFER_SIZE, 0);

        const auto now = std::chrono::steady_clock::now();

        if (not vrt_process(buffer, sizeof(buffer), &vrt_context, &vrt_packet)) {
            printf("Not a Vita49 packet?\n");
            continue;
        }

        if (vrt_packet.context) {
            if (not start_rx)
                vrt_print_context(&vrt_context);

            start_rx = true;

            const auto time_since_last_forward = now - last_gnuradio_forward;
            if (!start_rx or time_since_last_forward > std::chrono::seconds(1)) {

                last_gnuradio_forward = now;

                // Workaround to prevent GR 3.10 from crashing...
                if (alternate) {
                    pmt::pmt_t P_str = pmt::intern("freq");
                    pmt::pmt_t P_double = pmt::from_double(vrt_context.rf_freq);
                    pmt::pmt_t P_pair = pmt::cons(P_str, P_double);
                    std::string str = pmt::serialize_str(P_pair);
                    zmq_send (zmq_gr_freq, str.c_str(), str.size(), 0);
                } else {
                    pmt::pmt_t P_str = pmt::intern("sample_rate");
                    pmt::pmt_t P_double = pmt::from_double(vrt_context.sample_rate);
                    pmt::pmt_t P_pair = pmt::cons(P_str, P_double);
                    std::string str = pmt::serialize_str(P_pair);
                    zmq_send (zmq_gr_rate, str.c_str(), str.size(), 0);
                }
                alternate = !alternate;
                start_rx = true;
            }
        }

        if (start_rx and vrt_packet.data) {

            if (vrt_packet.lost_frame)
               if (not continue_on_bad_packet)
                    break;

            zmq_send(zmq_gr_data, (const char*)&buffer[vrt_packet.offset], sizeof(uint32_t)*vrt_packet.num_rx_samps, 0);

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
    zmq_close(zmq_gr_data);
    zmq_close(zmq_gr_rate);
    zmq_close(zmq_gr_freq);
    zmq_ctx_destroy(context);

    return 0;

}  
