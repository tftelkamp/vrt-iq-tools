#include <zmq.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <arpa/inet.h>

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/thread/thread.hpp>

#include <chrono>
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
#include "dt-extended-context.h"
#include "tracker-extended-context.h"

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
    std::string file, type, zmq_address1, zmq_address2;
    uint16_t pub_instance, instance1, main_port1, port1, instance2, main_port2, port2, pub_port;
    uint32_t channel;
    int hwm;

    size_t num_requested_samples;
    double total_time;

    // setup the program options
    po::options_description desc("Allowed options");
    // clang-format off

    desc.add_options()
        ("help", "help message")
        ("nsamps", po::value<size_t>(&num_requested_samples)->default_value(0), "total number of samples to receive")
        ("duration", po::value<double>(&total_time)->default_value(0), "total number of seconds to receive")
        ("channel", po::value<uint32_t>(&channel)->default_value(0), "VRT channel")
        ("progress", "periodically display short-term bandwidth")
        // ("int-second", "align start of reception to integer second")
        ("null", "run without writing to file")
        ("continue", "don't abort on a bad packet")
        // ("zmq-split", "create a ZeroMQ stream per VRT channel, increasing port number for additional streams")
        ("address1", po::value<std::string>(&zmq_address1)->default_value("localhost"), "VRT ZMQ address 1")
        ("instance1", po::value<uint16_t>(&instance1)->default_value(0), "VRT ZMQ instance 1")
        ("port1", po::value<uint16_t>(&port1), "VRT ZMQ port 1")
        ("address2", po::value<std::string>(&zmq_address2)->default_value("localhost"), "VRT ZMQ address 2")
        ("instance2", po::value<uint16_t>(&instance2)->default_value(0), "VRT ZMQ instance 2")
        ("port2", po::value<uint16_t>(&port2), "VRT ZMQ port 2")
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
        std::cout << boost::format("VRT merge. %s") % desc << std::endl;
        std::cout << std::endl
                  << "This application merges two VRT streams.\n"
                  << std::endl;
        return ~0;
    }

    bool progress               = vm.count("progress") > 0;
    bool stats                  = vm.count("stats") > 0;
    bool null                   = vm.count("null") > 0;
    bool continue_on_bad_packet = vm.count("continue") > 0;
    bool int_second             = vm.count("int-second");
    bool zmq_split              = vm.count("zmq-split") > 0;

    context_type vrt_context1;
    context_type vrt_context2;
    init_context(&vrt_context1);
    init_context(&vrt_context2);

    packet_type vrt_packet1;
    packet_type vrt_packet2;

    // tracker_ext_context_type tracker_ext_context;

    if (vm.count("port1") > 0) {
        main_port1 = port1;
    } else {
        main_port1 = DEFAULT_MAIN_PORT + MAX_CHANNELS*instance1;
    }
    if (vm.count("port2") > 0) {
        main_port2 = port2;
    } else {
        main_port2 = DEFAULT_MAIN_PORT + MAX_CHANNELS*instance2;
    }

    if (!(vm.count("pub-port") > 0)) {
        pub_port = DEFAULT_MAIN_PORT + MAX_CHANNELS*pub_instance;
    }

    // if (zmq_split) {
    //     main_port1 += channel;
    //     vrt_packet.channel_filt = 1;
    // } else {
    //     vrt_packet.channel_filt = 1<<channel;
    // }

    vrt_packet1.channel_filt = 1<<channel;
    vrt_packet2.channel_filt = 1<<channel;

    // ZMQ
    void *context = zmq_ctx_new();
    void *subscriber1 = zmq_socket(context, ZMQ_SUB);
    int rc = zmq_setsockopt (subscriber1, ZMQ_RCVHWM, &hwm, sizeof hwm);
    std::string connect_string = "tcp://" + zmq_address1 + ":" + std::to_string(main_port1);
    rc = zmq_connect(subscriber1, connect_string.c_str());
    assert(rc == 0);
    zmq_setsockopt(subscriber1, ZMQ_SUBSCRIBE, "", 0);

    void *subscriber2 = zmq_socket(context, ZMQ_SUB);
    rc = zmq_setsockopt (subscriber2, ZMQ_RCVHWM, &hwm, sizeof hwm);
    connect_string = "tcp://" + zmq_address2 + ":" + std::to_string(main_port2);
    rc = zmq_connect(subscriber2, connect_string.c_str());
    assert(rc == 0);
    zmq_setsockopt(subscriber2, ZMQ_SUBSCRIBE, "", 0);

    void *responder = zmq_socket(context, ZMQ_PUB);
    rc = zmq_setsockopt (responder, ZMQ_SNDHWM, &hwm, sizeof hwm);
    assert(rc == 0);
    connect_string = "tcp://*:" + std::to_string(pub_port);
    rc = zmq_bind(responder, connect_string.c_str());
    assert (rc == 0);

    // time keeping
    auto start_time = std::chrono::steady_clock::now();
    auto stop_time = start_time + std::chrono::milliseconds(int64_t(1000 * total_time));

    uint32_t rx_buffer[2][ZMQ_BUFFER_SIZE];
    uint32_t rx_stored[2][ZMQ_BUFFER_SIZE];

    uint32_t rx_stored_len[2] = {0};

    unsigned long long num_total_samps = 0;

    // Track time and samps between updating the BW summary
    auto last_update1                     = start_time;
    unsigned long long last_update_samps1 = 0;

    auto last_update2                     = start_time;
    unsigned long long last_update_samps2 = 0;

    bool first_frame = true;
    uint64_t last_fractional_seconds_timestamp = 0;

    // set to true to process data before context
    bool start_rx1 = false;
    bool start_rx2 = false;

    /* VRT init */
    struct vrt_packet p;
    vrt_init_packet(&p);
    vrt_init_data_packet(&p);
    p.fields.stream_id = 1;

    bool first_context = true;

    uint64_t t1_integer_seconds_timestamp = 0;
    uint64_t t1_fractional_seconds_timestamp = 0;
    uint64_t t2_integer_seconds_timestamp = 0;
    uint64_t t2_fractional_seconds_timestamp = 0;

    bool forwarded = false;

    double t1 = 0;
    double t2 = 0;

    int len1, len2;

    uint32_t frame_count = 0;

    while (not stop_signal_called
           and (num_requested_samples > num_total_samps or num_requested_samples == 0)
           and (total_time == 0.0 or std::chrono::steady_clock::now() <= stop_time)) {

        const auto now = std::chrono::steady_clock::now();

        len1 = 0;
        len2 = 0;

        if (t1 <= t2 or t1 == 0) {
            len1 = zmq_recv(subscriber1, rx_buffer[0], ZMQ_BUFFER_SIZE, ZMQ_NOBLOCK);
            // printf("wait 1\n");
        }

        if (t1 >= t2 or t2 == 0) {
            len2 = zmq_recv(subscriber2, rx_buffer[1], ZMQ_BUFFER_SIZE, ZMQ_NOBLOCK);
            // printf("wait 2\n");
        }

        if (len1 > 0) {
            if (not vrt_process(rx_buffer[0], sizeof(rx_buffer[0]), &vrt_context1, &vrt_packet1)) {
                printf("Not a Vita49 packet?\n");
                continue;
            }
            if (vrt_packet1.data) {
                if (vrt_packet1.lost_frame)
                    break;
                t1_integer_seconds_timestamp = vrt_packet1.integer_seconds_timestamp;
                t1_fractional_seconds_timestamp = vrt_packet1.fractional_seconds_timestamp;
                t1 = (double)vrt_packet1.integer_seconds_timestamp + (double)vrt_packet1.fractional_seconds_timestamp/1e12;
                forwarded = false;
                memcpy((char*)rx_stored[0], (char*)rx_buffer[0], len1);
                rx_stored_len[0] = len1;
            } else {
                zmq_msg_t msg;
                rx_buffer[0][1] = htonl(1u);
                zmq_msg_init_size (&msg, len1);
                memcpy (zmq_msg_data(&msg), rx_buffer[0], len1);
                zmq_msg_send(&msg, responder, 0);
                zmq_msg_close(&msg);
            }
            if (not start_rx1 and vrt_packet1.context) {
                vrt_print_context(&vrt_context1);
                start_rx1 = true;
            }
            if (progress && vrt_packet1.data)
                show_progress_stats(
                    now,
                    &last_update1,
                    &last_update_samps1,
                    &rx_buffer[0][vrt_packet1.offset],
                    vrt_packet1.num_rx_samps, 0
                );

        }

        if (len2 > 0) {
            if (not vrt_process(rx_buffer[1], sizeof(rx_buffer[1]), &vrt_context2, &vrt_packet2)) {
                printf("Not a Vita49 packet?\n");
                continue;
            }
            if (vrt_packet2.data) {
                if (vrt_packet2.lost_frame)
                    break;
                t2_integer_seconds_timestamp = vrt_packet2.integer_seconds_timestamp;
                t2_fractional_seconds_timestamp = vrt_packet2.fractional_seconds_timestamp;
                t2 = (double)vrt_packet2.integer_seconds_timestamp + (double)vrt_packet2.fractional_seconds_timestamp/1e12;
                forwarded = false;
                memcpy((char*)rx_stored[1], (char*)rx_buffer[1], len2);
                rx_stored_len[1] = len2;
            } else {
                zmq_msg_t msg;
                rx_buffer[1][1] = htonl(2u);
                zmq_msg_init_size (&msg, len2);
                memcpy (zmq_msg_data(&msg), rx_buffer[1], len2);
                zmq_msg_send(&msg, responder, 0);
                zmq_msg_close(&msg);
            }
            if (not start_rx2 and vrt_packet2.context) {
                vrt_print_context(&vrt_context2);
                start_rx2 = true;
            }
            if (progress && vrt_packet1.data)
                show_progress_stats(
                    now,
                    &last_update2,
                    &last_update_samps2,
                    &rx_buffer[1][vrt_packet2.offset],
                    vrt_packet2.num_rx_samps, 1
                );
        }

        if ((len1 >0 or len2>0) and !forwarded and t1 == t2 and t1 >0 and t2 >0){

            forwarded = true;
            if (first_frame) {
                printf("# Start forwarding at %llu full secs, %.09f frac secs\n", t1_integer_seconds_timestamp, (double)t1_fractional_seconds_timestamp/1e12);
                first_frame = false;
            }

            // channel 0
            uint32_t header = uint32_t ntohl(rx_stored[0][0]);
            header = header & 0xFFF0FFFF;
            header = header | ((uint32_t)frame_count%16) << 16;
            rx_stored[0][0] = htonl(header);
            rx_stored[0][1] = htonl(1u);
            frame_count++;

            zmq_msg_t msg;
            rx_stored[0][1] = htonl(1u);
            zmq_msg_init_size (&msg, rx_stored_len[0]);
            memcpy (zmq_msg_data(&msg), rx_stored[0], rx_stored_len[0]);
            zmq_msg_send(&msg, responder, 0);
            zmq_msg_close(&msg);

            // channel 1
            header = uint32_t ntohl(rx_stored[1][0]);
            header = header & 0xFFF0FFFF;
            header = header | ((uint32_t)frame_count%16) << 16;
            rx_stored[1][0] = htonl(header);
            rx_stored[1][1] = htonl(2u);
            frame_count++;

            zmq_msg_init_size (&msg, rx_stored_len[1]);
            memcpy (zmq_msg_data(&msg), rx_stored[1], rx_stored_len[1]);
            zmq_msg_send(&msg, responder, 0);
            zmq_msg_close(&msg);
        }

        // if (not vrt_packet.context and not vrt_packet.data)
        //     continue;

    }

    zmq_close(subscriber1);
    zmq_close(subscriber2);
    zmq_close(responder);
    zmq_ctx_destroy(context);

    return 0;

}
