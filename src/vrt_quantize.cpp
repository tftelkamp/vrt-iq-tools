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
    std::string file, type, zmq_address;
    uint16_t pub_instance, instance, main_port, port, pub_port;
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
        ("unpack", "unpack stream")
        // ("zmq-split", "create a ZeroMQ stream per VRT channel, increasing port number for additional streams")
        ("address", po::value<std::string>(&zmq_address)->default_value("localhost"), "VRT ZMQ address")
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
        std::cout << boost::format("VRT quantizes. %s") % desc << std::endl;
        std::cout << std::endl
                  << "This application 1-bit quantizes a VRT stream.\n"
                  << std::endl;
        return ~0;
    }

    bool progress               = vm.count("progress") > 0;
    bool stats                  = vm.count("stats") > 0;
    bool null                   = vm.count("null") > 0;
    bool continue_on_bad_packet = vm.count("continue") > 0;
    bool int_second             = vm.count("int-second");
    bool zmq_split              = vm.count("zmq-split") > 0;
    bool unpack                 = vm.count("unpack") > 0;

    context_type vrt_context;
    init_context(&vrt_context);

    packet_type vrt_packet;

    // tracker_ext_context_type tracker_ext_context;

    if (vm.count("port") > 0) {
        main_port = port;
    } else {
        main_port = DEFAULT_MAIN_PORT + MAX_CHANNELS*instance;
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

    vrt_packet.channel_filt = 1<<channel;
    // vrt_packet.channel_filt |= 1<<(channel+1);

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

    uint32_t buffer[ZMQ_BUFFER_SIZE];
    uint32_t data_buffer[VRT_SAMPLES_PER_PACKET];

    uint64_t num_total_samps = 0;

    // Track time and samps between updating the BW summary
    auto last_update = start_time;
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

    bool first_context = true;

    uint32_t frame_count = 0;

    int len;

    while (not stop_signal_called
           and (num_requested_samples > num_total_samps or num_requested_samples == 0)
           and (total_time == 0.0 or std::chrono::steady_clock::now() <= stop_time)) {

        const auto now = std::chrono::steady_clock::now();

        len = zmq_recv(subscriber, buffer, ZMQ_BUFFER_SIZE, 0);

        if (not vrt_process(buffer, sizeof(buffer), &vrt_context, &vrt_packet)) {
            printf("Not a Vita49 packet?\n");
            continue;
        }

        if (not start_rx and vrt_packet.context) {
            vrt_print_context(&vrt_context);
            start_rx = true;
            // Possibly do something with context here
            // vrt_context
        }

        if (vrt_packet.context) {
            zmq_msg_t msg;
            zmq_msg_init_size (&msg, len);
            memcpy (zmq_msg_data(&msg), buffer, len);
            zmq_msg_send(&msg, responder, 0);
            zmq_msg_close(&msg);
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
        }

        if (progress && vrt_packet.data)
            show_progress_stats(
                now,
                &last_update,
                &last_update_samps,
                &buffer[vrt_packet.offset],
                vrt_packet.num_rx_samps, 0
            );

        if (start_rx and vrt_packet.data) {
 
            if (!unpack) {

                memset(data_buffer, 0, VRT_SAMPLES_PER_PACKET * sizeof(uint32_t));

                size_t word_idx;           // Which uint32_t
                size_t bit_idx;            // Which bit in that uint32_t

                for (uint32_t i = 0; i < vrt_packet.num_rx_samps; i++) {

                    uint32_t word;
                    int16_t value;

                    // real

                    word = buffer[vrt_packet.offset+i];

                    value = ((int16_t)(word & 0xFFFF));

                    word_idx = (i*2) / 32;
                    bit_idx = (i*2) % 32;

                    // imag
                    data_buffer[word_idx] |= !(value >> 15) << bit_idx;

                    value = (int16_t)((word >> 16) & 0xFFFF);

                    word_idx = (i*2+1) / 32;
                    bit_idx = (i*2+1) % 32;

                    data_buffer[word_idx] |= !(value >> 15) << bit_idx;

                }

                size_t new_data_len_words = vrt_packet.num_rx_samps / 16;
                memcpy((char*)&buffer[vrt_packet.offset], (char*)data_buffer, new_data_len_words * sizeof(uint32_t));
                len = (vrt_packet.offset + new_data_len_words)*4;

            } else {

                memset(data_buffer, 0, VRT_SAMPLES_PER_PACKET * sizeof(uint32_t));

                size_t word_idx; // Which uint32_t
                size_t bit_idx;  // Which bit in that uint32_t

                for (uint32_t i = 0; i < vrt_packet.num_rx_samps; i++) {

                    // real
                    word_idx = (i*2) / 32;
                    bit_idx = (i*2) % 32;
                    int16_t re = (( buffer[vrt_packet.offset+word_idx] >> bit_idx) & 1) ? 1 : -1;
                    memcpy((char*)&data_buffer[i], &re, 2);

                    // imag
                    word_idx = (i*2+1) / 32;
                    bit_idx = (i*2+1) % 32;
                    int16_t img = (( buffer[vrt_packet.offset+word_idx] >> bit_idx) & 1) ? 1 : -1;
                    memcpy((char*)&data_buffer[i]+2, &img, 2);

                }

                size_t new_data_len_words = vrt_packet.num_rx_samps;
                memcpy((char*)&buffer[vrt_packet.offset], (char*)data_buffer, new_data_len_words * sizeof(uint32_t));
                len = (vrt_packet.offset + new_data_len_words)*4;
            }

            zmq_msg_t msg;
            zmq_msg_init_size (&msg, len);
            memcpy (zmq_msg_data(&msg), buffer, len);
            zmq_msg_send(&msg, responder, 0);
            zmq_msg_close(&msg);

        }
    }

    zmq_close(subscriber);
    zmq_close(responder);
    zmq_ctx_destroy(context);

    return 0;

}
