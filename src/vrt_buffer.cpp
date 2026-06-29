#include <zmq.h>
#include <stdio.h>
#include <assert.h>
#include <unistd.h>

#include <boost/format.hpp>
#include <boost/program_options.hpp>

#include <chrono>
#include <deque>
#include <csignal>
#include <iostream>

#include "vrt-tools.h"

namespace po = boost::program_options;

static bool stop_signal_called = false;
void sig_int_handler(int)
{
    stop_signal_called = true;
}

struct TimestampedMessage {
    std::chrono::time_point<std::chrono::steady_clock> timestamp;
    uint32_t data[VRT_DATA_PACKET_SIZE];
    int len;  // actual received length in bytes
};


int main(int argc, char* argv[])
{
    // variables to be set by po
    std::string file, type, zmq_address;
    uint16_t pub_instance, instance, main_port, port, pub_port;
    int hwm;

    uint16_t delay;

    // setup the program options
    po::options_description desc("Allowed options");
    // clang-format off

    desc.add_options()
        ("help", "help message")
        ("delay", po::value<uint16_t>(&delay)->default_value(1), "delay in seconds")
        ("address", po::value<std::string>(&zmq_address)->default_value("localhost"), "VRT ZMQ address")
        ("progress", "periodically display short-term bandwidth")
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
        std::cout << boost::format("VRT-buffer. %s") % desc << std::endl;
        std::cout << std::endl
                  << "This application replays a VRT stream with a delay.\n"
                  << std::endl;
        return ~0;
    }

    bool progress = vm.count("progress") > 0;

    if (vm.count("port") > 0) {
        main_port = port;
    } else {
        main_port = DEFAULT_MAIN_PORT + MAX_CHANNELS*instance;
    }

    if (vm.count("pub-port") == 0) {
        pub_port = DEFAULT_MAIN_PORT + MAX_CHANNELS*pub_instance;
    }

    std::deque<TimestampedMessage> delay_buffer;

    void *context = zmq_ctx_new();

    void *subscriber = zmq_socket(context, ZMQ_SUB);
    int rc = zmq_setsockopt(subscriber, ZMQ_RCVHWM, &hwm, sizeof hwm);
    assert(rc == 0);
    std::string connect_string = "tcp://" + zmq_address + ":" + std::to_string(main_port);
    rc = zmq_connect(subscriber, connect_string.c_str());
    assert(rc == 0);
    zmq_setsockopt(subscriber, ZMQ_SUBSCRIBE, "", 0);

    void *publisher = zmq_socket(context, ZMQ_PUB);
    rc = zmq_setsockopt(publisher, ZMQ_SNDHWM, &hwm, sizeof hwm);
    assert(rc == 0);
    connect_string = "tcp://*:" + std::to_string(pub_port);
    rc = zmq_bind(publisher, connect_string.c_str());
    assert(rc == 0);

    context_type vrt_context;
    init_context(&vrt_context);
    packet_type vrt_packet;

    vrt_packet.channel_filt = 1;

    bool start_rx = false;
    bool first_frame = true;

    // Track time and samps between updating the BW summary
    auto last_update = std::chrono::steady_clock::now();
    uint64_t last_update_samps = 0;

    std::signal(SIGINT, &sig_int_handler);
    std::cout << "Press Ctrl + C to stop streaming..." << std::endl;

    const auto delay_duration = std::chrono::milliseconds(int64_t(1000 * delay));

    while (not stop_signal_called) {
        while (true) {
            delay_buffer.emplace_back();
            TimestampedMessage& msg = delay_buffer.back();
            int len = zmq_recv(subscriber, msg.data, sizeof(msg.data), ZMQ_NOBLOCK);
            if (len < 0) { delay_buffer.pop_back(); break; }

            msg.timestamp = std::chrono::steady_clock::now();
            msg.len = len;

            if (not start_rx) {
                if (not vrt_process(msg.data, len, &vrt_context, &vrt_packet)) {
                    printf("Not a Vita49 packet?\n");
                    continue;
                }
                if (vrt_packet.context) {
                    vrt_print_context(&vrt_context);
                    start_rx = true;
                    // Possibly do something with context here
                }
            }

            if (progress) {
                if (not vrt_process(msg.data, len, &vrt_context, &vrt_packet)) {
                    printf("Not a Vita49 packet?\n");
                    continue;
                }
                if (vrt_packet.data) {
                    show_progress_stats(
                        msg.timestamp,
                        &last_update,
                        &last_update_samps,
                        &msg.data[vrt_packet.offset],
                        vrt_packet.num_rx_samps, 0
                    );
                }
            }

        }

        const auto now = std::chrono::steady_clock::now();
        while (!delay_buffer.empty()) {
            auto& front = delay_buffer.front();
            if (now - front.timestamp < delay_duration) break;  // rest are newer

            if (first_frame) {
                vrt_context.last_data_counter = 0;
                vrt_process(front.data, front.len, &vrt_context, &vrt_packet);
                if (vrt_packet.data) {
                    printf("# Start forwarding at %llu full secs, %.09f frac secs\n", vrt_packet.integer_seconds_timestamp, (double)vrt_packet.fractional_seconds_timestamp/1e12);
                    first_frame = false;
                    vrt_context.last_data_counter = 0;
                }
            }

            zmq_send(publisher, front.data, front.len, 0);
            delay_buffer.pop_front();
        }

        usleep(5);
    }

    zmq_close(subscriber);
    zmq_close(publisher);
    zmq_ctx_destroy(context);

    return 0;
}
