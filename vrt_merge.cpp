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

#include "vrt-tools.h"

#include <chrono>
#include <csignal>
#include <fstream>
#include <iostream>
#include <thread>

#include <complex.h>

namespace po = boost::program_options;

static bool stop_signal_called = false;
void sig_int_handler(int) {
    std::cout<<"Ctrl-C pressed"<<std::endl;
    stop_signal_called = true;
}

int main(int argc, char* argv[])
{
    std::signal(SIGINT, &sig_int_handler);

    // variables to be set by po
    std::string zmq_address, merge_address;
    uint16_t port, pub_port, merge_port;
    int hwm;

    bool merge;

    // setup the program options
    po::options_description desc("Allowed options");
    // clang-format off

    desc.add_options()
        ("help", "help message")
        ("address", po::value<std::string>(&zmq_address)->default_value("localhost"), "(VRT) ZMQ address")
        ("port", po::value<uint16_t>(&port)->default_value(50100), "(VRT) ZMQ SUB port")
        ("pub-port", po::value<uint16_t>(&pub_port)->default_value(50101), "VRT ZMQ PUB port")
        ("merge", po::value<bool>(&merge)->default_value(false), "Merge another VRT ZMQ stream (SUB connect)")
        ("merge-port", po::value<uint16_t>(&merge_port)->default_value(50011), "VRT ZMQ merge port")
        ("merge-address", po::value<std::string>(&merge_address)->default_value("localhost"), "VRT ZMQ merge address")
        ("hwm", po::value<int>(&hwm)->default_value(10000), "(VRT) ZMQ HWM")

    ;
    // clang-format on
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // print the help message
    if (vm.count("help")) {
        std::cout << boost::format("Forward and optionally merge streams. %s") % desc << std::endl;
        std::cout << std::endl
                  << "This application forwards and optionally merges ZMQ VRT streams.\n"
                  << std::endl;
        return ~0;
    }

    void *context = zmq_ctx_new();

    void *subscriber = zmq_socket(context, ZMQ_SUB);
    int rc = zmq_setsockopt (subscriber, ZMQ_RCVHWM, &hwm, sizeof hwm);
    std::string connect_string = "tcp://" + zmq_address + ":" + std::to_string(port);
    rc = zmq_connect(subscriber, connect_string.c_str());
    assert(rc == 0);
    zmq_setsockopt(subscriber, ZMQ_SUBSCRIBE, "", 0);

    void *merge_zmq = zmq_socket(context, ZMQ_SUB);
    if (merge) {
        connect_string = "tcp://" + merge_address + ":" + std::to_string(merge_port);
        rc = zmq_connect(merge_zmq, connect_string.c_str());
        assert(rc == 0);
        zmq_setsockopt(merge_zmq, ZMQ_SUBSCRIBE, "", 0);
    }

    void *publisher = zmq_socket(context, ZMQ_PUB);
    rc = zmq_setsockopt (publisher, ZMQ_SNDHWM, &hwm, sizeof hwm);
    connect_string = "tcp://*:" + std::to_string(pub_port);
    rc = zmq_bind(publisher, connect_string.c_str());
    assert(rc == 0);

    zmq_pollitem_t pollitems[2];
    int num_subscribers = 1;

    pollitems[0] = { subscriber, 0, ZMQ_POLLIN, 0 };

    if (merge) {
        pollitems[1] = { merge_zmq, 0, ZMQ_POLLIN, 0 };
        num_subscribers = 2;
    }

    zmq_msg_t message, message_meta;

    // Flush merge queue
    if (merge) {
        zmq_msg_init(&message_meta);
        while (zmq_msg_recv(&message, merge_zmq, ZMQ_DONTWAIT) > 0) {
            zmq_msg_close(&message);
            zmq_msg_init(&message);
        }
    }

    while (not stop_signal_called) {
        rc = zmq_poll(pollitems, num_subscribers, -1);
        assert(rc >= 0);

        if (pollitems[0].revents & ZMQ_POLLIN) {
            zmq_msg_init(&message);
            assert(rc != -1);
            rc = zmq_msg_recv(&message, subscriber, 0);
            assert(rc != -1);
            rc = zmq_msg_send(&message, publisher, 0);
            assert(rc != -1);
            zmq_msg_close(&message);
        }

        if (merge && (pollitems[1].revents & ZMQ_POLLIN)) {
            zmq_msg_init(&message_meta);
            rc = zmq_msg_recv(&message_meta, merge_zmq, 0);
            assert(rc != -1);
            rc = zmq_msg_send(&message_meta, publisher, 0);
            assert(rc != -1);
            zmq_msg_close(&message_meta);
        }
    }

    //zmq_proxy(subscriber, publisher, NULL);

    zmq_close(subscriber);

    if (merge)
        zmq_close(merge_zmq);

    zmq_close(publisher);
    zmq_ctx_destroy(context);

    return 0;
}
