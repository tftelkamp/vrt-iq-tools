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
#include <csignal>
#include <fstream>
#include <iostream>
#include <thread>

#include <complex.h>

namespace po = boost::program_options;

int main(int argc, char* argv[])
{

    // variables to be set by po
    std::string zmq_address;
    uint16_t port, pub_port;
    int hwm;

    // setup the program options
    po::options_description desc("Allowed options");
    // clang-format off

    desc.add_options()
        ("help", "help message")
        ("address", po::value<std::string>(&zmq_address)->default_value("localhost"), "(DIFI) ZMQ address")
        ("port", po::value<uint16_t>(&port)->default_value(50100), "(DIFI) ZMQ SUB port")
        ("pub-port", po::value<uint16_t>(&pub_port)->default_value(50101), "DIFI ZMQ PUB port")
        ("hwm", po::value<int>(&hwm)->default_value(10000), "(DIFI) ZMQ HWM")

    ;
    // clang-format on
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // print the help message
    if (vm.count("help")) {
        std::cout << boost::format("(DIFI) ZMQ forwarder %s") % desc << std::endl;
        std::cout << std::endl
                  << "This application forwards "
                     "(DIFI) ZMQ.\n"
                  << std::endl;
        return ~0;
    }

    void *context = zmq_ctx_new();

    void *subscriber = zmq_socket(context, ZMQ_XSUB);
    int rc = zmq_setsockopt (subscriber, ZMQ_RCVHWM, &hwm, sizeof hwm);
    std::string connect_string = "tcp://" + zmq_address + ":" + std::to_string(port);
    rc = zmq_connect(subscriber, connect_string.c_str());
    assert(rc == 0);

    void *publisher = zmq_socket(context, ZMQ_XPUB);
    rc = zmq_setsockopt (publisher, ZMQ_SNDHWM, &hwm, sizeof hwm);
    connect_string = "tcp://*:" + std::to_string(pub_port);
    rc = zmq_bind(publisher, connect_string.c_str());
    assert(rc == 0);

    zmq_proxy(subscriber, publisher, NULL);

    zmq_close(subscriber);
    zmq_close(publisher);
    zmq_ctx_destroy(context);

    return 0;

}
