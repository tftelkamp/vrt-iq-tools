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
// #include <fftw3.h>

#include "difi-tools.h"

namespace po = boost::program_options;

static bool stop_signal_called = false;
void sig_int_handler(int)
{
    stop_signal_called = true;
}

int main(int argc, char* argv[])
{

    // variables to be set by po
    std::string file, type, zmq_address;
    uint16_t port;
    int hwm;
    size_t num_requested_samples;
    double setup_time, freq, gain, lo_offset;

    // setup the program options
    po::options_description desc("Allowed options");
    // clang-format off

    desc.add_options()
        ("help", "help message")
        // ("file", po::value<std::string>(&file)->default_value("usrp_samples.dat"), "name of the file to write binary samples to")
        // ("type", po::value<std::string>(&type)->default_value("short"), "sample type: double, float, or short")
        // ("nsamps", po::value<size_t>(&num_requested_samples)->default_value(0), "total number of samples to receive")
        // ("duration", po::value<double>(&total_time)->default_value(0), "total number of seconds to receive")
        // ("progress", "periodically display short-term bandwidth")
        // ("int-second", "align start of reception to integer second")
        // ("null", "run without writing to file")
        ("setup", po::value<double>(&setup_time)->default_value(0.5), "seconds of setup time")
        ("freq", po::value<double>(&freq), "RF center frequency in Hz")
        ("gain", po::value<double>(&gain), "gain for the RF chain")
        ("lo-offset", po::value<double>(&lo_offset),"Offset for frontend LO in Hz (optional)")
        // ("continue", "don't abort on a bad packet")
        ("address", po::value<std::string>(&zmq_address)->default_value("localhost"), "DIFI ZMQ address")
        ("port", po::value<uint16_t>(&port)->default_value(50300), "DIFI ZMQ port")
        // ("hwm", po::value<int>(&hwm)->default_value(10000), "DIFI ZMQ HWM")
    ;
    // clang-format on
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // print the help message
    if (vm.count("help")) {
        std::cout << boost::format("DIFI samples to nothing. %s") % desc << std::endl;
        std::cout << std::endl
                  << "This application streams data from a DIFI stream "
                     "to nowhwere.\n"
                  << std::endl;
        return ~0;
    }

    bool progress               = vm.count("progress") > 0;
    bool stats                  = vm.count("stats") > 0;
    bool null                   = vm.count("null") > 0;
    bool continue_on_bad_packet = vm.count("continue") > 0;
    bool int_second             = (bool)vm.count("int-second");

    context_type difi_context;
    init_context(&difi_context);

    // ZMQ
    void *context = zmq_ctx_new();
    void *subscriber = zmq_socket(context, ZMQ_PUB);
    // int rc = zmq_setsockopt (subscriber, ZMQ_RCVHWM, &hwm, sizeof hwm);
    std::string connect_string = "tcp://" + zmq_address + ":" + std::to_string(port);
    int rc = zmq_connect(subscriber, connect_string.c_str());
    assert(rc == 0);
    // zmq_setsockopt(subscriber, ZMQ_SUBSCRIBE, "", 0);

    // Sleep setup time
    std::this_thread::sleep_for(std::chrono::milliseconds(int64_t(1000 * setup_time)));

    // VITA 49.2
    /* Initialize to reasonable values */
    struct vrt_packet pc;
    vrt_init_packet(&pc);

    uint32_t buffer[DIFI_DATA_PACKET_SIZE];

    /* DIFI Configure */
    difi_init_context_packet(&pc);

    pc.fields.stream_id = 0;

    struct timeval time_now{};
    gettimeofday(&time_now, nullptr);
    pc.fields.integer_seconds_timestamp = time_now.tv_sec;
    pc.fields.fractional_seconds_timestamp = 1e3*time_now.tv_usec;

    pc.if_context.has.bandwidth   = false;
    pc.if_context.has.sample_rate = false;
    pc.if_context.has.reference_point_identifier = true;
    pc.if_context.has.if_reference_frequency = false;
    pc.if_context.has.rf_reference_frequency = false;
    pc.if_context.has.if_band_offset = false;
    pc.if_context.has.reference_level = false;
    pc.if_context.has.gain = false;
    pc.if_context.has.timestamp_adjustment = false;
    pc.if_context.has.timestamp_calibration_time = false;
    pc.if_context.has.state_and_event_indicators = true;
    pc.if_context.has.data_packet_payload_format = true;
    pc.if_context.state_and_event_indicators.has.reference_lock = false;
    pc.if_context.state_and_event_indicators.has.calibrated_time = false;

    if (vm.count("gain")>0) {
        pc.if_context.has.gain = true;
        pc.if_context.gain.stage1 = gain;
        pc.if_context.gain.stage2 = 0;
    }

    if (vm.count("freq")>0) {
        pc.if_context.has.rf_reference_frequency = true;
        pc.if_context.rf_reference_frequency = freq;
    }

   if (vm.count("lo-offset")>0) {
        pc.if_context.has.if_band_offset = true;
        pc.if_context.if_band_offset = lo_offset;
    }

    int32_t rv = vrt_write_packet(&pc, buffer, DIFI_DATA_PACKET_SIZE, true);
    if (rv < 0) {
        fprintf(stderr, "Failed to write packet: %s\n", vrt_string_error(rv));
    }

    // ZMQ
    zmq_send (subscriber, buffer, rv*4, 0);

    zmq_close(subscriber);
    zmq_ctx_destroy(context);

    return 0;

}  
