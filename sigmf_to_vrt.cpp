//
// Copyright 2021/2022 by Thomas Telkamp 
//
// SPDX-License-Identifier: GPL-3.0-or-later
//

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/circular_buffer.hpp>

#include <chrono>
#include <complex>
#include <csignal>
#include <fstream>
#include <iostream>
#include <thread>

#include <sys/time.h>

#include <zmq.h>
#include <assert.h>

// VRT
#include <vrt/vrt_init.h>
#include <vrt/vrt_string.h>
#include <vrt/vrt_types.h>
#include <vrt/vrt_util.h>
#include <vrt/vrt_write.h>
#include <vrt/vrt_read.h>

// TCP/UDP
#include <sys/socket.h> 
#include <arpa/inet.h> 
#include <netinet/in.h> 
#include <netinet/tcp.h>
#include <netinet/udp.h>
#include <netdb.h>

// #include <boost/thread/thread.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
// #include <boost/date_time.hpp>

// Short alias for this namespace
namespace pt = boost::property_tree;

// DIFI tools functions
#include "difi-tools.h"

unsigned long long num_total_samps = 0;

int sockfd, connfd;

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
    std::string udp_forward, ref, file, time_cal, type, start_time_str;
    size_t total_num_samps;
    uint16_t port;
    uint32_t stream_id;
    int hwm;
    int16_t gain;
    double datarate;
    double rate, freq, bw, total_time, setup_time, lo_offset;

    FILE *read_ptr;
    FILE *read_ptr_2;

    // setup the program options
    po::options_description desc("Allowed options");
    // clang-format off
    desc.add_options()
        ("help", "help message")
        ("file", po::value<std::string>(&file)->default_value("samples.sigmf-meta"), "name of the SigMF meta file")
        ("setup", po::value<double>(&setup_time)->default_value(1.0), "seconds of setup time")
        ("datarate", po::value<double>(&datarate)->default_value(0), "rate of outgoing samples")
        ("dual-chan", "use two SigMF files for dual channel stream (chan0+chan1)")
        ("progress", "periodically display short-term bandwidth")
        ("stats", "show average bandwidth on exit")
        ("null", "run without streaming")
        ("continue", "don't abort on a bad packet")
        ("port", po::value<uint16_t>(&port)->default_value(50100), "DIFI ZMQ port")
        ("hwm", po::value<int>(&hwm)->default_value(10000), "DIFI ZMQ HWM")
    ;
    // clang-format on
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // print the help message
    if (vm.count("help")) {
        std::cout << boost::format("SigMF samples to Vita49 %s") % desc << std::endl;
        std::cout << std::endl
                  << "This application streams data from a SigMF file "
                     "to Vita49.\n"
                  << std::endl;
        return ~0;
    }

    bool bw_summary             = vm.count("progress") > 0;
    bool stats                  = vm.count("stats") > 0;
    bool null                   = vm.count("null") > 0;
    bool continue_on_bad_packet = vm.count("continue") > 0;
    bool dual_chan              = vm.count("dual-chan") > 0;

    struct timeval time_now{};
    gettimeofday(&time_now, nullptr);

    // seed random generator with seconds and microseconds
    srand(time_now.tv_usec + time_now.tv_sec);

    // Create ptree root
    pt::ptree root;

    // Load the json file in this ptree
    pt::read_json(file, root);

    rate = root.get<double>("global.core:sample_rate", 0);
    bw = root.get<double>("global.difi:bandwidth", 0);
    gain = root.get<int>("global.difi:rx_gain", 0);
    ref = root.get<std::string>("global.difi:reference", "");
    time_cal = root.get<std::string>("global.difi:time_source", "");
    type = root.get<std::string>("global.core:datatype", "");
    stream_id = root.get<uint32_t>("global.difi:stream_id", 0);

    for (auto& item : root.get_child("captures")) {
        freq = item.second.get<double>("core:frequency");
        start_time_str = item.second.get<std::string>("core:datetime");
    }

    printf("SigMF meta data:\n");
    printf("    Start time: %s\n", start_time_str.c_str());
    printf("    Sample rate: %i\n", (int)rate);
    printf("    Frequency: %i\n", (int)freq);
    printf("    Reference: %s\n", ref.c_str());
    printf("    Time calibration: %s\n", time_cal.c_str());
    printf("    Data type: %s\n", type.c_str());
    printf("    Stream ID: %u\n", stream_id);

    // Some Checks
    
    if (rate == 0 || freq == 0) {
        printf("Frequency and sample rate need to be set.\n");
        exit(1);
    }

    if (type != "ci16_le") {
        printf("Only 16 bit complex int data formmat supported (\"ci16_le\")\n");
        exit(1);
    }

    if (datarate == 0)
        datarate = rate;

    // Open data file

    std::string data_filename;
    boost::filesystem::path base_fn_fp(file);
    base_fn_fp.replace_extension(".sigmf-data");
    data_filename = base_fn_fp.string();

    printf("SigMF Data Filename: %s\n", data_filename.c_str());

    if (data_filename.c_str())
        read_ptr = fopen(data_filename.c_str(),"rb");  // r for read, b for binary

    if (dual_chan) {
        std::string data_filename_2(data_filename);
        boost::replace_all(data_filename_2,"chan0","chan1");
        printf("Second SigMF Data Filename: %s\n", data_filename_2.c_str());
        read_ptr_2 = fopen(data_filename_2.c_str(),"rb");  // r for read, b for binary
    }

	size_t samps_per_buff = DIFI_SAMPLES_PER_PACKET;

	unsigned long long num_requested_samples = total_num_samps;
    double time_requested = total_time;

    uint32_t buffer[DIFI_DATA_PACKET_SIZE];
   
    bool first_frame = true;

    struct vrt_packet p;
    vrt_init_packet(&p);

    /* Warn if not standards compliant */
    if (vrt_is_platform_little_endian()) {
        printf("Warning: little endian support is work in progress.\n");
    }

    /* DIFI init */
    difi_init_data_packet(&p);
    
    // p.fields.stream_id = stream_id;
    
    // ZMQ
    void *zmq_server;
  
    void *context = zmq_ctx_new();
    void *responder = zmq_socket(context, ZMQ_PUB);
    int rc = zmq_setsockopt (responder, ZMQ_SNDHWM, &hwm, sizeof hwm);
    assert(rc == 0);

    std::string connect_string = "tcp://*:" + std::to_string(port);
    rc = zmq_bind(responder, connect_string.c_str());
    assert (rc == 0);
    zmq_server = responder;

    // Sleep setup time
    std::this_thread::sleep_for(std::chrono::milliseconds(int64_t(1000 * setup_time)));

    if (total_num_samps == 0) {
        std::signal(SIGINT, &sig_int_handler);
        std::cout << "Press Ctrl + C to stop streaming..." << std::endl;
    }

    // time keeping
    auto start_time = std::chrono::steady_clock::now();

    // Track time and samps between updating the BW summary
    auto last_update                     = start_time;
    auto last_context                    = start_time;
    unsigned long long last_update_samps = 0;

    // Run this loop until either time expired (if a duration was given), until
    // the requested number of samples were collected (if such a number was
    // given), or until Ctrl-C was pressed.

    uint32_t frame_count = 0;
    uint32_t num_words_read=0;

    std::complex<short> samples[samps_per_buff];

    timeval time_first_sample;

    boost::posix_time::ptime t1(boost::posix_time::from_iso_extended_string(start_time_str));

    time_t integer_time_first_sample = boost::posix_time::to_time_t(t1);
    boost::posix_time::ptime t2(boost::posix_time::from_time_t(integer_time_first_sample));

    boost::posix_time::time_duration fractional_sec = t1-t2;
    time_first_sample.tv_sec = integer_time_first_sample;
    time_first_sample.tv_usec = fractional_sec.total_microseconds();

    auto difi_time = time_first_sample;

    // trigger context update
    last_context -= std::chrono::seconds(4*VRT_CONTEXT_INTERVAL);

    int update_interval = 1e6*samps_per_buff/datarate;

    printf("Update interval: %i\n", update_interval);

    while (not stop_signal_called) {
 
        const auto now = std::chrono::steady_clock::now();

        // wait
        auto wait_time = start_time + std::chrono::microseconds(frame_count*update_interval) - now;

        std::this_thread::sleep_for(wait_time);   

        const auto time_since_last_context = now - last_context;
        if (time_since_last_context > std::chrono::milliseconds(VRT_CONTEXT_INTERVAL)) {

            last_context = now;

            // VITA 49.2
            /* Initialize to reasonable values */
            struct vrt_packet pc;
            vrt_init_packet(&pc);

            /* DIFI Configure. Note that context packets cannot have a trailer word. */
            difi_init_context_packet(&pc);

            pc.fields.integer_seconds_timestamp = difi_time.tv_sec;
            pc.fields.fractional_seconds_timestamp = 1e3*difi_time.tv_usec;

            pc.fields.stream_id = 1;

            pc.if_context.bandwidth                         = bw;
            pc.if_context.sample_rate                       = rate;
            pc.if_context.rf_reference_frequency            = freq;
            pc.if_context.rf_reference_frequency_offset     = 0;
            pc.if_context.if_reference_frequency            = 0; // Zero-IF
            pc.if_context.gain.stage1                       = gain;
            pc.if_context.gain.stage2                       = 0;

            pc.if_context.state_and_event_indicators.reference_lock = (bool)(ref=="external") ? true : false;

            pc.if_context.state_and_event_indicators.calibrated_time = (bool)(time_cal=="external" || time_cal=="pps") ? true : false;

            int32_t rv = vrt_write_packet(&pc, buffer, DIFI_DATA_PACKET_SIZE, true);
            if (rv < 0) {
                fprintf(stderr, "Failed to write packet: %s\n", vrt_string_error(rv));
            }
            zmq_send (zmq_server, buffer, rv*4, 0);

            if (dual_chan) {
                // duplicate context of channel 0 on channel 1
                pc.fields.stream_id = 2;
                rv = vrt_write_packet(&pc, buffer, DIFI_DATA_PACKET_SIZE, true);
                if (rv < 0) {
                    fprintf(stderr, "Failed to write packet: %s\n", vrt_string_error(rv));
                }
                zmq_send (zmq_server, buffer, rv*4, 0);   
            }

        }

        // Read

        if (fread(samples, sizeof(samples), 1, read_ptr) == 1) {

            num_words_read = samps_per_buff;

            struct timeval interval_time;
            int64_t first_sample = frame_count*samps_per_buff;

            double interval = (double)first_sample/(double)rate;
            interval_time.tv_sec = (time_t)interval;
            interval_time.tv_usec = (interval-(time_t)interval)*1e6;

            timeradd(&time_first_sample, &interval_time, &difi_time);

            if (first_frame) {
                std::cout << boost::format(
                                 "First frame: %u samples, %u full secs, %.09f frac secs")
                                 % (num_words_read) % difi_time.tv_sec
                                 % (difi_time.tv_usec/1e6)
                          << std::endl;
                first_frame = false;
            }

            p.fields.stream_id = 1;
            p.body = samples;
            p.header.packet_count = (uint8_t)frame_count%16;
            p.fields.integer_seconds_timestamp = difi_time.tv_sec;
            p.fields.fractional_seconds_timestamp = 1e6*difi_time.tv_usec;
    
            zmq_msg_t msg;
            int rc = zmq_msg_init_size (&msg, DIFI_DATA_PACKET_SIZE*4);

            int32_t rv = vrt_write_packet(&p, zmq_msg_data(&msg), DIFI_DATA_PACKET_SIZE, true);

            zmq_msg_send(&msg, zmq_server, 0);
            zmq_msg_close(&msg);

            if (dual_chan) {
                if (fread(samples, sizeof(samples), 1, read_ptr_2) == 1) {
                    p.fields.stream_id = 2;
                    p.body = samples;
                    p.header.packet_count = (uint8_t)frame_count%16;
                    p.fields.integer_seconds_timestamp = difi_time.tv_sec;
                    p.fields.fractional_seconds_timestamp = 1e6*difi_time.tv_usec;
            
                    zmq_msg_t msg;
                    int rc = zmq_msg_init_size (&msg, DIFI_DATA_PACKET_SIZE*4);

                    int32_t rv = vrt_write_packet(&p, zmq_msg_data(&msg), DIFI_DATA_PACKET_SIZE, true);

                    zmq_msg_send(&msg, zmq_server, 0);
                    zmq_msg_close(&msg);
                } else {
                    break;
                }

            }

            frame_count++;

            if (bw_summary) {
                last_update_samps += num_words_read;
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

                    for (int i=0; i<samps_per_buff; i++ ) {
                        auto sample_i = get_abs_val(samples[i]);
                        sum_i += sample_i;
                        if (sample_i > datatype_max*0.99)
                            clip_i++;
                    }
                    sum_i = sum_i/10000;
                    std::cout << boost::format("%.0f") % (100.0*log2(sum_i)/log2(datatype_max)) << "% I (";
                    std::cout << boost::format("%.0f") % ceil(log2(sum_i)+1) << " of ";
                    std::cout << (int)ceil(log2(datatype_max)+1) << " bits), ";
                    std::cout << "" << boost::format("%.0f") % (100.0*clip_i/10000) << "% I clip.";
                    std::cout << std::endl;

                }
            }
        } else {
            // printf("no more samples in data file\n");
            break;
        }
    }

    const auto actual_stop_time = std::chrono::steady_clock::now();


    if (stats) {
        std::cout << std::endl;
        const double actual_duration_seconds =
            std::chrono::duration<float>(actual_stop_time - start_time).count();

        std::cout << boost::format("Received %d samples in %f seconds.") % num_total_samps
                         % actual_duration_seconds
                  << std::endl;
        const double rate = (double)num_total_samps / actual_duration_seconds;
        std::cout << (rate / 1e6) << " Msps." << std::endl;
    }
  
    /* clean up */
    fclose(read_ptr);

    // Sleep setup time
    std::this_thread::sleep_for(std::chrono::milliseconds(int64_t(1000 * setup_time)));
  
    // finished
    std::cout << std::endl << "Done!" << std::endl << std::endl;

    return EXIT_SUCCESS;
}
