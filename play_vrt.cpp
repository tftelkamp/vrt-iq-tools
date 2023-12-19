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
#include <boost/date_time/posix_time/posix_time.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

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

// Short alias for this namespace
namespace pt = boost::property_tree;

// VRT tools functions
#include "vrt-tools.h"

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

unsigned long long num_total_samps = 0;

namespace po = boost::program_options;

static bool stop_signal_called = false;
static bool last_frame = false;

void sig_int_handler(int)
{
    stop_signal_called = true;
    last_frame = true;
}

int main(int argc, char* argv[])
{
    // variables to be set by po
    std::string ref, file, time_cal, type, start_time_str, zmq_address;
    size_t total_num_samps, tx_int;
    uint16_t port, tx_buffer_size;
    uint32_t stream_id;
    int hwm;
    uint16_t gain, tx_gain;
    double datarate;
    double rate, freq, bw, total_time, setup_time, lo_offset, tx_freq;

    FILE *read_ptr;
    FILE *read_ptr_2;

    datarate = 0;
    tx_freq = 0;
    tx_gain = 0;

    // setup the program options
    po::options_description desc("Allowed options");
    // clang-format off
    desc.add_options()
        ("help", "help message")
        ("file", po::value<std::string>(&file)->default_value("samples.sigmf-meta"), "name of the SigMF meta file")
        ("setup", po::value<double>(&setup_time)->default_value(1.0), "seconds of setup time")
        ("datarate", po::value<double>(&datarate), "rate of outgoing samples")
        ("tx-buffer", po::value<uint16_t>(&tx_buffer_size)->default_value(10), "VRT ZMQ transmit buffer size")
        ("tx-freq", po::value<double>(&tx_freq), "TX RF center frequency in Hz")
        ("tx-gain", po::value<uint16_t>(&tx_gain), "gain for the TX RF chain")
        // ("dual-chan", "use two SigMF files for dual channel stream (chan0+chan1)")
        ("progress", "periodically display short-term bandwidth")
        ("timed-tx", "Start transmission at given time (SigMF)")
        ("tx-interval", po::value<size_t>(&tx_int)->default_value(0), "start transmission at multiple of interval")
        ("stats", "show average bandwidth on exit")
        ("null", "run without streaming")
        ("continue", "don't abort on a bad packet")
        ("stdin", "read stream from stdin")
        ("repeat", "repeat the input file")
        ("port", po::value<uint16_t>(&port)->default_value(50500), "VRT ZMQ transmit port")
        ("address", po::value<std::string>(&zmq_address)->default_value("localhost"), "VRT ZMQ transmit address")
        ("hwm", po::value<int>(&hwm)->default_value(10000), "VRT ZMQ HWM")
    ;

    // clang-format on
    po::positional_options_description parser_positional;
    parser_positional.add("file", -1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(parser_positional).run(), vm);
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
    bool repeat                 = vm.count("repeat") > 0;
    bool read_stdin             = vm.count("stdin") > 0;
    bool timed_tx               = vm.count("timed-tx") > 0;
    bool send_context           = true;

    size_t filesize = 0;

    struct timeval time_now{};
    gettimeofday(&time_now, nullptr);

    // seed random generator with seconds and microseconds
    srand(time_now.tv_usec + time_now.tv_sec);

    if (not read_stdin) {
        // Create ptree root
        pt::ptree root;

        // Load the json file in this ptree
        std::string meta_filename;
        boost::filesystem::path base_fn_fp(file);
        base_fn_fp.replace_extension(".sigmf-meta");
        meta_filename = base_fn_fp.string();
        pt::read_json(meta_filename, root);

        rate = root.get<double>("global.core:sample_rate", datarate);
        bw = root.get<double>("global.vrt:bandwidth", 0);
        gain = root.get<int>("global.vrt:tx_gain", tx_gain);
        type = root.get<std::string>("global.core:datatype", "");
        stream_id = root.get<uint32_t>("global.vrt:stream_id", 0);

        for (auto& item : root.get_child("captures")) {
            freq = item.second.get<double>("core:frequency", tx_freq);
            start_time_str = item.second.get<std::string>("core:datetime");
        }

        printf(" SigMF meta data:\n");
        printf("      Start time: %s\n", start_time_str.c_str());
        printf("     Sample rate: %i\n", (int)rate);
        printf("            Gain: %i\n", (int)gain);
        printf("       Frequency: %i\n", (int)freq);
        printf("       Data type: %s\n", type.c_str());
        printf("       Stream ID: %u\n", stream_id);

        // Some Checks

        if (type != "ci16_le") {
            printf("Only 16 bit complex int data format supported (\"ci16_le\")\n");
            exit(1);
        }

        if (datarate == 0)
            datarate = rate;

        // Open data file

        std::string data_filename;
        base_fn_fp.replace_extension(".sigmf-data");
        data_filename = base_fn_fp.string();

        printf("SigMF Data Filename: %s\n", data_filename.c_str());

        if (data_filename.c_str()) {
            read_ptr = fopen(data_filename.c_str(),"rb");  // r for read, b for binary
            fseek(read_ptr, 0L, SEEK_END);    // seek to the EOF
            filesize = ftell(read_ptr);       // get the current position
            rewind(read_ptr);                 // rewind to the beginning of file
        }
    } else {
        read_ptr = stdin;
    }

    if (vm.count("tx-freq"))
        freq = tx_freq;

    if (vm.count("tx-gain"))
        gain = tx_gain;

    if (vm.count("datarate"))
        rate = datarate;

    if (rate == 0 || freq == 0) {
            printf("Frequency and sample rate need to be specified.\n");
            exit(1);
    }
    
    // if (dual_chan) {
    //     std::string data_filename_2(data_filename);
    //     boost::replace_all(data_filename_2,"chan0","chan1");
    //     printf("Second SigMF Data Filename: %s\n", data_filename_2.c_str());
    //     read_ptr_2 = fopen(data_filename_2.c_str(),"rb");  // r for read, b for binary
    // }

	size_t samps_per_buff = VRT_SAMPLES_PER_PACKET;

	unsigned long long num_requested_samples = total_num_samps;
    double time_requested = total_time;

    uint32_t buffer[VRT_DATA_PACKET_SIZE];
   
    bool first_frame = true;
    bool context_changed = true;

    struct vrt_packet p;
    vrt_init_packet(&p);

    /* Warn if not standards compliant */
    if (vrt_is_platform_little_endian()) {
        printf("Warning: little endian support is work in progress.\n");
    }

    /* VRT init */
    vrt_init_data_packet(&p);
    
    // p.fields.stream_id = stream_id;
    
    // ZMQ
    void *context = zmq_ctx_new();
    void *subscriber = zmq_socket(context, ZMQ_PUB);
    std::string connect_string = "tcp://" + zmq_address + ":" + std::to_string(port);
    int rc = zmq_connect(subscriber, connect_string.c_str());
    assert(rc == 0);

    // stdin binary
    if (read_stdin)
        freopen(NULL, "rb", stdin);
    // _setmode(_fileno(stdin), _O_BINARY);

    // Sleep setup time
    std::this_thread::sleep_for(std::chrono::milliseconds(int64_t(1000 * setup_time)));

    total_num_samps = 0;

    if (total_num_samps == 0) {
        std::signal(SIGINT, &sig_int_handler);
        std::cout << "Press Ctrl + C to stop streaming..." << std::endl;
    }

    // Run this loop until either time expired (if a duration was given), until
    // the requested number of samples were collected (if such a number was
    // given), or until Ctrl-C was pressed.

    uint32_t frame_count = 0;
    uint32_t num_words_read = 0;

    std::complex<short> samples[samps_per_buff];

    timeval time_first_sample;

    boost::posix_time::ptime t1;

    if (read_stdin) {
        // now
        t1 = boost::posix_time::microsec_clock::universal_time();
        timed_tx = false;
    } else {
        // from SigMF
        t1 = boost::posix_time::from_iso_extended_string(start_time_str);
    }

    if (tx_int > 0) {
        timed_tx = true;
        t1 = boost::posix_time::microsec_clock::universal_time();
        printf("    now: %li\n", boost::posix_time::to_time_t(t1));
        t1 = t1 + boost::posix_time::milliseconds(200);
        time_t integer_time_tx = tx_int*(boost::posix_time::to_time_t(t1) / tx_int) + tx_int;
        printf("tx time: %li\n", integer_time_tx);
        t1 = boost::posix_time::from_time_t(integer_time_tx);
    }

    time_t integer_time_first_sample = boost::posix_time::to_time_t(t1);
    boost::posix_time::ptime t2(boost::posix_time::from_time_t(integer_time_first_sample));

    boost::posix_time::time_duration fractional_sec = t1-t2;
    time_first_sample.tv_sec = integer_time_first_sample;
    time_first_sample.tv_usec = fractional_sec.total_microseconds();

    auto vrt_time = time_first_sample;

    int update_interval = 1e6*samps_per_buff/datarate;

    printf("Update interval: %i\n", update_interval);

    unsigned long long last_update_samps = 0;

    if (timed_tx) {
        tx_buffer_size = 0;

        boost::posix_time::time_duration const time_since_epoch=t1-boost::posix_time::from_time_t(0); 
        std::chrono::time_point<std::chrono::system_clock> t_temp = std::chrono::system_clock::from_time_t(time_since_epoch.total_seconds()); 
        long nsec=time_since_epoch.fractional_seconds()*(1000000000/time_since_epoch.ticks_per_second()); 
        auto chrono_t1 = t_temp + std::chrono::nanoseconds(nsec); 

        auto wait_time = chrono_t1 - std::chrono::system_clock::now() - std::chrono::milliseconds(150);

        if (wait_time > std::chrono::microseconds(0))
            std::this_thread::sleep_for(wait_time);
    }

    // time keeping
    auto start_time = std::chrono::steady_clock::now();

    // Track time and samps between updating the BW summary
    auto last_update                     = start_time;
    auto last_context                    = start_time;

    // trigger context update
    last_context -= std::chrono::seconds(4*VRT_CONTEXT_INTERVAL);


    while (not stop_signal_called or last_frame) {
 
        const auto now = std::chrono::steady_clock::now();

        if (frame_count < tx_buffer_size) {
            // don't wait
            start_time = now;
        } else {
            // wait
            auto wait_time = start_time + std::chrono::microseconds((frame_count-tx_buffer_size)*update_interval) - now;
            if (wait_time > std::chrono::microseconds(0))
                std::this_thread::sleep_for(wait_time);
        }

        struct timeval interval_time;
        int64_t first_sample = frame_count*samps_per_buff;

        double interval = (double)first_sample/(double)datarate;
        interval_time.tv_sec = (time_t)interval;
        interval_time.tv_usec = (interval-(time_t)interval)*1e6;

        timeradd(&time_first_sample, &interval_time, &vrt_time);

        const auto time_since_last_context = now - last_context;
        if (last_frame or (send_context and time_since_last_context > std::chrono::milliseconds(VRT_CONTEXT_INTERVAL))) {

            last_context = now;

            // VITA 49.2
            /* Initialize to reasonable values */
            struct vrt_packet pc;
            vrt_init_packet(&pc);

            /* VRT Configure. Note that context packets cannot have a trailer word. */
            vrt_init_context_packet(&pc);

            pc.fields.integer_seconds_timestamp = vrt_time.tv_sec;
            pc.fields.fractional_seconds_timestamp = 1e6*vrt_time.tv_usec;

            pc.fields.stream_id = 1;

            if (freq != 0) {
                pc.if_context.has.rf_reference_frequency = true;
                pc.if_context.rf_reference_frequency            = freq;
                pc.if_context.rf_reference_frequency_offset     = 0;
                pc.if_context.if_reference_frequency            = 0; // Zero-IF
            } else {
                pc.if_context.has.rf_reference_frequency = false;
            }

            if (gain != 0) {
                pc.if_context.has.gain = true;
                pc.if_context.gain.stage1                       = gain;
                pc.if_context.gain.stage2                       = 0;
            } else {
                pc.if_context.has.gain = false;
            }

            pc.if_context.bandwidth                         = bw;
            pc.if_context.sample_rate                       = rate;
            
            if (not context_changed)
                pc.if_context.context_field_change_indicator = false;
            else {
                pc.if_context.context_field_change_indicator = true;
                context_changed = false;
            }

            if (timed_tx) {
                pc.if_context.state_and_event_indicators.has.calibrated_time = true;
                pc.if_context.state_and_event_indicators.calibrated_time = true;
            }

            if (first_frame) {
                pc.if_context.state_and_event_indicators.user_defined = 0x1;
            } else if (last_frame) {
                pc.if_context.state_and_event_indicators.user_defined = 0x2;
            } else {
                pc.if_context.state_and_event_indicators.user_defined = 0x0;
            }

            int32_t rv = vrt_write_packet(&pc, buffer, VRT_DATA_PACKET_SIZE, true);
            if (rv < 0) {
                fprintf(stderr, "Failed to write packet: %s\n", vrt_string_error(rv));
            }
            zmq_send (subscriber, buffer, rv*4, 0);

            // if (dual_chan) {
            //     // duplicate context of channel 0 on channel 1
            //     pc.fields.stream_id = 2;
            //     rv = vrt_write_packet(&pc, buffer, VRT_DATA_PACKET_SIZE, true);
            //     if (rv < 0) {
            //         fprintf(stderr, "Failed to write packet: %s\n", vrt_string_error(rv));
            //     }
            //     zmq_send (subscriber, buffer, rv*4, 0);   
            // }
            last_frame = false;

        }

        // Data
        if (fread(samples, sizeof(samples), 1, read_ptr) == 1) {

            num_words_read = samps_per_buff;

            if (first_frame) {
                std::cout << boost::format(
                                 "First frame: %u samples, %u full secs, %.09f frac secs")
                                 % (num_words_read) % vrt_time.tv_sec
                                 % (vrt_time.tv_usec/1e6)
                          << std::endl;
                first_frame = false;
            }

            if (not repeat && filesize > 0 && ftell(read_ptr) > filesize-sizeof(samples)) {
                last_frame = true;
                // trigger context
                last_context -= std::chrono::seconds(4*VRT_CONTEXT_INTERVAL);
            }

            p.fields.stream_id = 1;
            p.body = samples;
            p.header.packet_count = (uint8_t)frame_count%16;
            p.fields.integer_seconds_timestamp = vrt_time.tv_sec;
            p.fields.fractional_seconds_timestamp = 1e6*vrt_time.tv_usec;

            zmq_msg_t msg;
            int rc = zmq_msg_init_size (&msg, VRT_DATA_PACKET_SIZE*4);

            int32_t rv = vrt_write_packet(&p, zmq_msg_data(&msg), VRT_DATA_PACKET_SIZE, true);

            zmq_msg_send(&msg, subscriber, 0);
            zmq_msg_close(&msg);

            // if (dual_chan) {
            //     if (fread(samples, sizeof(samples), 1, read_ptr_2) == 1) {
            //         p.fields.stream_id = 2;
            //         p.body = samples;
            //         p.header.packet_count = (uint8_t)frame_count%16;
            //         p.fields.integer_seconds_timestamp = vrt_time.tv_sec;
            //         p.fields.fractional_seconds_timestamp = 1e6*vrt_time.tv_usec;
            
            //         zmq_msg_t msg;
            //         int rc = zmq_msg_init_size (&msg, VRT_DATA_PACKET_SIZE*4);

            //         int32_t rv = vrt_write_packet(&p, zmq_msg_data(&msg), VRT_DATA_PACKET_SIZE, true);

            //         zmq_msg_send(&msg, zmq_server, 0);
            //         zmq_msg_close(&msg);
            //     } else {
            //         if (repeat)
            //             rewind(read_ptr_2);
            //         else 
            //             break;
            //     }
            // }

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
            printf("no more samples in data file\n");
            if (not read_stdin and repeat) 
                rewind(read_ptr);
            else
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
