//
// Copyright 2025 by Thomas Telkamp
//
// SPDX-License-Identifier: MIT
//

#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/circular_buffer.hpp>

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

void sig_int_handler(int)
{
    stop_signal_called = true;
}

// Sleep until tp in short slices, so Ctrl+C is not held up by a long wait for
// the start time. Returns false when stopped.
static bool wait_until(std::chrono::steady_clock::time_point tp)
{
    while (not stop_signal_called) {
        const auto now = std::chrono::steady_clock::now();
        if (now >= tp)
            return true;
        std::this_thread::sleep_for(
            std::min<std::chrono::steady_clock::duration>(tp - now, std::chrono::milliseconds(100)));
    }
    return false;
}

int main(int argc, char* argv[])
{
    // variables to be set by po
    std::string ref, file, time_cal, type, start_time_str, zmq_address, start_tx_str;
    size_t tx_int;
    uint16_t port;
    uint32_t stream_id = 0;
    int hwm;
    uint16_t gain = 0, tx_gain;
    double datarate, tx_lead;
    double rate = 0, freq = 0, bw = 0, setup_time, tx_freq, tx_lo_offset;

    FILE *read_ptr;
    FILE *read_ptr_2;

    datarate = 0;
    tx_freq = 0;
    tx_gain = 0;
    tx_lo_offset = 0;

    // setup the program options
    po::options_description desc("Allowed options");
    // clang-format off
    desc.add_options()
        ("help", "help message")
        ("file", po::value<std::string>(&file)->default_value("samples.sigmf-meta"), "name of the SigMF meta file")
        ("setup", po::value<double>(&setup_time)->default_value(1.0), "seconds of setup time")
        ("datarate", po::value<double>(&datarate), "rate of outgoing samples")
        ("tx-lead", po::value<double>(&tx_lead)->default_value(200), "time the stream runs ahead of transmission (ms), must exceed the usrp_to_vrt GPIO delay")
        ("tx-freq", po::value<double>(&tx_freq), "TX RF center frequency in Hz")
        ("tx-gain", po::value<uint16_t>(&tx_gain), "gain for the TX RF chain")
        ("tx-lo-offset", po::value<double>(&tx_lo_offset), "offset for frontend TX LO in Hz")
        // ("dual-chan", "use two SigMF files for dual channel stream (chan0+chan1)")
        ("progress", "periodically display short-term bandwidth")
        ("timed-tx", "Start transmission at given time (SigMF)")
        ("start-time", po::value<std::string>(&start_tx_str), "start transmission at given timestamp")
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
        std::cout << boost::format("SigMF samples to VRT. %s") % desc << std::endl;
        std::cout << std::endl
                  << "This application streams data from a SigMF file "
                     "to VRT.\n"
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
    bool start_at_timestamp     = vm.count("start-time") > 0;
    bool send_context           = true;

    if (tx_lead <= 0) {
        printf("--tx-lead needs to be positive.\n");
        exit(1);
    }

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
        printf("       Frequency: %lli\n", (long long int)freq);
        printf("       Data type: %s\n", type.c_str());
        printf("       Stream ID: %u\n", stream_id);

        // Some Checks

        if (type != "ci16_le") {
            printf("Only 16 bit complex int data format supported (\"ci16_le\")\n");
            exit(1);
        }

        // Open data file

        std::string data_filename;
        base_fn_fp.replace_extension(".sigmf-data");
        data_filename = base_fn_fp.string();

        printf("SigMF Data Filename: %s\n", data_filename.c_str());

        read_ptr = fopen(data_filename.c_str(),"rb");  // r for read, b for binary
        if (read_ptr == nullptr) {
            perror(data_filename.c_str());
            exit(1);
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
    else
        datarate = rate;

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

    const size_t samps_per_buff = VRT_SAMPLES_PER_PACKET;

    // header, stream ID, class ID and timestamp words in front of the samples
    const uint32_t packet_overhead = VRT_DATA_PACKET_SIZE - VRT_SAMPLES_PER_PACKET;

    uint32_t buffer[VRT_DATA_PACKET_SIZE];

    struct vrt_packet p;
    vrt_init_packet(&p);

    /* Warn if not standards compliant */
    if (vrt_is_platform_little_endian()) {
        printf("Warning: little endian support is work in progress.\n");
    }

    /* VRT init */
    vrt_init_data_packet(&p);

    // p.fields.stream_id = stream_id;

    // The lead is queued in ZMQ while the transmitter waits for its start time,
    // a PUB socket drops what does not fit
    const double lead_packets = tx_lead / 1000.0 * datarate / samps_per_buff;
    printf("TX lead: %.0f ms (%.0f packets)\n", tx_lead, lead_packets);
    if (lead_packets > hwm) {
        printf("Warning: TX lead exceeds the ZMQ HWM of %i packets, samples will be dropped.\n", hwm);
    }

    // ZMQ
    void *context = zmq_ctx_new();
    void *publisher = zmq_socket(context, ZMQ_PUB);
    int rc = zmq_setsockopt(publisher, ZMQ_SNDHWM, &hwm, sizeof hwm);
    assert(rc == 0);
    // deliver what is still queued on exit, including the end of transmit
    int linger = 2000;
    rc = zmq_setsockopt(publisher, ZMQ_LINGER, &linger, sizeof linger);
    assert(rc == 0);
    std::string connect_string = "tcp://" + zmq_address + ":" + std::to_string(port);
    rc = zmq_connect(publisher, connect_string.c_str());
    assert(rc == 0);

    // stdin binary
    if (read_stdin)
        freopen(NULL, "rb", stdin);
    // _setmode(_fileno(stdin), _O_BINARY);

    std::signal(SIGINT, &sig_int_handler);

    // Sleep setup time
    std::this_thread::sleep_for(std::chrono::milliseconds(int64_t(1000 * setup_time)));

    std::cout << "Press Ctrl + C to stop streaming..." << std::endl;

    std::complex<short> samples[VRT_SAMPLES_PER_PACKET];

    struct vrt_time_ps t1 = {0, 0};

    if (read_stdin) {
        // now
        t1 = vrt_time_now();
        timed_tx = false;
    } else {
        // from SigMF
        if (not vrt_parse_iso_datetime_ns(start_time_str, &t1)) {
            std::cerr << "Failed to parse core:datetime: " << start_time_str << std::endl;
            exit(1);
        }
    }

    if (tx_int > 0) {
        timed_tx = true;
        t1 = vrt_time_now();
        printf("    now: %li\n", (long int)t1.seconds);
        // leave room for the full lead
        t1.frac_ps += (uint64_t)((tx_lead + 100) * 1e9);
        vrt_time_normalize(&t1);
        time_t integer_time_tx = tx_int*(t1.seconds / tx_int) + tx_int;
        printf("tx time: %li\n", integer_time_tx);
        t1.seconds = integer_time_tx;
        t1.frac_ps = 0;
    }

    if (start_at_timestamp) {

        // Unix time or ISO 8601, kept to picosecond resolution
        struct vrt_time_ps utc_time = {0, 0};
        if (not vrt_parse_time_arg(start_tx_str, &utc_time)) {
            std::cerr << "Failed to parse --start-time: " << start_tx_str << std::endl;
            exit(1);
        }
        // Print parsed time
        std::cout << "UTC start time: " << vrt_iso_datetime_ns(utc_time.seconds, utc_time.frac_ps)
                  << std::endl;

        t1 = utc_time;
        timed_tx = true;
    }

    const struct vrt_time_ps time_first_sample = t1;

    // Pacing: the samples at stream time s are sent at anchor + s. In a timed
    // transmission the anchor is the lead before the start time. Otherwise the
    // transmission starts when the first samples arrive, and the lead is sent
    // right away: the anchor is set at the first read, stdin may be slow to start.
    const auto lead = std::chrono::duration_cast<std::chrono::steady_clock::duration>(
        std::chrono::duration<double>(tx_lead / 1000.0));
    std::chrono::steady_clock::time_point anchor;

    if (timed_tx) {
        const struct vrt_time_ps now_ps = vrt_time_now();
        const double time_to_start = (double)(t1.seconds - now_ps.seconds)
                                     + ((double)t1.frac_ps - (double)now_ps.frac_ps) / 1e12;
        if (time_to_start <= 0) {
            printf("Start time in the past\n");
            exit(1);
        }
        if (time_to_start < tx_lead / 1000.0) {
            printf("Warning: start time is only %.0f ms away, less than the TX lead.\n",
                   time_to_start * 1000.0);
        }
        anchor = std::chrono::steady_clock::now()
                 + std::chrono::duration_cast<std::chrono::steady_clock::duration>(
                     std::chrono::duration<double>(time_to_start))
                 - lead;
    }

    // time keeping
    auto start_time = std::chrono::steady_clock::now();

    // Track time and samps between updating the BW summary
    auto last_update  = start_time;
    auto last_context = start_time;

    unsigned long long last_update_samps = 0;

    uint32_t frame_count = 0;
    uint64_t samples_sent = 0;
    uint64_t late_packets = 0;
    bool first_frame = true;
    bool context_changed = true;
    bool read_since_rewind = false;

    // VITA 49.2 context, timestamped with the time of the next sample
    auto send_context_packet = [&](uint8_t user_defined) {
        const struct vrt_time_ps vrt_time =
            vrt_time_add_samples(time_first_sample, samples_sent, datarate);

        /* Initialize to reasonable values */
        struct vrt_packet pc;
        vrt_init_packet(&pc);

        /* VRT Configure. Note that context packets cannot have a trailer word. */
        vrt_init_context_packet(&pc);

        pc.fields.integer_seconds_timestamp = vrt_time.seconds;
        pc.fields.fractional_seconds_timestamp = vrt_time.frac_ps;

        pc.fields.stream_id = 1;

        if (freq != 0) {
            pc.if_context.has.rf_reference_frequency = true;
            pc.if_context.rf_reference_frequency            = freq;
            pc.if_context.rf_reference_frequency_offset     = 0;
            pc.if_context.if_reference_frequency            = 0; // Zero-IF
            if (vm.count("tx-lo-offset")) {
                pc.if_context.has.if_band_offset = true;
                pc.if_context.if_band_offset = tx_lo_offset;
            } else {
                pc.if_context.has.if_band_offset = false;
            }
        } else {
            pc.if_context.has.rf_reference_frequency = false;
            pc.if_context.has.if_band_offset = false;
        }

        pc.if_context.has.gain = true;
        pc.if_context.gain.stage1                       = gain;
        pc.if_context.gain.stage2                       = 0;

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

        // 0x1 start of transmission, 0x2 end of transmission
        pc.if_context.state_and_event_indicators.user_defined = user_defined;

        int32_t rv = vrt_write_packet(&pc, buffer, VRT_DATA_PACKET_SIZE, true);
        if (rv < 0) {
            fprintf(stderr, "Failed to write packet: %s\n", vrt_string_error(rv));
            return;
        }
        zmq_send (publisher, buffer, rv*4, 0);

        // if (dual_chan) {
        //     // duplicate context of channel 0 on channel 1
        //     pc.fields.stream_id = 2;
        //     rv = vrt_write_packet(&pc, buffer, VRT_DATA_PACKET_SIZE, true);
        //     if (rv < 0) {
        //         fprintf(stderr, "Failed to write packet: %s\n", vrt_string_error(rv));
        //     }
        //     zmq_send (publisher, buffer, rv*4, 0);
        // }
    };

    // Run this loop until the input ends or Ctrl-C was pressed. Every way out
    // of the loop is followed by the end of transmission context.
    while (not stop_signal_called) {

        // Data, a short read is the tail of the input
        size_t num_words_read = fread(samples, sizeof(samples[0]), samps_per_buff, read_ptr);

        if (num_words_read == 0) {
            if (ferror(read_ptr)) {
                perror("Failed to read samples");
                break;
            }
            printf("no more samples in data file\n");
            if (not read_stdin and repeat and read_since_rewind) {
                rewind(read_ptr);
                read_since_rewind = false;
                continue;
            }
            break;
        }
        read_since_rewind = true;

        if (first_frame and not timed_tx)
            anchor = std::chrono::steady_clock::now() - lead;

        // wait for the send time of these samples
        const auto send_time = anchor
                               + std::chrono::duration_cast<std::chrono::steady_clock::duration>(
                                   std::chrono::duration<double>((double)samples_sent / datarate));
        if (not wait_until(send_time))
            break;

        const auto now = std::chrono::steady_clock::now();

        // after the lead has passed, these samples are due for transmission already
        if (not first_frame and now > send_time + lead) {
            if (late_packets == 0) {
                printf("Warning: stream is behind by %.1f ms, transmitter may underflow.\n",
                       std::chrono::duration<double, std::milli>(now - send_time - lead).count());
            }
            late_packets++;
        }

        if (first_frame) {
            start_time = now;
            last_update = now;
            last_context = now;
            send_context_packet(0x1);
        } else if (send_context and now - last_context > std::chrono::milliseconds(VRT_CONTEXT_INTERVAL)) {
            last_context = now;
            send_context_packet(0x0);
        }

        const struct vrt_time_ps vrt_time =
            vrt_time_add_samples(time_first_sample, samples_sent, datarate);

        if (first_frame) {
            std::cout << boost::format(
                             "First frame: %u samples, %u full secs, %.09f frac secs")
                             % (num_words_read) % vrt_time.seconds
                             % ((double)vrt_time.frac_ps/1e12)
                      << std::endl;
            first_frame = false;
        }

        p.fields.stream_id = 1;
        p.body = samples;
        p.words_body = num_words_read;
        p.header.packet_size = num_words_read + packet_overhead;
        p.header.packet_count = (uint8_t)frame_count%16;
        p.fields.integer_seconds_timestamp = vrt_time.seconds;
        p.fields.fractional_seconds_timestamp = vrt_time.frac_ps;

        zmq_msg_t msg;
        zmq_msg_init_size (&msg, p.header.packet_size*4);

        int32_t rv = vrt_write_packet(&p, zmq_msg_data(&msg), p.header.packet_size, true);
        if (rv < 0) {
            fprintf(stderr, "Failed to write packet: %s\n", vrt_string_error(rv));
            zmq_msg_close(&msg);
            break;
        }

        zmq_msg_send(&msg, publisher, 0);
        zmq_msg_close(&msg);

        // if (dual_chan) {
        //     if (fread(samples, sizeof(samples), 1, read_ptr_2) == 1) {
        //         p.fields.stream_id = 2;
        //         p.body = samples;
        //         p.header.packet_count = (uint8_t)frame_count%16;
        //         p.fields.integer_seconds_timestamp = vrt_time.seconds;
        //         p.fields.fractional_seconds_timestamp = vrt_time.frac_ps;

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
        samples_sent += num_words_read;
        num_total_samps += num_words_read;

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

                for (size_t i=0; i<num_words_read; i++ ) {
                    auto sample_i = get_abs_val(samples[i]);
                    sum_i += sample_i;
                    if (sample_i > datatype_max*0.99)
                        clip_i++;
                }
                sum_i = sum_i/num_words_read;
                std::cout << boost::format("%.0f") % (100.0*log2(sum_i)/log2(datatype_max)) << "% I (";
                std::cout << boost::format("%.0f") % ceil(log2(sum_i)+1) << " of ";
                std::cout << (int)ceil(log2(datatype_max)+1) << " bits), ";
                std::cout << "" << boost::format("%.0f") % (100.0*clip_i/num_words_read) << "% I clip.";
                std::cout << std::endl;

            }
        }
    }

    // End of transmission, timestamped with the end of the last sample sent.
    // Nothing is sent after it. Without any samples sent there is nothing to end.
    if (not first_frame) {
        send_context_packet(0x2);
        printf("End of transmission after %llu samples.\n", (unsigned long long)samples_sent);
    }

    const auto actual_stop_time = std::chrono::steady_clock::now();

    if (late_packets > 0) {
        printf("Warning: %llu packets were sent later than the TX lead allows.\n",
               (unsigned long long)late_packets);
    }

    if (stats) {
        std::cout << std::endl;
        const double actual_duration_seconds =
            std::chrono::duration<float>(actual_stop_time - start_time).count();

        std::cout << boost::format("Sent %d samples in %f seconds.") % num_total_samps
                         % actual_duration_seconds
                  << std::endl;
        const double rate = (double)num_total_samps / actual_duration_seconds;
        std::cout << (rate / 1e6) << " Msps." << std::endl;
    }

    /* clean up */
    fclose(read_ptr);

    // blocks until the queue is delivered, or the linger time expired
    zmq_close(publisher);
    zmq_ctx_term(context);

    // finished
    std::cout << std::endl << "Done!" << std::endl << std::endl;

    return EXIT_SUCCESS;
}
