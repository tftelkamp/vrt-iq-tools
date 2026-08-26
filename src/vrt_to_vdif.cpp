#include <zmq.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>

#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/thread/thread.hpp>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/posix_time_io.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <csignal>
#include <fstream>
#include <iostream>
#include <thread>
#include <vector>

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
#include "vdif.h"
#include "dt-extended-context.h"
#include "tracker-extended-context.h"

/* Fallback values in case the defines are not provided */
#ifndef GIT_COMMIT
#define GIT_COMMIT "unknown"
#endif

namespace po = boost::program_options;

/* Interval between progress lines, in seconds. */
#define VDIF_STATS_INTERVAL 1.0

static bool stop_signal_called = false;
void sig_int_handler(int)
{
    stop_signal_called = true;
}

/* Return whether a .vdif of this base name exists */
bool vdif_data_file_exists(std::string base_filename)
{
    return boost::filesystem::exists(boost::filesystem::path(base_filename + ".vdif"));
}

/* Get suffix _1, _2 so that base_filename_suffix does not overwrite existing data */
std::string generate_nonexisting_base_filename_suffix(std::string base_filename)
{
    if (not vdif_data_file_exists(base_filename)) {
        return "";
    } else {
        for (int i = 1;; i++) {
            if (not vdif_data_file_exists(base_filename + "_" + std::to_string(i)))
                return "_" + std::to_string(i);
        }
    }
}

int main(int argc, char* argv[])
{

    // variables to be set by po
    std::string file, auto_file, zmq_address, station, start_reception;
    uint16_t    instance, main_port, port, thread_id;
    uint32_t    channel, bits_option, frame_bytes_option;
    int         hwm;
    size_t      num_requested_samples;
    double      total_time, max_fill, jitter;

    // setup the program options
    po::options_description desc("Allowed options");
    // clang-format off

    desc.add_options()
        ("help", "help message")
        ("file", po::value<std::string>(&file)->default_value("vrt_samples"), "name of the file to write VDIF frames to")
        ("auto-file", po::value<std::string>(&auto_file), "prefix of the auto generated filename to write VDIF frames to")
        ("nsamps", po::value<size_t>(&num_requested_samples)->default_value(0), "total number of samples to receive")
        ("duration", po::value<double>(&total_time)->default_value(0), "total number of seconds to receive")
        ("channel", po::value<uint32_t>(&channel)->default_value(0), "VRT channel")
        ("station", po::value<std::string>(&station)->default_value("0"), "VDIF station ID: number or two characters")
        ("thread", po::value<uint16_t>(&thread_id)->default_value(0), "VDIF thread ID")
        ("frame-bytes", po::value<uint32_t>(&frame_bytes_option)->default_value(VDIF_DEFAULT_PAYLOAD_BYTES), "VDIF data array size [bytes]")
        ("bits", po::value<uint32_t>(&bits_option)->default_value(2), "bits per component, if the context does not signal it")
        ("max-fill", po::value<double>(&max_fill)->default_value(1.0), "gap that is filled with invalid frames before resyncing [seconds]")
        ("jitter", po::value<double>(&jitter)->default_value(1.0), "timestamp jitter that is taken up by the sample count [microseconds]")
        ("int-second", "align start of reception to integer second")
        ("start-time", po::value<std::string>(&start_reception), "start reception at given timestamp")
        ("progress", "periodically display short-term bandwidth")
        ("null", "run without writing to file")
        ("continue", "don't abort on a bad packet")
        ("dt-trace", "add DT trace data")
        ("tracking", "add tracking context data")
        ("address", po::value<std::string>(&zmq_address)->default_value("localhost"), "VRT ZMQ address")
        ("instance", po::value<uint16_t>(&instance)->default_value(0), "VRT ZMQ instance")
        ("port", po::value<uint16_t>(&port), "VRT ZMQ port")
        ("hwm", po::value<int>(&hwm)->default_value(10000), "VRT ZMQ HWM")
    ;
    // clang-format on
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // print the help message
    if (vm.count("help")) {
        std::cout << boost::format("VRT samples to VDIF. %s") % desc << std::endl;
        std::cout << std::endl
                  << "This application writes a single channel of a VRT stream to a file of VDIF\n"
                     "(VLBI Data Interchange Format, release 1.1.1) data frames, as one complex\n"
                     "single-channel data thread with 32 byte headers and EDV 0.\n"
                     "\n"
                     "Streams of 1 or 2 bits per component, as produced by vrt_quantize, and plain\n"
                     "16 bit streams are supported. The sample size is taken from the payload format\n"
                     "of the first context packet, falling back on --bits for streams that do not\n"
                     "signal it.\n"
                     "\n"
                     "The data array size is chosen such that there is an integral number of data\n"
                     "frames per second, taking the largest size that does not exceed --frame-bytes,\n"
                     "or else the smallest one above it. Recording starts at the first data frame\n"
                     "boundary of the stream, which may be in the middle of a second, or with\n"
                     "--int-second at the start of the next second. A timestamp that falls between\n"
                     "two samples is rounded to the nearest one, and the remainder is reported as\n"
                     "start_residual_seconds in the metadata.\n"
                     "\n"
                     "Gaps in the stream are filled with frames flagged invalid, up to --max-fill\n"
                     "seconds, beyond which the frame count is resynchronised to the incoming\n"
                     "timestamps. Timestamps that are off by no more than --jitter are taken up by\n"
                     "the sample count, so that a source timestamping in whole microseconds does not\n"
                     "fragment the recording.\n"
                     "\n"
                     "The metadata of the VRT stream that has no place in a VDIF header, such as the\n"
                     "RF frequency and the sample rate, is written to a separate .vdif-meta file.\n"
                  << std::endl;
        return ~0;
    }

    bool progress               = vm.count("progress") > 0;
    bool null                   = vm.count("null") > 0;
    bool continue_on_bad_packet = vm.count("continue") > 0;
    bool dt_trace               = vm.count("dt-trace") > 0;
    bool tracking               = vm.count("tracking") > 0;
    bool do_auto_file           = vm.count("auto-file") > 0;
    bool int_second             = vm.count("int-second") > 0;
    bool start_at_timestamp     = vm.count("start-time") > 0;
    bool frame_bytes_given      = not vm["frame-bytes"].defaulted();

    uint16_t station_id;
    bool     station_is_ascii;
    if (not vdif_parse_station_id(station, &station_id, &station_is_ascii)) {
        std::cerr << "Station ID \"" << station
                  << "\" is neither a 16-bit number nor two characters." << std::endl;
        return ~0;
    }

    if (thread_id > 1023) {
        std::cerr << "Thread ID must be in the range 0 to 1023." << std::endl;
        return ~0;
    }

    if (max_fill < 0) {
        std::cerr << "--max-fill cannot be negative." << std::endl;
        return ~0;
    }

    if (jitter < 0) {
        std::cerr << "--jitter cannot be negative." << std::endl;
        return ~0;
    }

    int64_t start_second_requested = 0;
    if (start_at_timestamp) {
        // Unix time or ISO 8601, kept to picosecond resolution
        struct vrt_time_ps utc_time = {0, 0};
        if (not vrt_parse_time_arg(start_reception, &utc_time)) {
            std::cerr << "Failed to parse --start-time: " << start_reception << std::endl;
            return ~0;
        }
        std::cout << "UTC start time: " << vrt_iso_datetime_ns(utc_time.seconds, utc_time.frac_ps)
                  << std::endl;
        /* VDIF frames start on a second boundary, so the fraction only selects
         * which second is waited for. */
        start_second_requested = utc_time.seconds;
    }

    context_type              vrt_context;
    dt_ext_context_type       dt_ext_context;
    tracker_ext_context_type  tracker_ext_context;
    init_context(&vrt_context);

    packet_type vrt_packet;
    vrt_packet.first_frame = true;
    vrt_packet.lost_frame  = false;

    if (vm.count("port") > 0) {
        main_port = port;
    } else {
        main_port = DEFAULT_MAIN_PORT + MAX_CHANNELS * instance;
    }

    vrt_packet.channel_filt = 1 << channel;

    // ZMQ
    void* context    = zmq_ctx_new();
    void* subscriber = zmq_socket(context, ZMQ_SUB);
    int   rc         = zmq_setsockopt(subscriber, ZMQ_RCVHWM, &hwm, sizeof hwm);
    std::string connect_string = "tcp://" + zmq_address + ":" + std::to_string(main_port);
    rc                         = zmq_connect(subscriber, connect_string.c_str());
    assert(rc == 0);
    zmq_setsockopt(subscriber, ZMQ_SUBSCRIBE, "", 0);

    std::signal(SIGINT, &sig_int_handler);

    // time keeping
    auto start_time = std::chrono::steady_clock::now();

    uint32_t buffer[ZMQ_BUFFER_SIZE];

    // Files
    file += generate_nonexisting_base_filename_suffix(file);
    std::string data_filename = file + ".vdif";
    std::string meta_filename = file + ".vdif-meta";
    std::ofstream datafile;

    // Frame geometry, known once the first context packet has been received
    uint32_t bits              = 0;
    uint32_t bits_per_sample   = 0;
    uint32_t sample_rate       = 0;
    uint32_t payload_bytes     = 0;
    uint32_t payload_words     = 0;
    uint32_t samples_per_frame = 0;
    uint32_t frames_per_second = 0;

    struct vdif_header header;
    uint8_t            reference_epoch = 0;
    int64_t            epoch_start     = 0;

    // Frame under construction
    std::vector<uint32_t> frame;
    uint64_t              frame_fill    = 0;
    bool                  frame_invalid = false;

    // Position of the next sample to be written, in samples since 1970
    bool     started    = false;
    bool     have_start = false;
    uint64_t abs_pos    = 0;
    uint64_t start_abs  = 0;

    // Rounding of the first timestamp onto the sample grid, in seconds
    double start_residual = 0;

    // Counters for the summary
    uint64_t frames_written    = 0;
    uint64_t frames_invalid    = 0;
    uint64_t num_total_samps   = 0;
    uint64_t fill_samps        = 0;
    uint64_t late_packets      = 0;
    uint64_t lost_vrt_frames   = 0;
    uint64_t resyncs           = 0;
    uint64_t jitter_packets    = 0;
    uint64_t context_changes   = 0;

    // Timestamp error that is taken as jitter rather than as missing data
    uint64_t jitter_samples = 0;

    // State occupancy of the quantizer, as a histogram over payload bytes so
    // that it costs a single increment per byte
    uint64_t byte_hist[256]       = {0};
    uint64_t byte_hist_total[256] = {0};

    // Level statistics of a 16 bit stream, index 0 is the real component
    double   sample_peak[2]  = {0, 0};
    double   sample_sumsq[2] = {0, 0};
    uint64_t sample_count    = 0;

    // Context values at the start of the recording, to detect changes
    int64_t  start_rf_freq   = 0;
    int32_t  start_gain      = 0;
    uint32_t start_bandwidth = 0;

    auto last_stats        = start_time;
    uint64_t last_stats_frames = 0;

    // set to true to process data before context
    bool start_rx = false;

    // set when the recording is given up on, to be told apart from a clean end
    bool fatal = false;

    int len;

    /* Expand the byte histogram into per-component state counts. One byte holds
     * an even number of components, alternating between real and imaginary. */
    auto expand_hist = [&](const uint64_t* hist, uint64_t counts[2][4]) {
        const unsigned fields = 8 / bits;
        const uint32_t mask   = (1u << bits) - 1;
        memset(counts, 0, 2 * 4 * sizeof(uint64_t));
        for (unsigned b = 0; b < 256; b++) {
            if (hist[b] == 0)
                continue;
            for (unsigned j = 0; j < fields; j++)
                counts[j % 2][(b >> (j * bits)) & mask] += hist[b];
        }
    };

    /* Write the frame under construction, using the timestamp of its first
     * sample. Data frames never straddle a second, since the number of samples
     * per frame divides the sample rate. A frame that is only partly filled is
     * padded with zeros and flagged invalid by the caller. */
    auto flush_frame = [&]() {
        const uint64_t start   = abs_pos - frame_fill;
        const uint64_t second  = start / sample_rate;
        const uint32_t frame_n = (uint32_t)((start % sample_rate) / samples_per_frame);

        vdif_set_time(&header, (uint32_t)((int64_t)second - epoch_start), frame_n, frame_invalid);

        if (not null) {
            datafile.write((const char*)header.word, VDIF_HEADER_BYTES);
            datafile.write((const char*)frame.data(), payload_bytes);
        }

        frames_written++;
        if (frame_invalid)
            frames_invalid++;

        memset(frame.data(), 0, frame.size() * sizeof(uint32_t));
        frame_fill    = 0;
        frame_invalid = false;
    };

    /* Add n samples to the stream, taken from sample src_sample of payload, or
     * zeros with the invalid flag set when payload is NULL. */
    auto append_samples = [&](const uint32_t* payload, uint64_t src_sample, uint64_t n) {
        while (n > 0) {
            const uint64_t take = std::min((uint64_t)samples_per_frame - frame_fill, n);

            if (payload == NULL) {
                /* The frame buffer is already zeroed */
                frame_invalid = true;
            } else if (bits == 16) {
                /* A 16 bit complex sample is exactly one word, so this is a
                 * straight copy apart from the conversion of both components
                 * from two's complement to the offset binary of VDIF */
                for (uint64_t i = 0; i < take; i++)
                    frame[frame_fill + i] = payload[src_sample + i] ^ 0x80008000u;
            } else {
                vdif_copy_bits(frame.data(), frame_fill * bits_per_sample, payload,
                               src_sample * bits_per_sample, take * bits_per_sample);
            }

            frame_fill += take;
            src_sample += take;
            abs_pos += take;
            n -= take;

            if (frame_fill == samples_per_frame)
                flush_frame();
        }
    };

    /* Everything of the VRT stream that a VDIF header has no place for, plus the
     * VDIF parameters in use, so that the recording is self describing. Written
     * once at the start of the recording and rewritten with the summary added
     * when it ends. */
    auto write_meta = [&](bool final) {
        if (null)
            return;

        std::ofstream metafile(meta_filename.c_str());
        if (not metafile.is_open()) {
            std::cerr << "# Error: could not write " << meta_filename << std::endl;
            return;
        }

        const uint64_t start_second = start_abs / sample_rate;

        std::string json = str(boost::format("{\n"
            "    \"vdif\": {\n"
            "        \"recorder\": \"vrt_to_vdif\",\n"
            "        \"git_commit\": \"%s\",\n"
            "        \"version\": %u,\n"
            "        \"edv\": 0,\n"
            "        \"header_bytes\": %u,\n"
            "        \"frame_bytes\": %u,\n"
            "        \"payload_bytes\": %u,\n"
            "        \"samples_per_frame\": %u,\n"
            "        \"frames_per_second\": %u,\n"
            "        \"sample_rate\": %u,\n"
            "        \"complex\": true,\n"
            "        \"bits_per_component\": %u,\n"
            "        \"bits_per_sample\": %u,\n"
            "        \"log2_channels\": 0,\n"
            "        \"thread_id\": %u,\n"
            "        \"station_id\": %u,\n")
            % GIT_COMMIT
            % VDIF_VERSION
            % VDIF_HEADER_BYTES
            % (payload_bytes + VDIF_HEADER_BYTES)
            % payload_bytes
            % samples_per_frame
            % frames_per_second
            % sample_rate
            % bits
            % bits_per_sample
            % thread_id
            % station_id);

        if (station_is_ascii)
            json += str(boost::format("        \"station_id_ascii\": \"%c%c\",\n")
                        % (char)(station_id >> 8) % (char)(station_id & 0xFF));

        json += str(boost::format(
            "        \"reference_epoch\": %u,\n"
            "        \"reference_epoch_datetime\": \"%s\",\n"
            "        \"start_seconds_from_epoch\": %lld,\n"
            "        \"start_datetime\": \"%s\",\n"
            "        \"start_frame\": %u,\n"
            "        \"start_residual_seconds\": %.6e,\n"
            "        \"data_file\": \"%s\"\n"
            "    },\n")
            % (unsigned)reference_epoch
            % (boost::posix_time::to_iso_extended_string(boost::posix_time::from_time_t(epoch_start)))
            % (long long int)((int64_t)start_second - epoch_start)
            % (boost::posix_time::to_iso_extended_string(boost::posix_time::from_time_t(start_second)))
            % (uint32_t)((start_abs % sample_rate) / samples_per_frame)
            % start_residual
            % data_filename);

        json += str(boost::format("    \"vrt\": {\n"
            "        \"sample_rate\": %u,\n"
            "        \"frequency\": %.0f,\n"
            "        \"frac_frequency\": %.6e,\n"
            "        \"bandwidth\": %u,\n"
            "        \"rx_gain\": %i,\n"
            "        \"reference\": \"%s\",\n"
            "        \"time_source\": \"%s\",\n"
            "        \"time_adjust\": %.9f,\n"
            "        \"cal_time\": %u,\n"
            "        \"stream_id\": %u,\n"
            "        \"channel\": %u\n"
            "    }")
            % vrt_context.sample_rate
            % (double)vrt_context.rf_freq
            % vrt_context.rf_frac_freq
            % vrt_context.bandwidth
            % vrt_context.gain
            % (vrt_context.reflock ? "external" : "internal")
            % (vrt_context.time_cal ? "pps" : "internal")
            % ((double)vrt_context.timestamp_adjustment / 1e12)
            % vrt_context.timestamp_calibration_time
            % vrt_context.stream_id
            % channel);

        if (dt_trace and dt_ext_context.dt_ext_context_received) {
            char const* trackerStrings[] = {"idle",       "azel",       "j2000tracker", "moontracker",
                                            "suntracker", "sattracker", "manual"};
            const uint32_t tracker = (dt_ext_context.active_tracker <
                                      sizeof(trackerStrings) / sizeof(trackerStrings[0]))
                                         ? dt_ext_context.active_tracker
                                         : 0;
            json += str(boost::format(",\n    \"dt\": {\n"
                "        \"datetime\": \"%s\",\n"
                "        \"pointing:active_tracker\": \"%s\",\n"
                "        \"pointing:tracking_enabled\": %s,\n"
                "        \"pointing:refraction\": %s,\n"
                "        \"pointing:dt_model\": %s,\n"
                "        \"pointing:refraction_j2000\": %s,\n"
                "        \"pointing:dt_model_j2000\": %s,\n"
                "        \"pointing:current:az_deg\": %.3f,\n"
                "        \"pointing:current:el_deg\": %.3f,\n"
                "        \"pointing:error:az_deg\": %.3f,\n"
                "        \"pointing:error:el_deg\": %.3f,\n"
                "        \"pointing:offset:az_deg\": %.3f,\n"
                "        \"pointing:offset:el_deg\": %.3f,\n"
                "        \"pointing:setpoint:ra_h\": %.3f,\n"
                "        \"pointing:setpoint:dec_deg\": %.3f,\n"
                "        \"pointing:current:ra_h\": %.3f,\n"
                "        \"pointing:current:dec_deg\": %.3f,\n"
                "        \"focusbox_position_mm\": %.0f\n"
                "    }")
                % (vrt_iso_datetime_ns(dt_ext_context.integer_seconds_timestamp,
                                       dt_ext_context.fractional_seconds_timestamp))
                % (trackerStrings[tracker])
                % (dt_ext_context.tracking_enabled ? "true" : "false")
                % (dt_ext_context.refraction ? "true" : "false")
                % (dt_ext_context.dt_model ? "true" : "false")
                % (dt_ext_context.refraction_j2000 ? "true" : "false")
                % (dt_ext_context.dt_model_j2000 ? "true" : "false")
                % ((180.0 / M_PI) * dt_ext_context.azimuth)
                % ((180.0 / M_PI) * dt_ext_context.elevation)
                % ((180.0 / M_PI) * dt_ext_context.azimuth_error)
                % ((180.0 / M_PI) * dt_ext_context.elevation_error)
                % ((180.0 / M_PI) * dt_ext_context.azimuth_offset)
                % ((180.0 / M_PI) * dt_ext_context.elevation_offset)
                % ((12.0 / M_PI) * dt_ext_context.ra_setpoint)
                % ((180.0 / M_PI) * dt_ext_context.dec_setpoint)
                % ((12.0 / M_PI) * dt_ext_context.ra_current)
                % ((180.0 / M_PI) * dt_ext_context.dec_current)
                % dt_ext_context.focusbox);
        }

        if (tracking and tracker_ext_context.tracker_ext_context_received) {
            json += str(boost::format(",\n    \"tracker\": {\n"
                "        \"datetime\": \"%s\",\n"
                "        \"object_name\": \"%.32s\",\n"
                "        \"tracking_source\": \"%.32s\",\n"
                "        \"object_id\": %i,\n"
                "        \"az_deg\": %.3f,\n"
                "        \"el_deg\": %.3f,\n"
                "        \"ra_h\": %.3f,\n"
                "        \"dec_deg\": %.3f,\n"
                "        \"distance\": %.2f,\n"
                "        \"speed\": %.2f,\n"
                "        \"frequency\": %.0f,\n"
                "        \"doppler\": %.4f,\n"
                "        \"doppler_rate\": %.4f\n"
                "    }")
                % (vrt_iso_datetime_ns(tracker_ext_context.integer_seconds_timestamp,
                                       tracker_ext_context.fractional_seconds_timestamp))
                % (tracker_ext_context.object_name)
                % (tracker_ext_context.tracking_source)
                % (tracker_ext_context.object_id)
                % (tracker_ext_context.azimuth)
                % (tracker_ext_context.elevation)
                % ((12.0 / 180.0) * tracker_ext_context.ra)
                % (tracker_ext_context.dec)
                % (tracker_ext_context.distance)
                % (tracker_ext_context.speed)
                % (tracker_ext_context.frequency)
                % (tracker_ext_context.doppler)
                % (tracker_ext_context.doppler_rate));
        }

        if (final) {
            const uint64_t end_second = abs_pos / sample_rate;

            json += str(boost::format(",\n    \"summary\": {\n"
                "        \"end_datetime\": \"%s\",\n"
                "        \"duration_seconds\": %.6f,\n"
                "        \"frames_written\": %llu,\n"
                "        \"frames_invalid\": %llu,\n"
                "        \"samples_written\": %llu,\n"
                "        \"fill_samples\": %llu,\n"
                "        \"vrt_frames_lost\": %llu,\n"
                "        \"vrt_packets_late\": %llu,\n"
                "        \"resyncs\": %llu,\n"
                "        \"jitter_packets\": %llu,\n"
                "        \"context_changes\": %llu")
                % (vrt_iso_datetime_ns(end_second,
                                       (uint64_t)(((unsigned __int128)(abs_pos % sample_rate) *
                                                   VRT_PS_PER_SECOND) / sample_rate)))
                % ((double)(abs_pos - start_abs) / sample_rate)
                % (unsigned long long)frames_written
                % (unsigned long long)frames_invalid
                % (unsigned long long)num_total_samps
                % (unsigned long long)fill_samps
                % (unsigned long long)lost_vrt_frames
                % (unsigned long long)late_packets
                % (unsigned long long)resyncs
                % (unsigned long long)jitter_packets
                % (unsigned long long)context_changes);

            /* Occupancy over the data received, which is what says whether the
             * quantizer was set up well: near 50/50 for 1 bit and near
             * 16/34/34/16 for 2 bit */
            if (bits <= 2) {
                uint64_t counts[2][4];
                expand_hist(byte_hist_total, counts);
                for (int c = 0; c < 2; c++) {
                    uint64_t total = 0;
                    for (uint32_t s = 0; s < (1u << bits); s++)
                        total += counts[c][s];
                    json += str(boost::format(",\n        \"state_occupancy_percent_%s\": [")
                                % (c == 0 ? "i" : "q"));
                    for (uint32_t s = 0; s < (1u << bits); s++)
                        json += str(boost::format("%s%.3f") % (s > 0 ? ", " : "")
                                    % (total > 0 ? 100.0 * counts[c][s] / total : 0.0));
                    json += "]";
                }
            } else if (sample_count > 0) {
                for (int c = 0; c < 2; c++)
                    json += str(boost::format(",\n        \"%s_peak\": %.0f,\n"
                                              "        \"%s_rms\": %.1f")
                                % (c == 0 ? "i" : "q") % sample_peak[c]
                                % (c == 0 ? "i" : "q") % sqrt(sample_sumsq[c] / sample_count));
            }

            json += "\n    }";
        }

        json += "\n}\n";

        metafile << json;
        metafile.close();
    };

    while (not stop_signal_called and
           (num_requested_samples > num_total_samps or num_requested_samples == 0)) {

        const auto now = std::chrono::steady_clock::now();

        len = zmq_recv(subscriber, buffer, ZMQ_BUFFER_SIZE, 0);
        if (len < 0)
            continue;

        if (not vrt_process(buffer, sizeof(buffer), &vrt_context, &vrt_packet)) {
            printf("Not a Vita49 packet?\n");
            continue;
        }

        if (vrt_packet.context) {

            if (not start_rx) {

                vrt_print_context(&vrt_context);

                if (vrt_context.sample_rate == 0) {
                    std::cerr << "# Error: the context does not give a sample rate." << std::endl;
                    fatal = true;
                    break;
                }

                /* The sample size is signalled in the payload format of the
                 * context. Streams that do not carry one fall back on --bits. */
                if (vrt_context.has_payload_format) {
                    bits = vrt_context.data_item_size;
                } else {
                    bits = bits_option;
                    fprintf(stderr,
                            "# Warning: context does not signal a sample size, "
                            "assuming %u bit(s) per component\n",
                            bits);
                }

                /* vrt_quantize emits 1 and 2 bit components in the offset binary
                 * coding that VDIF asks for (section 10), so those go in as they
                 * are, and 16 bit components only need their sign bit flipped.
                 * Other sizes would need a coding this tool does not know. */
                if (bits != 1 and bits != 2 and bits != 16) {
                    std::cerr << "# Error: " << bits
                              << " bits per component is not supported, expected 1, 2 or 16."
                              << std::endl;
                    fatal = true;
                    break;
                }

                bits_per_sample = 2 * bits;
                sample_rate     = vrt_context.sample_rate;

                if (vrt_context.has_payload_format and
                    vrt_context.item_packing_field_size != bits_per_sample) {
                    std::cerr << "# Error: the payload packs samples in "
                              << (unsigned)vrt_context.item_packing_field_size
                              << " bits instead of " << bits_per_sample
                              << ", padded payloads are not supported." << std::endl;
                    fatal = true;
                    break;
                }

                if (not vdif_solve_geometry(sample_rate, bits_per_sample, frame_bytes_option,
                                            &payload_bytes, &samples_per_frame, &frames_per_second)) {
                    std::cerr << boost::format("# Error: no data frame size gives an integral number "
                                               "of frames per second at %u samples per second and %u "
                                               "bits per sample.\n")
                                     % sample_rate % bits_per_sample;
                    fatal = true;
                    break;
                }

                if (frame_bytes_given and payload_bytes != frame_bytes_option) {
                    std::cerr << boost::format("# Error: a data array of %u bytes does not give an "
                                               "integral number of frames per second, %u bytes does.\n")
                                     % frame_bytes_option % payload_bytes;
                    fatal = true;
                    break;
                }

                jitter_samples = (uint64_t)llround(jitter * 1e-6 * sample_rate);

                /* A duration is a number of samples to write, so that it counts
                 * the recording and not the time spent waiting for a stream */
                if (total_time > 0)
                    num_requested_samples = (size_t)llround(total_time * sample_rate);

                payload_words = payload_bytes / 4;
                /* One word of slack, so that a copy ending on the last word can
                 * write the part that is shifted out */
                frame.assign(payload_words + 1, 0);

                vdif_init_header(&header, payload_bytes, bits, true, 0, station_id, thread_id);

                std::cout << "# VDIF:\n";
                std::cout << boost::format("#    Data frame: %u bytes (%u byte data array)\n") %
                                 (payload_bytes + VDIF_HEADER_BYTES) % payload_bytes;
                std::cout << boost::format("#    Samples per frame: %u\n") % samples_per_frame;
                std::cout << boost::format("#    Frames per second: %u\n") % frames_per_second;
                std::cout << boost::format("#    Bits per component: %u (complex)\n") % bits;
                if (station_is_ascii)
                    std::cout << boost::format("#    Station ID: %u (%c%c)\n") % station_id %
                                     (char)(station_id >> 8) % (char)(station_id & 0xFF);
                else
                    std::cout << boost::format("#    Station ID: %u\n") % station_id;
                std::cout << boost::format("#    Thread ID: %u\n") % thread_id;

                if (not null) {
                    datafile.open(data_filename.c_str(), std::ofstream::binary);
                    if (not datafile.is_open()) {
                        std::cerr << "# Error: could not open " << data_filename << std::endl;
                        fatal = true;
                        break;
                    }
                }

                start_rf_freq   = vrt_context.rf_freq;
                start_gain      = vrt_context.gain;
                start_bandwidth = vrt_context.bandwidth;

                start_rx = true;

            } else {

                /* The frame geometry follows from the sample rate and the sample
                 * size, so neither can change within a recording */
                if (vrt_context.sample_rate != sample_rate) {
                    std::cerr << boost::format("# Error: the sample rate changed from %u to %u, "
                                               "stopping.\n")
                                     % sample_rate % vrt_context.sample_rate;
                    fatal = true;
                    break;
                }
                if (vrt_context.has_payload_format and vrt_context.data_item_size != bits) {
                    std::cerr << boost::format("# Error: the sample size changed from %u to %u bits, "
                                               "stopping.\n")
                                     % bits % (unsigned)vrt_context.data_item_size;
                    fatal = true;
                    break;
                }

                if (vrt_context.rf_freq != start_rf_freq or vrt_context.gain != start_gain or
                    vrt_context.bandwidth != start_bandwidth) {
                    if (context_changes == 0)
                        fprintf(stderr, "# Warning: the context changed during the recording, the "
                                        "metadata describes its state at the start\n");
                    context_changes++;
                    start_rf_freq   = vrt_context.rf_freq;
                    start_gain      = vrt_context.gain;
                    start_bandwidth = vrt_context.bandwidth;
                }
            }
        }

        if (vrt_packet.extended_context) {
            if (dt_trace)
                dt_process(buffer, sizeof(buffer), &vrt_packet, &dt_ext_context);
            if (tracking)
                tracker_process(buffer, sizeof(buffer), &vrt_packet, &tracker_ext_context);
        }

        if (not (start_rx and vrt_packet.data and vrt_packet.num_rx_samps > 0))
            continue;

        if (vrt_packet.lost_frame) {
            lost_vrt_frames++;
            if (not continue_on_bad_packet)
                break;
        }

        /* Take the payload length from the received message rather than from
         * packet_size, which older versions of vrt_quantize left at the
         * unquantized length after shrinking the payload */
        if (len / 4 <= (int)vrt_packet.offset) {
            fprintf(stderr, "# Error: packet of %d bytes has no payload\n", len);
            continue;
        }

        const uint32_t* payload      = &buffer[vrt_packet.offset];
        const uint64_t payload_words_in = len / 4 - vrt_packet.offset;
        uint64_t       num_samps     = payload_words_in * 32 / bits_per_sample;

        if (num_samps == 0)
            continue;

        /* Occupancy of the quantizer states, or the levels of a 16 bit stream */
        if (bits <= 2) {
            const uint8_t* p      = (const uint8_t*)payload;
            const size_t   nbytes = num_samps * bits_per_sample / 8;
            for (size_t i = 0; i < nbytes; i++)
                byte_hist[p[i]]++;
        } else {
            for (uint64_t i = 0; i < num_samps; i++) {
                const double v[2] = {(double)(int16_t)(payload[i] & 0xFFFF),
                                     (double)(int16_t)(payload[i] >> 16)};
                for (int c = 0; c < 2; c++) {
                    sample_peak[c] = std::max(sample_peak[c], std::fabs(v[c]));
                    sample_sumsq[c] += v[c] * v[c];
                }
            }
            sample_count += num_samps;
        }

        /* Position of the first sample of this packet, in samples since 1970.
         * The fractional timestamp is in picoseconds, the rounding error of the
         * conversion is far below one sample. */
        uint64_t pkt_abs =
            (uint64_t)vrt_packet.integer_seconds_timestamp * sample_rate +
            (uint64_t)llround((double)vrt_packet.fractional_seconds_timestamp * 1e-12 * sample_rate);

        bool just_started = false;

        if (not started) {

            /* The earliest sample the recording may begin at. A stream that
             * starts just after a second boundary would lose almost a whole
             * second waiting for the next one, so by default the recording
             * starts as soon as a data frame can begin, which may be in the
             * middle of a second. Frame numbers are counted within the second
             * either way, so the frames say for themselves where they belong. */
            if (not have_start) {
                start_abs = int_second ? ((uint64_t)vrt_packet.integer_seconds_timestamp + 1) * sample_rate
                                       : pkt_abs;
                if (start_at_timestamp)
                    start_abs = std::max(start_abs, (uint64_t)start_second_requested * sample_rate);
                have_start = true;
            }

            if (pkt_abs + num_samps <= start_abs)
                continue;

            /* A data frame has to begin on its own position in the second */
            uint64_t begin = std::max(start_abs, pkt_abs);
            if (begin % samples_per_frame != 0)
                begin += samples_per_frame - (begin % samples_per_frame);

            if (begin >= pkt_abs + num_samps)
                continue;

            /* VDIF places every sample on the grid that starts at the second,
             * so a timestamp that falls between two samples is rounded to the
             * nearest one. The remainder is the same for the whole recording
             * and is reported in the metadata, for a clock model to take up. */
            const double frac_samples =
                (double)vrt_packet.fractional_seconds_timestamp * 1e-12 * sample_rate;
            start_residual = ((double)llround(frac_samples) - frac_samples) / sample_rate;

            abs_pos      = begin;
            start_abs    = begin;
            started      = true;
            just_started = true;

            if (not vdif_epoch_from_unix((int64_t)(start_abs / sample_rate), &reference_epoch,
                                         &epoch_start)) {
                std::cerr << "# Error: could not determine the VDIF reference epoch." << std::endl;
                fatal = true;
                break;
            }
            vdif_set_epoch(&header, reference_epoch);

            std::cout << boost::format("#    Reference epoch: %u (%s)\n") % (unsigned)reference_epoch %
                             boost::posix_time::to_iso_extended_string(
                                 boost::posix_time::from_time_t(epoch_start));
            std::cout << boost::format("#    Start: %s (second %lld of epoch), frame %u\n") %
                             boost::posix_time::to_iso_extended_string(boost::posix_time::from_time_t(
                                 start_abs / sample_rate)) %
                             (long long int)((int64_t)(start_abs / sample_rate) - epoch_start) %
                             (uint32_t)((start_abs % sample_rate) / samples_per_frame);
            if (start_residual != 0)
                std::cout << boost::format("#    Timing residual: %.1f ns\n") % (start_residual * 1e9);

            write_meta(false);
        }

        /* A source that timestamps in whole microseconds, or one whose clock
         * drifts slowly, puts a packet a sample or two off the position the
         * sample count gives. Data really goes missing in whole packets, so a
         * difference this small is taken as jitter and the count is kept. */
        if (not just_started and jitter_samples > 0) {
            const uint64_t diff = (pkt_abs > abs_pos) ? (pkt_abs - abs_pos) : (abs_pos - pkt_abs);
            if (diff > 0 and diff <= jitter_samples) {
                pkt_abs = abs_pos;
                jitter_packets++;
            }
        }

        uint64_t src_sample = 0;

        if (pkt_abs > abs_pos) {

            const uint64_t missing = pkt_abs - abs_pos;

            if (missing <= (uint64_t)(max_fill * sample_rate)) {
                /* Keep the frame sequence continuous over a gap by filling it
                 * with frames flagged invalid */
                append_samples(NULL, 0, missing);
                fill_samps += missing;
            } else {
                /* Too large to fill, follow the incoming timestamps instead */
                if (frame_fill > 0) {
                    frame_invalid = true;
                    flush_frame();
                }
                uint64_t begin = pkt_abs;
                if (begin % samples_per_frame != 0)
                    begin += samples_per_frame - (begin % samples_per_frame);
                fprintf(stderr,
                        "# Warning: gap of %.3f seconds is larger than --max-fill, "
                        "resynchronising\n",
                        (double)missing / sample_rate);
                abs_pos = begin;
                resyncs++;
            }
        }

        if (abs_pos > pkt_abs) {

            /* A packet that arrives after its position has been written, or one
             * whose first samples were skipped at the start of the recording */
            const uint64_t skip = abs_pos - pkt_abs;

            if (not just_started)
                late_packets++;

            if (skip >= num_samps)
                continue;

            src_sample = skip;
            num_samps -= skip;
        }

        append_samples(payload, src_sample, num_samps);
        num_total_samps += num_samps;

        if (progress and
            std::chrono::duration<double>(now - last_stats).count() > VDIF_STATS_INTERVAL) {

            const double elapsed = std::chrono::duration<double>(now - last_stats).count();

            std::cout << "\t" << boost::format("%.1f fps, %llu invalid") %
                                     ((frames_written - last_stats_frames) / elapsed) %
                                     (unsigned long long)frames_invalid;

            if (bits <= 2) {
                uint64_t counts[2][4];
                expand_hist(byte_hist, counts);
                for (int c = 0; c < 2; c++) {
                    uint64_t total = 0;
                    for (uint32_t s = 0; s < (1u << bits); s++)
                        total += counts[c][s];
                    std::cout << (c == 0 ? ", I: " : "  Q: ");
                    for (uint32_t s = 0; s < (1u << bits); s++)
                        std::cout << boost::format("%s%2.0f") % (s > 0 ? "/" : "") %
                                         (total > 0 ? 100.0 * counts[c][s] / total : 0.0);
                    std::cout << "%";
                }
            } else if (sample_count > 0) {
                for (int c = 0; c < 2; c++)
                    std::cout << boost::format("%s%3.0f dBFS") % (c == 0 ? ", I: " : "  Q: ") %
                                     (20 * log10(std::max(1.0, sample_peak[c]) / 32767.0));
            }

            std::cout << std::endl;

            for (unsigned b = 0; b < 256; b++) {
                byte_hist_total[b] += byte_hist[b];
                byte_hist[b] = 0;
            }

            last_stats        = now;
            last_stats_frames = frames_written;
        }
    }

    if (started) {

        /* Pad the frame that is under way, so that the recording ends on a
         * complete data frame */
        if (frame_fill > 0) {
            frame_invalid = true;
            flush_frame();
        }

        for (unsigned b = 0; b < 256; b++)
            byte_hist_total[b] += byte_hist[b];

        write_meta(true);
    }

    if (not null and datafile.is_open())
        datafile.close();

    if (do_auto_file and started and not null) {

        boost::posix_time::ptime starttime =
            boost::posix_time::from_time_t(start_abs / sample_rate);
        std::string timestring = boost::posix_time::to_iso_extended_string(starttime);
        std::replace(timestring.begin(), timestring.end(), ':', '_');
        std::replace(timestring.begin(), timestring.end(), '-', '_');
        std::replace(timestring.begin(), timestring.end(), 'T', '_');

        boost::format auto_format = boost::format("%s_%s_%.3fMHz_%.2fMsps_ci%u")
                    % (auto_file)
                    % (timestring)
                    % (vrt_context.rf_freq / 1e6)
                    % (vrt_context.sample_rate / 1e6)
                    % bits;

        boost::filesystem::rename(data_filename, auto_format.str() + ".vdif");
        boost::filesystem::rename(meta_filename, auto_format.str() + ".vdif-meta");
    }

    /* Clean up files of a run that never got to write anything */
    if (not null and not started) {
        for (const std::string& name : {data_filename, meta_filename})
            if (boost::filesystem::exists(name) and boost::filesystem::is_empty(name))
                boost::filesystem::remove(name);
    }

    if (started)
        std::cout << boost::format("# %llu frames written, %llu flagged invalid\n") %
                         (unsigned long long)frames_written % (unsigned long long)frames_invalid;

    zmq_close(subscriber);
    zmq_ctx_destroy(context);

    return fatal ? ~0 : 0;
}
