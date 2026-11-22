//
// Copyright 2026 by Thomas Telkamp
//
// SPDX-License-Identifier: MIT
//

#include <boost/format.hpp>
#include <boost/program_options.hpp>

#include <atomic>
#include <chrono>
#include <cmath>
#include <complex>
#include <csignal>
#include <cstring>
#include <iostream>
#include <thread>
#include <vector>

#include <zmq.h>

// VRT
#include <vrt/vrt_init.h>
#include <vrt/vrt_string.h>
#include <vrt/vrt_types.h>
#include <vrt/vrt_util.h>
#include <vrt/vrt_write.h>
#include <vrt/vrt_read.h>

#include <inttypes.h>

#include <libairspyhf/airspyhf.h>

// VRT tools functions
#include "vrt-tools.h"
#include "vrt-rx-queue.h"

#define HF_ATT_MAX (8)
#define HF_MAX_FREQ (31e6)      // the HF input, and its gain controls, end here

// Attenuation of the HF attenuator steps and gain of the HF preamplifier, in
// dB. Measured with a signal generator at 10 MHz on one Airspy HF+ (firmware
// R3.0.6-CD), so the gain reported in the context is approximate. The driver
// documents steps of 6 dB and a preamplifier that is compensated digitally;
// the measurement matched neither.
static const double hf_att_db[HF_ATT_MAX + 1] = { 0.0, 3.2, 7.8, 13.3, 17.6, 23.7, 30.6, 36.7, 41.3 };
#define HF_LNA_DB (14.0)

unsigned long long num_total_samps = 0;

namespace po = boost::program_options;

// The receive queue, and what the callback needs to fill it.
struct hf_rx_state {
    vrt_rx_queue queue;

    // Float full scale to int16, including the --scale gain.
    float scale = 32768.0f;

    // Used by the callback thread only.
    std::vector<int16_t> converted;

    std::atomic<uint64_t> clipped{0};
};

static bool stop_signal_called = false;
void sig_int_handler(int)
{
    stop_signal_called = true;
}

template <typename samp_type> inline float get_abs_val(samp_type t)
{
    return std::fabs(t);
}

int parse_u64(const char* s, uint64_t* const value) {
    uint_fast8_t base = 10;
    char* s_end;
    uint64_t u64_value;

    if( strlen(s) > 2 ) {
        if( s[0] == '0' ) {
            if( (s[1] == 'x') || (s[1] == 'X') ) {
                base = 16;
                s += 2;
            } else if( (s[1] == 'b') || (s[1] == 'B') ) {
                base = 2;
                s += 2;
            }
        }
    }

    s_end = (char*)s;
    u64_value = strtoull(s, &s_end, base);
    if( (s != s_end) && (*s_end == 0) ) {
        *value = u64_value;
        return AIRSPYHF_SUCCESS;
    } else {
        return AIRSPYHF_ERROR;
    }
}

// libairspyhf delivers float samples with a full scale of 1.0.
//
// The device itself sends 16 bit samples, but the library scales them down by
// the gain of the firmware's decimation filter, which depends on the sample
// rate (measured with firmware R3.0.6-CD: 9 dB at 384 kS/s, 13 dB at 256 kS/s,
// 16 dB at 192 kS/s). Without a signal the receiver noise is then less than one
// int16 step. The library does not report that gain; --scale makes up for it.
static inline int16_t float_to_int16(float value, float scale, uint64_t* clipped)
{
    const float scaled = value * scale;

    if (scaled >= 32767.0f) {
        (*clipped)++;
        return 32767;
    }
    if (scaled <= -32768.0f) {
        (*clipped)++;
        return -32768;
    }
    return (int16_t)lrintf(scaled);
}

int rx_callback(airspyhf_transfer_t* transfer)
{
    const struct vrt_time_ps arrival = vrt_time_now();

    hf_rx_state* rx = (hf_rx_state*)transfer->ctx;

    if (not rx->queue.enabled.load() or transfer->samples == nullptr or transfer->sample_count <= 0)
        return 0;

    const size_t count = (size_t)transfer->sample_count;
    uint64_t clipped = 0;

    rx->converted.resize(2*count);
    for (size_t i = 0; i < count; i++) {
        rx->converted[2*i]   = float_to_int16(transfer->samples[i].re, rx->scale, &clipped);
        rx->converted[2*i+1] = float_to_int16(transfer->samples[i].im, rx->scale, &clipped);
    }
    if (clipped > 0)
        rx->clipped += clipped;

    vrt_rx_push(&rx->queue, arrival, rx->converted.data(), count, transfer->dropped_samples);

    return 0;
}

int main(int argc, char* argv[])
{
    // variables to be set by po
    std::string merge_address, dev_given, agc_threshold;
    uint16_t instance, port, merge_port;
    int hwm;
    int32_t ppb = 0;
    uint32_t att = 0;
    double rate, freq, setup_time, if_freq, scale_db;
    bool merge;

    // Airspy HF+
    airspyhf_device_t* device = 0;
    char airspyhf_version[MAX_VERSION_STRING_SIZE + 1] = "";
    uint64_t serial_number_val = 0;
    uint32_t sample_rate_u32 = 0;
    uint32_t nsrates = 0;

    // setup the program options
    po::options_description desc("Allowed options");
    // clang-format off
    desc.add_options()
        ("help", "help message")
        ("rate", po::value<double>(&rate)->default_value(768000), "rate of incoming samples")
        ("freq", po::value<double>(&freq)->default_value(10e6), "RF center frequency in Hz")
        ("serial", po::value<std::string>(&dev_given), "device serial")
        ("if-freq", po::value<double>(&if_freq)->default_value(0.0), "IF center frequency in Hz")
        ("setup", po::value<double>(&setup_time)->default_value(1.0), "seconds of setup time")
        ("att", po::value<uint32_t>(&att)->default_value(0), "HF attenuator step, 0 to 8, about 41 dB in all (AGC off)")
        ("lna", "enable the HF preamplifier")
        ("agc", "enable the HF AGC (default is off, for a fixed gain)")
        ("agc-threshold", po::value<std::string>(&agc_threshold)->default_value("low"), "HF AGC threshold, low or high (AGC on)")
        ("scale", po::value<double>(&scale_db)->default_value(0.0), "digital gain in dB applied before the conversion to 16 bits")
        ("ppb", po::value<int32_t>(&ppb), "frequency calibration in ppb, for this run (default is the value stored in the device)")
        ("progress", "periodically display short-term bandwidth")
        ("stats", "show average bandwidth on exit")
        ("int-second", "align start of reception to integer second")
        ("instance", po::value<uint16_t>(&instance)->default_value(0), "VRT ZMQ instance")
        ("port", po::value<uint16_t>(&port), "VRT ZMQ port")
        ("merge", po::value<bool>(&merge)->default_value(true), "Merge another VRT ZMQ stream (SUB connect)")
        ("merge-port", po::value<uint16_t>(&merge_port)->default_value(50011), "VRT ZMQ merge port")
        ("merge-address", po::value<std::string>(&merge_address)->default_value("localhost"), "VRT ZMQ merg address")
        ("hwm", po::value<int>(&hwm)->default_value(10000), "VRT ZMQ HWM")
    ;
    // clang-format on
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // print the help message
    if (vm.count("help")) {
        std::cout << boost::format("Airspy HF+ RX samples to VRT. %s") % desc << std::endl;
        std::cout << std::endl
                  << "This application streams data from an Airspy HF+ "
                     "device to VRT.\n"
                  << std::endl;
        return ~0;
    }

    bool bw_summary    = vm.count("progress") > 0;
    bool stats         = vm.count("stats") > 0;
    bool int_second    = vm.count("int-second") > 0;
    bool serial_number = vm.count("serial") > 0;
    bool hf_lna        = vm.count("lna") > 0;
    bool hf_agc        = vm.count("agc") > 0;
    bool set_ppb       = vm.count("ppb") > 0;

    if (serial_number) {
        if (parse_u64(dev_given.c_str(), &serial_number_val) != AIRSPYHF_SUCCESS) {
            fprintf(stderr, "argument error: invalid serial number: %s\n", dev_given.c_str());
            return EXIT_FAILURE;
        }
    }

    if (att > HF_ATT_MAX) {
        fprintf(stderr, "argument error: att out of range (0 to %d)\n", HF_ATT_MAX);
        return EXIT_FAILURE;
    }

    if (agc_threshold != "low" and agc_threshold != "high") {
        fprintf(stderr, "argument error: agc-threshold is low or high\n");
        return EXIT_FAILURE;
    }

    // The tuner takes whole Hz.
    freq = std::round(freq - if_freq);

    if (freq < 0 or freq > 4294967295.0) {
        fprintf(stderr, "argument error: frequency %.0f Hz out of range\n", freq);
        return EXIT_FAILURE;
    }
    if (not ((freq >= 9e3 and freq <= 31e6) or (freq >= 60e6 and freq <= 260e6)))
        fprintf(stderr, "Warning: frequency %.0f Hz is outside the specified ranges of 9 kHz to 31 MHz and 60 to 260 MHz.\n", freq);

    /* VRT init */
    struct vrt_packet p;
    vrt_init_packet(&p);
    vrt_init_data_packet(&p);

    p.fields.stream_id = 1;

    uint16_t main_port;

    if (vm.count("port") > 0) {
        main_port = port;
    } else {
        main_port = DEFAULT_MAIN_PORT + MAX_CHANNELS*instance;
    }

    // ZMQ
    void *zmq_server;

    // Checked explicitly: assert() is compiled out in release builds (NDEBUG),
    // so a failed bind (e.g. port already in use) would go unnoticed.
    auto zmq_fail = [](const char* what, const std::string& endpoint) {
        fprintf(stderr, "Error: unable to %s %s: %s\n", what, endpoint.c_str(), zmq_strerror(errno));
        return EXIT_FAILURE;
    };

    void *context = zmq_ctx_new();
    void *responder = zmq_socket(context, ZMQ_PUB);
    if (zmq_setsockopt (responder, ZMQ_SNDHWM, &hwm, sizeof hwm) != 0)
        return zmq_fail("set the high water mark on", "the VRT socket");

    std::string connect_string = "tcp://*:" + std::to_string(main_port);
    if (zmq_bind(responder, connect_string.c_str()) != 0)
        return zmq_fail("bind", connect_string);
    zmq_server = responder;

    // Merge
    void *merge_zmq = zmq_socket(context, ZMQ_SUB);
    if (merge) {
        connect_string = "tcp://" + merge_address + ":" + std::to_string(merge_port);
        if (zmq_connect(merge_zmq, connect_string.c_str()) != 0)
            return zmq_fail("connect to", connect_string);
        zmq_setsockopt(merge_zmq, ZMQ_SUBSCRIBE, "", 0);
    }

    // open the device

    if (serial_number) {
        if (airspyhf_open_sn(&device, serial_number_val) != AIRSPYHF_SUCCESS) {
            fprintf(stderr, "airspyhf_open_sn() failed\n");
            return EXIT_FAILURE;
        }
    } else {
        if (airspyhf_open(&device) != AIRSPYHF_SUCCESS) {
            fprintf(stderr, "airspyhf_open() failed\n");
            return EXIT_FAILURE;
        }
    }

    // check sample rates
    airspyhf_get_samplerates(device, &nsrates, 0);
    std::vector<uint32_t> supported_samplerates(nsrates);
    if (nsrates > 0)
        airspyhf_get_samplerates(device, supported_samplerates.data(), nsrates);

    fprintf(stderr, "Supported sample rates:");
    for (uint32_t s = 0; s < nsrates; s++) {
        fprintf(stderr, " %d kS/s", supported_samplerates[s]/1000);
        if (rate == supported_samplerates[s])
            sample_rate_u32 = supported_samplerates[s];
    }
    fprintf(stderr, "\n");

    if (sample_rate_u32 == 0) {
        fprintf(stderr, "unsupported sample rate: %f\n", rate);
        airspyhf_close(device);
        return EXIT_FAILURE;
    }

    // set sample rate
    if (airspyhf_set_samplerate(device, sample_rate_u32) != AIRSPYHF_SUCCESS) {
        fprintf(stderr, "airspyhf_set_samplerate() failed: %d\n", sample_rate_u32);
        airspyhf_close(device);
        return EXIT_FAILURE;
    }

    rate = sample_rate_u32;

    // At a zero-IF rate the library tunes the LO 5 kHz above the requested
    // frequency and shifts the samples back, so an LO spur sits 5 kHz off centre.
    fprintf(stderr, "Sample rate: %f MS/s IQ (%s)\n", sample_rate_u32 * 0.000001f,
        airspyhf_is_low_if(device) ? "low IF" : "zero IF");

    airspyhf_read_partid_serialno_t read_partid_serialno;

    if (airspyhf_board_partid_serialno_read(device, &read_partid_serialno) != AIRSPYHF_SUCCESS) {
        fprintf(stderr, "airspyhf_board_partid_serialno_read() failed\n");
        airspyhf_close(device);
        return EXIT_FAILURE;
    } else
        fprintf(stderr, "Device Serial Number: 0x%08X%08X\n",
                read_partid_serialno.serial_no[0],
                read_partid_serialno.serial_no[1]);

    if (airspyhf_version_string_read(device, &airspyhf_version[0], MAX_VERSION_STRING_SIZE) != AIRSPYHF_SUCCESS)
        fprintf(stderr, "airspyhf_version_string_read() failed\n");

    fprintf(stderr, "Firmware Version: %s\n", airspyhf_version);

    // Frequency calibration. The library applies it when it tunes; it is only
    // stored in the device with airspyhf_calibrate.
    if (set_ppb) {
        if (airspyhf_set_calibration(device, ppb) != AIRSPYHF_SUCCESS)
            fprintf(stderr, "airspyhf_set_calibration() failed\n");
    }
    int32_t calibration_ppb = 0;
    if (airspyhf_get_calibration(device, &calibration_ppb) == AIRSPYHF_SUCCESS)
        fprintf(stderr, "Calibration: %d ppb%s\n", calibration_ppb, set_ppb ? "" : " (stored in the device)");

    // Gain. These controls act on the HF input, and older firmware does not
    // have them, so a failure is reported but is not fatal.
    if (airspyhf_set_hf_agc(device, hf_agc ? 1 : 0) != AIRSPYHF_SUCCESS)
        fprintf(stderr, "airspyhf_set_hf_agc() failed\n");

    if (hf_agc) {
        if (airspyhf_set_hf_agc_threshold(device, agc_threshold == "high" ? 1 : 0) != AIRSPYHF_SUCCESS)
            fprintf(stderr, "airspyhf_set_hf_agc_threshold() failed\n");
    } else {
        if (airspyhf_set_hf_att(device, (uint8_t)att) != AIRSPYHF_SUCCESS)
            fprintf(stderr, "airspyhf_set_hf_att() failed\n");
    }

    if (airspyhf_set_hf_lna(device, hf_lna ? 1 : 0) != AIRSPYHF_SUCCESS)
        fprintf(stderr, "airspyhf_set_hf_lna() failed\n");

    // Gain relative to the HF input with no attenuation, including the digital
    // scale. The attenuator, preamplifier and AGC only act on the HF input, and
    // with the AGC on the analog gain is not known.
    const bool hf_input = freq <= HF_MAX_FREQ;
    double gain = scale_db;

    if (hf_input) {
        if (not hf_agc)
            gain -= hf_att_db[att];
        if (hf_lna)
            gain += HF_LNA_DB;

        if (hf_agc)
            fprintf(stderr, "HF AGC: on (threshold %s), HF preamplifier: %s\n", agc_threshold.c_str(), hf_lna ? "on" : "off");
        else
            fprintf(stderr, "HF AGC: off, HF attenuator: step %u (about %.1f dB), HF preamplifier: %s\n", att, hf_att_db[att], hf_lna ? "on" : "off");
    } else if (hf_agc or hf_lna or att > 0) {
        fprintf(stderr, "Warning: the AGC, attenuator and preamplifier only act on the HF input, below %.0f MHz.\n", HF_MAX_FREQ/1e6);
    }

    // Samples received during setup are discarded: the callback only queues
    // blocks once the queue is enabled, just before the main loop.
    fprintf(stderr, "Digital scale: %+.1f dB\n", scale_db);

    hf_rx_state rx;
    rx.scale            = 32768.0f * (float)std::pow(10.0, scale_db/20.0);
    rx.queue.rate       = rate;
    rx.queue.max_blocks = 256; // transfers are only a few thousand samples

    if (airspyhf_start(device, rx_callback, &rx) != AIRSPYHF_SUCCESS) {
        fprintf(stderr, "airspyhf_start() failed.\n");
        airspyhf_close(device);
        return EXIT_FAILURE;
    }

    // The device tunes in steps of 1 kHz; the library takes out the rest, and
    // the calibration, by rotating the samples.
    if (airspyhf_set_freq(device, (uint32_t)freq) != AIRSPYHF_SUCCESS) {
        fprintf(stderr, "airspyhf_set_freq() failed.\n");
        airspyhf_stop(device);
        airspyhf_close(device);
        return EXIT_FAILURE;
    }

    // Sleep setup time
    std::this_thread::sleep_for(std::chrono::milliseconds(int64_t(1000 * setup_time)));

    struct vrt_time_ps time_now = vrt_time_now();

    // Receive

    std::signal(SIGINT, &sig_int_handler);
    std::cout << "Press Ctrl + C to stop streaming..." << std::endl;

    size_t samps_per_buff = VRT_SAMPLES_PER_PACKET;

    uint32_t buffer[SIZE];

    // Merged messages come from other tools, so they get a buffer of their
    // own that is large enough for any VRT packet.
    std::vector<uint8_t> merge_buffer(ZMQ_BUFFER_SIZE);

    bool first_frame = true;
    bool context_changed = true;

    /* Warn if not standards compliant */
    if (vrt_is_platform_little_endian()) {
        printf("Warning: little endian support is work in progress.\n");
    }

    // time keeping
    auto start_time = std::chrono::steady_clock::now();

    // Track time and samps between updating the BW summary
    auto last_update                     = start_time;
    auto last_context                    = start_time;
    unsigned long long last_update_samps = 0;

    uint32_t frame_count = 0;

    vrt_rx_packetizer packetizer(samps_per_buff);
    std::vector<int16_t>& bodydata = packetizer.body;

    struct vrt_time_ps time_first_sample = {0, 0};

    uint64_t dropped_reported = 0;
    auto last_drop_report = start_time;
    uint64_t clipped_reported = 0;
    auto last_clip_report = start_time;
    bool stream_failed = false;

    // flush merge queue
    if (merge)
        while ( zmq_recv(merge_zmq, merge_buffer.data(), merge_buffer.size(), ZMQ_NOBLOCK) > 0 ) { }

    // Start accepting samples. With --int-second the callback starts at the
    // sample that falls on the next integer second.
    rx.queue.int_second = int_second;
    rx.queue.enabled = true;

    while (not stop_signal_called) {

        const auto now = std::chrono::steady_clock::now();

        // libairspyhf stops for good on the first failed USB transfer.
        if (not airspyhf_is_streaming(device)) {
            fprintf(stderr, "Error: Airspy HF+ stopped streaming (USB transfer failed).\n");
            stream_failed = true;
            break;
        }

        vrt_rx_block block;

        if (vrt_rx_pop(&rx.queue, &block, &time_first_sample, std::chrono::milliseconds(1))) {

            packetizer.feed(block, [&](uint64_t packet_first_sample, bool after_gap) {

                // Samples were lost ahead of this packet: skip a packet count
                // so that consumers see a lost frame.
                if (after_gap)
                    frame_count++;

                // Timestamp of the first sample in the packet, counted in
                // samples from the start of the stream.
                time_now = vrt_time_add_samples(time_first_sample, packet_first_sample, rate);

                if (first_frame) {
                    std::cout << boost::format(
                                     "First frame: %u full secs, %.09f frac secs")
                                     % time_now.seconds
                                     % ((double)time_now.frac_ps/1e12)
                              << std::endl;
                    first_frame = false;
                }

                num_total_samps += samps_per_buff;

                p.body = bodydata.data();
                p.header.packet_count = (uint8_t)frame_count%16;
                p.fields.integer_seconds_timestamp = time_now.seconds;
                p.fields.fractional_seconds_timestamp = time_now.frac_ps;

                zmq_msg_t msg;
                zmq_msg_init_size (&msg, SIZE*4);

                int32_t rv = vrt_write_packet(&p, zmq_msg_data(&msg), SIZE, true);
                if (rv < 0) {
                    fprintf(stderr, "Failed to write packet: %s\n", vrt_string_error(rv));
                }

                frame_count++;

                // VRT
                zmq_msg_send(&msg, zmq_server, 0);
                zmq_msg_close(&msg);

                const auto time_since_last_context = now - last_context;
                if (time_since_last_context > std::chrono::milliseconds(VRT_CONTEXT_INTERVAL)) {

                    last_context = now;

                    // VITA 49.2
                    /* Initialize to reasonable values */
                    struct vrt_packet pc;
                    vrt_init_packet(&pc);

                    /* VRT Configure. Note that context packets cannot have a trailer word. */
                    vrt_init_context_packet(&pc);

                    pc.fields.integer_seconds_timestamp = time_now.seconds;
                    pc.fields.fractional_seconds_timestamp = time_now.frac_ps;

                    pc.fields.stream_id = p.fields.stream_id;

                    pc.if_context.bandwidth                         = rate;
                    pc.if_context.sample_rate                       = rate;
                    pc.if_context.rf_reference_frequency            = freq+if_freq;
                    pc.if_context.rf_reference_frequency_offset     = 0;
                    pc.if_context.if_reference_frequency            = if_freq; // 0 for Zero-IF
                    pc.if_context.gain.stage1                       = gain;
                    pc.if_context.gain.stage2                       = 0;

                    // The Airspy HF+ has no external reference input.
                    pc.if_context.state_and_event_indicators.reference_lock = false;

                    pc.if_context.state_and_event_indicators.calibrated_time = false;

                    pc.if_context.context_field_change_indicator = context_changed;

                    int32_t rv = vrt_write_packet(&pc, buffer, SIZE, true);
                    if (rv < 0) {
                        fprintf(stderr, "Failed to write packet: %s\n", vrt_string_error(rv));
                    } else {
                        // ZMQ
                        zmq_send (zmq_server, buffer, rv*4, 0);
                    }

                    context_changed = false;
                }

                if (bw_summary) {
                    last_update_samps += samps_per_buff;
                    const auto time_since_last_update = now - last_update;
                    if (time_since_last_update > std::chrono::seconds(1)) {

                        const double time_since_last_update_s =
                            std::chrono::duration<double>(time_since_last_update).count();
                        const double rate = double(last_update_samps) / time_since_last_update_s;
                        last_update_samps = 0;
                        last_update       = now;

                        vrt_rx_print_progress(rate, 0, bodydata.data(), samps_per_buff);

                    }
                }
            });

            vrt_rx_recycle(&rx.queue, &block);
        }

        // Report lost samples, at most once a second.
        const uint64_t dropped_total = rx.queue.dropped_driver + rx.queue.dropped_queue + packetizer.dropped_partial;
        if (dropped_total != dropped_reported and now - last_drop_report > std::chrono::seconds(1)) {
            fprintf(stderr, "Warning: lost %llu samples (%llu in total).\n",
                (unsigned long long)(dropped_total - dropped_reported),
                (unsigned long long)dropped_total);
            dropped_reported = dropped_total;
            last_drop_report = now;
        }

        // Report clipping in the conversion to 16 bits, at most once a second.
        const uint64_t clipped_total = rx.clipped;
        if (clipped_total != clipped_reported and now - last_clip_report > std::chrono::seconds(1)) {
            fprintf(stderr, "Warning: %llu sample values clipped in the conversion to 16 bits (%llu in total), reduce --scale.\n",
                (unsigned long long)(clipped_total - clipped_reported),
                (unsigned long long)clipped_total);
            clipped_reported = clipped_total;
            last_clip_report = now;
        }

        // Merge
        if (merge) {
            int mergelen;
            while ( (mergelen = zmq_recv(merge_zmq, merge_buffer.data(), merge_buffer.size(), ZMQ_NOBLOCK)) > 0  ) {
                // zmq_recv() truncates a longer message but returns its full length.
                if ((size_t)mergelen > merge_buffer.size()) {
                    fprintf(stderr, "Warning: merged message of %d bytes is too long, dropped.\n", mergelen);
                    continue;
                }
                zmq_msg_t msg;
                zmq_msg_init_size (&msg, mergelen);
                memcpy (zmq_msg_data(&msg), merge_buffer.data(), mergelen);
                zmq_msg_send(&msg, zmq_server, 0);
                zmq_msg_close(&msg);
            }
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

        const uint64_t dropped_driver = rx.queue.dropped_driver;
        const uint64_t dropped_queue  = rx.queue.dropped_queue;
        std::cout << boost::format("Lost %d samples (%d in the driver, %d in the queue, %d in cut packets).")
                         % (dropped_driver + dropped_queue + packetizer.dropped_partial)
                         % dropped_driver % dropped_queue % packetizer.dropped_partial
                  << std::endl;

        const uint64_t clipped = rx.clipped;
        std::cout << boost::format("Clipped %d sample values in the conversion to 16 bits.") % clipped
                  << std::endl;
    }

    /* clean up */
    if (airspyhf_stop(device) != AIRSPYHF_SUCCESS) {
        fprintf(stderr, "airspyhf_stop() failed\n");
    }

    if (airspyhf_close(device) != AIRSPYHF_SUCCESS) {
        fprintf(stderr, "airspyhf_close() failed\n");
    }

    // finished
    std::cout << std::endl << "Done!" << std::endl << std::endl;

    return stream_failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
