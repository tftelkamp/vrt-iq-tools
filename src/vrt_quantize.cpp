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

#include <algorithm>
#include <chrono>
#include <cmath>
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

/* Optimal threshold of the 4-level (2-bit) quantizer, in units of sigma, for
 * Gaussian noise. */
#define DEFAULT_THRESHOLD 0.9816

/* sqrt(pi/2): converts the mean absolute deviation of Gaussian noise to sigma. */
#define MAD_TO_SIGMA 1.2533141373155003

/* Interval between quantizer statistics lines, in seconds. */
#define QUANTIZE_STATS_INTERVAL 1.0

/* Scratch space for rewriting a context packet. */
#define VRT_CONTEXT_BUFFER_WORDS 1024

/* Reconstruction levels used when unpacking back to 16 bit. */
static const int16_t levels_1_bit[2] = {-1, 1};
static const int16_t levels_2_bit[4] = {-3, -1, 1, 3};

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

/* Accumulate the per-component mean and dispersion of one packet of 16-bit
 * complex samples. The dispersion is the mean square in RMS mode, or the mean
 * absolute deviation around ref in MAD mode. Index 0 is the real component,
 * index 1 the imaginary one; they are treated as independent samplers. */
static void accumulate_levels(const uint32_t* payload,
                              uint32_t num_samps,
                              bool mad,
                              const double* ref,
                              double* mean,
                              double* dispersion)
{
    double sum[2] = {0, 0};
    double acc[2] = {0, 0};

    for (uint32_t i = 0; i < num_samps; i++) {
        const uint32_t word = payload[i];
        const double   v[2] = {(double)(int16_t)(word & 0xFFFF),
                               (double)(int16_t)((word >> 16) & 0xFFFF)};
        for (int c = 0; c < 2; c++) {
            sum[c] += v[c];
            acc[c] += mad ? std::fabs(v[c] - ref[c]) : v[c] * v[c];
        }
    }

    for (int c = 0; c < 2; c++) {
        mean[c]       = sum[c] / num_samps;
        dispersion[c] = acc[c] / num_samps;
    }
}

/* Update packet_size in the header of a data packet whose payload length has
 * changed. Returns the new packet length in bytes, or -1 on failure. */
static int update_packet_size(uint32_t* buffer, uint32_t offset, size_t payload_words)
{
    struct vrt_header h;

    int32_t rv = vrt_read_header(buffer, ZMQ_BUFFER_SIZE, &h, true);
    if (rv < 0) {
        fprintf(stderr, "Failed to parse header: %s\n", vrt_string_error(rv));
        return -1;
    }

    h.packet_size = (uint16_t)(offset + payload_words);

    rv = vrt_write_header(&h, buffer, ZMQ_BUFFER_SIZE, true);
    if (rv < 0) {
        fprintf(stderr, "Failed to write header: %s\n", vrt_string_error(rv));
        return -1;
    }

    return (int)(offset + payload_words) * 4;
}

/* Copy a context packet, setting the payload format to describe the sample
 * size the data packets now carry, so that downstream tools (and --unpack)
 * can tell how many bits per component there are. Returns the length of the
 * rewritten packet in bytes, or -1 if it could not be rewritten. */
static int rewrite_payload_format(const uint32_t* in,
                                  int len,
                                  uint32_t bits,
                                  uint32_t* out,
                                  int32_t out_words)
{
    struct vrt_header     h;
    struct vrt_fields     f;
    struct vrt_if_context c;

    const int32_t words = len / 4;

    int32_t rv = vrt_read_header(in, words, &h, true);
    if (rv < 0) {
        fprintf(stderr, "Failed to parse header: %s\n", vrt_string_error(rv));
        return -1;
    }
    int32_t offset = rv;

    rv = vrt_read_fields(&h, in + offset, words - offset, &f, true);
    if (rv < 0) {
        fprintf(stderr, "Failed to parse fields section: %s\n", vrt_string_error(rv));
        return -1;
    }
    offset += rv;

    rv = vrt_read_if_context(in + offset, words - offset, &c, true);
    if (rv < 0) {
        fprintf(stderr, "Failed to parse IF context section: %s\n", vrt_string_error(rv));
        return -1;
    }

    /* One item packing field holds a complex sample, one data item a single
     * component. Unlike context_type, these are the raw Vita49 fields, which
     * hold one less than the actual size. */
    c.has.data_packet_payload_format                     = true;
    c.data_packet_payload_format.packing_method          = VRT_PM_LINK_EFFICIENT;
    c.data_packet_payload_format.real_or_complex         = VRT_ROC_COMPLEX_CARTESIAN;
    c.data_packet_payload_format.data_item_format        = VRT_DIF_SIGNED_FIXED_POINT;
    c.data_packet_payload_format.sample_component_repeat = false;
    c.data_packet_payload_format.item_packing_field_size = (uint8_t)(2 * bits - 1);
    c.data_packet_payload_format.data_item_size          = (uint8_t)(bits - 1);

    struct vrt_packet pc;
    vrt_init_packet(&pc);
    pc.header     = h;
    pc.fields     = f;
    pc.if_context = c;

    rv = vrt_write_packet(&pc, out, out_words, true);
    if (rv < 0) {
        fprintf(stderr, "Failed to write context packet: %s\n", vrt_string_error(rv));
        return -1;
    }

    return rv * 4;
}

int main(int argc, char* argv[])
{

    // variables to be set by po
    std::string file, type, zmq_address, level_mode;
    uint16_t pub_instance, instance, main_port, port, pub_port;
    uint32_t channel, bits;
    int hwm;
    double threshold, tau, sigma_re, sigma_im;

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
        ("bits", po::value<uint32_t>(&bits)->default_value(1), "bits per component (1 or 2)")
        ("level-mode", po::value<std::string>(&level_mode)->default_value("rms"), "level estimator: rms or mad")
        ("threshold", po::value<double>(&threshold)->default_value(DEFAULT_THRESHOLD), "2-bit threshold, in units of sigma")
        ("tau", po::value<double>(&tau)->default_value(1.0), "time constant of the level estimator [seconds]")
        ("sigma-re", po::value<double>(&sigma_re)->default_value(0), "fixed sigma for the real component")
        ("sigma-im", po::value<double>(&sigma_im)->default_value(0), "fixed sigma for the imaginary component")
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
                  << "This application 1-bit or 2-bit quantizes a VRT stream, and unpacks it\n"
                     "back to 16 bit with --unpack. The sample size is signalled downstream in\n"
                     "the data_item_size field of the context packet payload format.\n"
                     "\n"
                     "2-bit quantization uses the 4-level VLBI encoding, offset binary with\n"
                     "00 = -high, 01 = -low, 10 = +low, 11 = +high, thresholded at 0.9816 sigma.\n"
                     "Sigma and the DC offset are estimated per component over a sliding window,\n"
                     "either from the RMS or from the mean absolute deviation, the latter being\n"
                     "less sensitive to impulsive RFI. 1-bit quantization is a pure sign bit and\n"
                     "needs no level estimate.\n"
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

    if (bits != 1 and bits != 2) {
        std::cerr << "Only 1 or 2 bits per component are supported." << std::endl;
        return ~0;
    }

    if (level_mode != "rms" and level_mode != "mad") {
        std::cerr << "Unknown level mode \"" << level_mode << "\", expected rms or mad." << std::endl;
        return ~0;
    }

    bool mad = (level_mode == "mad");

    // Fixed sigma per component, zero when it is to be estimated
    double sigma_fixed[2] = {sigma_re, sigma_im};

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
    uint32_t context_buffer[VRT_CONTEXT_BUFFER_WORDS];

    uint64_t num_total_samps = 0;

    // Track time and samps between updating the BW summary
    auto last_update = start_time;
    uint64_t last_update_samps = 0;

    // Level estimator state, index 0 is the real, index 1 the imaginary component
    bool   levels_valid = false;
    double level_mean[2]       = {0, 0};
    double level_dispersion[2] = {0, 0};
    double sigma[2]            = {0, 0};
    double level_threshold[2]  = {0, 0};

    // Quantizer state occupancy since the last statistics line
    uint64_t state_count[2][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}};
    uint64_t state_samps       = 0;
    auto     last_stats        = start_time;

    // Bits per component of the incoming stream, when unpacking
    uint32_t unpack_bits    = bits;
    bool     unpack_warned  = false;

    // Reset the level estimate when the receiver gain changes
    int32_t last_gain = 0;

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
            last_gain = vrt_context.gain;
            // Possibly do something with context here
            // vrt_context
        }

        if (vrt_packet.context) {

            if (unpack) {
                // The sample size of the incoming stream. Streams quantized by
                // an older version do not signal it, fall back on --bits then.
                if (vrt_context.has_payload_format and vrt_context.data_item_size < 16) {
                    unpack_bits = vrt_context.data_item_size;
                } else {
                    unpack_bits = bits;
                    if (not unpack_warned) {
                        fprintf(stderr,
                                "# Warning: context does not signal a quantized sample size, "
                                "assuming %u bit(s) per component\n",
                                unpack_bits);
                        unpack_warned = true;
                    }
                }
            } else if (vrt_context.gain != last_gain) {
                // A gain change invalidates the level estimate
                levels_valid = false;
                last_gain    = vrt_context.gain;
            }

            // Signal the sample size the data packets now carry
            const uint32_t bits_out = unpack ? 16 : bits;
            const int      ctx_len =
                rewrite_payload_format(buffer, len, bits_out, context_buffer, VRT_CONTEXT_BUFFER_WORDS);

            zmq_msg_t msg;
            if (ctx_len > 0) {
                zmq_msg_init_size (&msg, ctx_len);
                memcpy (zmq_msg_data(&msg), context_buffer, ctx_len);
            } else {
                // Could not be rewritten, forward it unchanged
                zmq_msg_init_size (&msg, len);
                memcpy (zmq_msg_data(&msg), buffer, len);
            }
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

        if (progress && vrt_packet.data && !unpack)
            show_progress_stats(
                now,
                &last_update,
                &last_update_samps,
                &buffer[vrt_packet.offset],
                vrt_packet.num_rx_samps, 0
            );

        if (start_rx and vrt_packet.data and vrt_packet.num_rx_samps > 0) {

            size_t new_data_len_words;
            size_t num_samps;

            if (!unpack) {

                // For a 16-bit complex stream one word holds one sample
                const uint32_t num_samps      = vrt_packet.num_rx_samps;
                const uint32_t bits_per_samp  = 2 * bits;

                // Estimate the level of this packet. Not needed for a plain
                // sign bit, but still worth having for the statistics line.
                if (bits > 1 or progress) {

                    double packet_mean[2], packet_dispersion[2];

                    accumulate_levels(&buffer[vrt_packet.offset], num_samps, mad, level_mean,
                                      packet_mean, packet_dispersion);

                    if (not levels_valid) {
                        // The first packet has no earlier mean to deviate from,
                        // so redo the pass now that its own mean is known
                        if (mad) {
                            double redo_mean[2];
                            accumulate_levels(&buffer[vrt_packet.offset], num_samps, mad, packet_mean,
                                              redo_mean, packet_dispersion);
                        }
                        for (int c = 0; c < 2; c++) {
                            level_mean[c]       = packet_mean[c];
                            level_dispersion[c] = packet_dispersion[c];
                        }
                        levels_valid = true;
                    } else {
                        // Track the level slowly, so that the estimate follows
                        // gain and system temperature drift but not the signal
                        double alpha = 1.0;
                        if (vrt_context.sample_rate > 0 and tau > 0)
                            alpha = 1.0 - exp(-((double)num_samps / vrt_context.sample_rate) / tau);
                        for (int c = 0; c < 2; c++) {
                            level_mean[c]       += alpha * (packet_mean[c] - level_mean[c]);
                            level_dispersion[c] += alpha * (packet_dispersion[c] - level_dispersion[c]);
                        }
                    }

                    for (int c = 0; c < 2; c++) {
                        if (mad)
                            sigma[c] = MAD_TO_SIGMA * level_dispersion[c];
                        else
                            sigma[c] = sqrt(std::max(0.0, level_dispersion[c]
                                                          - level_mean[c] * level_mean[c]));
                        if (sigma_fixed[c] > 0)
                            sigma[c] = sigma_fixed[c];
                        level_threshold[c] = threshold * sigma[c];
                    }
                }

                memset(data_buffer, 0, VRT_SAMPLES_PER_PACKET * sizeof(uint32_t));

                for (uint32_t i = 0; i < num_samps; i++) {

                    const uint32_t word = buffer[vrt_packet.offset+i];
                    const int16_t  v[2] = {(int16_t)(word & 0xFFFF),
                                           (int16_t)((word >> 16) & 0xFFFF)};

                    for (int c = 0; c < 2; c++) {

                        uint32_t code;

                        if (bits == 1) {
                            // 2-level, 0 is negative
                            code = (v[c] >= 0) ? 1 : 0;
                        } else {
                            // 4-level, offset binary as used in VLBI:
                            // 00 = -high, 01 = -low, 10 = +low, 11 = +high
                            const double d = (double)v[c] - level_mean[c];
                            if (d >= 0)
                                code = (d > level_threshold[c]) ? 3 : 2;
                            else
                                code = (-d > level_threshold[c]) ? 0 : 1;
                        }

                        state_count[c][code]++;

                        // A field never straddles a word boundary, both 2 and 4
                        // bits per sample divide 32
                        const size_t bit_pos = i*bits_per_samp + c*bits;
                        data_buffer[bit_pos / 32] |= code << (bit_pos % 32);
                    }
                }

                state_samps += num_samps;

                new_data_len_words = (num_samps * bits_per_samp + 31) / 32;

            } else {

                const uint32_t bits_per_samp = 2 * unpack_bits;

                // Take the payload length from the received message rather than
                // from packet_size, which older versions left at the unquantized
                // length after shrinking the payload
                if (len/4 <= (int)vrt_packet.offset) {
                    fprintf(stderr, "# Error: packet of %d bytes has no payload\n", len);
                    continue;
                }

                const size_t payload_words = len/4 - vrt_packet.offset;
                num_samps     = payload_words * 32 / bits_per_samp;

                if (num_samps > VRT_SAMPLES_PER_PACKET) {
                    fprintf(stderr, "# Error: %zu samples in packet, truncating to %u\n",
                            num_samps, VRT_SAMPLES_PER_PACKET);
                    num_samps = VRT_SAMPLES_PER_PACKET;
                }

                const int16_t* levels = (unpack_bits == 1) ? levels_1_bit : levels_2_bit;
                const uint32_t mask   = (1u << unpack_bits) - 1;

                memset(data_buffer, 0, VRT_SAMPLES_PER_PACKET * sizeof(uint32_t));

                for (uint32_t i = 0; i < num_samps; i++) {

                    uint16_t v[2];

                    for (int c = 0; c < 2; c++) {
                        const size_t   bit_pos = i*bits_per_samp + c*unpack_bits;
                        const uint32_t code    =
                            (buffer[vrt_packet.offset + bit_pos/32] >> (bit_pos % 32)) & mask;
                        v[c] = (uint16_t)levels[code];
                    }

                    data_buffer[i] = ((uint32_t)v[1] << 16) | v[0];
                }

                // 16-bit complex, one word per sample
                new_data_len_words = num_samps;
            }

            memcpy((char*)&buffer[vrt_packet.offset], (char*)data_buffer,
                   new_data_len_words * sizeof(uint32_t));

            len = update_packet_size(buffer, vrt_packet.offset, new_data_len_words);
            if (len < 0)
                break;

            zmq_msg_t msg;
            zmq_msg_init_size (&msg, len);
            memcpy (zmq_msg_data(&msg), buffer, len);
            zmq_msg_send(&msg, responder, 0);
            zmq_msg_close(&msg);

            if (progress && vrt_packet.data && unpack)
            show_progress_stats(
                now,
                &last_update,
                &last_update_samps,
                &buffer[vrt_packet.offset],
                num_samps, 0
            );

            if (progress and not unpack and state_samps > 0
                and std::chrono::duration<double>(now - last_stats).count() > QUANTIZE_STATS_INTERVAL) {

                std::cout << "\t";
                for (int c = 0; c < 2; c++) {
                    std::cout << (c == 0 ? ". I: " : "  Q: ");
                    std::cout << boost::format("sigma %7.1f, dc %6.1f, states ") % sigma[c] % level_mean[c];
                    for (uint32_t s = 0; s < (1u << bits); s++)
                        std::cout << boost::format("%s%2.0f") % (s > 0 ? "/" : "")
                                         % (100.0*state_count[c][s]/state_samps);
                    std::cout << "%";
                }
                std::cout << std::endl;

                memset(state_count, 0, sizeof(state_count));
                state_samps = 0;
                last_stats  = now;
            }

        }
    }

    zmq_close(subscriber);
    zmq_close(responder);
    zmq_ctx_destroy(context);

    return 0;

}
