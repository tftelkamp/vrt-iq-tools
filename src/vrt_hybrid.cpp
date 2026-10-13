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

// #include <complex.h>
#include <complex>

#include "vrt-tools.h"
#include "dt-extended-context.h"
#include "tracker-extended-context.h"

namespace po = boost::program_options;

/* Interval between equalizer statistics lines, in seconds. */
#define HYBRID_STATS_INTERVAL 1.0

static bool stop_signal_called = false;
void sig_int_handler(int)
{
    stop_signal_called = true;
}

/* A data packet of one input channel, held until the packet of the other
 * input channel with the same timestamp has arrived. */
struct pending_packet {
    bool                  have;
    uint64_t              integer_seconds_timestamp;
    uint64_t              fractional_seconds_timestamp;
    uint32_t              offset;
    uint32_t              num_samps;
    std::vector<uint32_t> words;  // header, fields and payload
};

/* Smoothed level of one input channel. The power is the variance, so that a
 * DC offset of the receiver does not count as signal. */
struct level_estimate {
    double mean_re;
    double mean_im;
    double power;  // mean |x|^2, including DC
};

static void accumulate_level(const uint32_t* payload, uint32_t num_samps, level_estimate* l)
{
    double sum_re = 0, sum_im = 0, sum_pow = 0;

    for (uint32_t i = 0; i < num_samps; i++) {
        const uint32_t word = payload[i];
        const double   re   = (double)(int16_t)(word & 0xFFFF);
        const double   im   = (double)(int16_t)((word >> 16) & 0xFFFF);
        sum_re  += re;
        sum_im  += im;
        sum_pow += re * re + im * im;
    }

    l->mean_re = sum_re / num_samps;
    l->mean_im = sum_im / num_samps;
    l->power   = sum_pow / num_samps;
}

static double level_variance(const level_estimate* l)
{
    return std::max(0.0, l->power - l->mean_re * l->mean_re - l->mean_im * l->mean_im);
}

/* Delay a packet of samples by the length of hist, which holds the tail of the
 * previous packet and is updated to the tail of this one. */
static void delay_samples(const uint32_t* payload,
                          uint32_t num_samps,
                          std::vector<std::complex<float>>& hist,
                          std::vector<std::complex<float>>& scratch,
                          std::complex<float>* out)
{
    const size_t d = hist.size();

    scratch.resize(d + num_samps);
    std::copy(hist.begin(), hist.end(), scratch.begin());
    for (uint32_t i = 0; i < num_samps; i++) {
        const uint32_t word = payload[i];
        scratch[d + i] = std::complex<float>((float)(int16_t)(word & 0xFFFF),
                                             (float)(int16_t)((word >> 16) & 0xFFFF));
    }

    std::copy(scratch.begin(), scratch.begin() + num_samps, out);
    std::copy(scratch.end() - d, scratch.end(), hist.begin());
}

static inline int16_t saturate(float v, uint64_t* clipped)
{
    const float r = std::nearbyint(v);
    if (r > 32767.0f) {
        (*clipped)++;
        return 32767;
    }
    if (r < -32768.0f) {
        (*clipped)++;
        return -32768;
    }
    return (int16_t)r;
}

/* Send a copy of a packet with its stream id replaced, and for a data packet
 * its packet count, so that the output channels count in step. */
static void send_packet(void* responder,
                        const uint32_t* words,
                        int len,
                        uint32_t stream_id,
                        int32_t packet_count)
{
    zmq_msg_t msg;
    zmq_msg_init_size(&msg, len);
    uint32_t* out = (uint32_t*)zmq_msg_data(&msg);
    memcpy(out, words, len);

    if (packet_count >= 0) {
        uint32_t header = ntohl(out[0]);
        header          = (header & 0xFFF0FFFF) | ((uint32_t)packet_count % 16) << 16;
        out[0]          = htonl(header);
    }
    out[1] = htonl(stream_id);

    zmq_msg_send(&msg, responder, 0);
    zmq_msg_close(&msg);
}

int main(int argc, char* argv[])
{

    // variables to be set by po
    std::string zmq_address, channel_list;
    uint16_t pub_instance, instance, main_port, port, pub_port;
    int hwm;
    double angle, gain_db, tau;
    int32_t delay;

    size_t num_requested_samples;
    double total_time;

    // setup the program options
    po::options_description desc("Allowed options");
    // clang-format off

    desc.add_options()
        ("help", "help message")
        ("nsamps", po::value<size_t>(&num_requested_samples)->default_value(0), "total number of samples to receive")
        ("duration", po::value<double>(&total_time)->default_value(0), "total number of seconds to receive")
        ("channel", po::value<std::string>(&channel_list)->default_value("0,1"), "the two VRT channels to combine, e.g. H and V (specify \"0,1\")")
        ("angle", po::value<double>(&angle)->default_value(90.0), "phase rotation applied to the second channel before combining (degrees)")
        ("gain", po::value<double>(&gain_db), "fixed amplitude correction of the second channel relative to the first (dB), disables automatic equalization")
        ("tau", po::value<double>(&tau)->default_value(1.0), "time constant of the equalizer level estimator [seconds]")
        ("delay", po::value<int32_t>(&delay)->default_value(0), "delay of the second channel relative to the first (samples), negative delays the first")
        ("single", "only output the first combination (channel 0)")
        ("progress", "periodically display short-term bandwidth and equalizer levels")
        ("continue", "don't abort on a bad packet")
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
        std::cout << boost::format("VRT hybrid. %s") % desc << std::endl;
        std::cout << std::endl
                  << "This application is a software hybrid coupler: it coherently combines the\n"
                     "two channels A and B of a dual channel VRT stream, for example the H and V\n"
                     "channels of a dual linear feed, into\n"
                     "\n"
                     "    channel 0 = (gA*A + gB*exp(j*angle)*B) / sqrt(2)\n"
                     "    channel 1 = (gA*A - gB*exp(j*angle)*B) / sqrt(2)\n"
                     "\n"
                     "With A = H and B = V, --angle 90 makes channel 0 RHCP and channel 1 LHCP\n"
                     "(IEEE convention), --angle -90 the other way around. This assumes H points\n"
                     "right and V up as seen from the antenna looking at the source, and a\n"
                     "receiver that does not invert the spectrum. Any instrumental phase\n"
                     "difference between the channels is to be included in the angle.\n"
                     "\n"
                     "Before combining, the channels are equalized in amplitude: the gains gA and\n"
                     "gB bring both channels to the geometric mean of their powers, estimated\n"
                     "(without DC) over a sliding window with time constant --tau. This assumes\n"
                     "the power is dominated by noise or unpolarized emission; for a strongly\n"
                     "polarized source set a fixed correction with --gain instead.\n"
                     "\n"
                     "The input is expected to be a 16-bit stream with packets of both channels\n"
                     "sharing a timestamp, such as made by vrt_merge or a dual channel receiver.\n"
                  << std::endl;
        return ~0;
    }

    bool progress               = vm.count("progress") > 0;
    bool continue_on_bad_packet = vm.count("continue") > 0;
    bool single                 = vm.count("single") > 0;
    bool fixed_gain             = vm.count("gain") > 0;

    // detect which channels to use
    std::vector<std::string> channel_strings;
    std::vector<uint32_t> channel_nums;
    boost::split(channel_strings, channel_list, boost::is_any_of("\"',"));
    for (size_t ch = 0; ch < channel_strings.size(); ch++) {
        if (channel_strings[ch].empty())
            continue;
        channel_nums.push_back(std::stoi(channel_strings[ch]));
    }

    if (channel_nums.size() != 2 or channel_nums[0] == channel_nums[1]
        or channel_nums[0] >= 32 or channel_nums[1] >= 32) {
        std::cerr << "Specify two different channels with --channel, e.g. \"0,1\"." << std::endl;
        return ~0;
    }

    if (tau < 0) {
        std::cerr << "Error: --tau must not be negative." << std::endl;
        return ~0;
    }

    // Equalizer gains, fixed or estimated
    double gain[2] = {1.0, 1.0};
    if (fixed_gain) {
        gain[0] = pow(10.0, -gain_db / 40.0);
        gain[1] = pow(10.0, gain_db / 40.0);
    }

    // Combining rotation of the second channel
    const std::complex<double> rotation = std::polar(1.0, angle * M_PI / 180.0);

    context_type vrt_context;
    init_context(&vrt_context);

    context_type input_context[2];
    init_context(&input_context[0]);
    init_context(&input_context[1]);

    packet_type vrt_packet {};
    vrt_packet.channel_filt = (1u << channel_nums[0]) | (1u << channel_nums[1]);

    if (vm.count("port") > 0) {
        main_port = port;
    } else {
        main_port = DEFAULT_MAIN_PORT + MAX_CHANNELS*instance;
    }

    if (!(vm.count("pub-port") > 0)) {
        pub_port = DEFAULT_MAIN_PORT + MAX_CHANNELS*pub_instance;
    }

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
    uint32_t out_buffer[2][ZMQ_BUFFER_SIZE];

    uint64_t num_total_samps = 0;

    // Track time and samps between updating the BW summary
    auto last_update = start_time;
    uint64_t last_update_samps = 0;

    // Contexts of the two input channels, bit k set once channel k has one
    uint32_t contexts_received = 0;
    int32_t  last_gain[2]      = {0, 0};

    // Packets waiting for their counterpart on the other channel
    pending_packet pending[2];
    for (int k = 0; k < 2; k++) {
        pending[k].have = false;
        pending[k].words.resize(ZMQ_BUFFER_SIZE);
    }
    uint64_t unpaired = 0;
    bool unpaired_reported = false;

    // Equalizer state
    bool           levels_valid = false;
    level_estimate level[2]     = {{0, 0, 0}, {0, 0, 0}};
    auto           last_stats   = start_time;
    uint64_t       clipped      = 0;  // components of both combinations
    uint64_t       stats_samps  = 0;

    // Delay lines, one of which is empty
    std::vector<std::complex<float>> hist[2];
    hist[0].assign(delay < 0 ? -delay : 0, 0);
    hist[1].assign(delay > 0 ?  delay : 0, 0);
    std::vector<std::complex<float>> scratch;
    std::vector<std::complex<float>> samples[2];
    samples[0].resize(VRT_SAMPLES_PER_PACKET);
    samples[1].resize(VRT_SAMPLES_PER_PACKET);

    uint32_t frame_count = 0;

    // set to true to process data before context
    bool start_rx = false;

    int len;

    std::signal(SIGINT, &sig_int_handler);

    while (not stop_signal_called
           and (num_requested_samples > num_total_samps or num_requested_samples == 0)
           and (total_time == 0.0 or std::chrono::steady_clock::now() <= stop_time)) {

        len = zmq_recv(subscriber, buffer, ZMQ_BUFFER_SIZE, 0);
        if (len <= 0)
            continue;

        const auto now = std::chrono::steady_clock::now();

        if (not vrt_process(buffer, sizeof(buffer), &vrt_context, &vrt_packet)) {
            printf("Not a Vita49 packet?\n");
            continue;
        }

        if (not vrt_packet.context and not vrt_packet.data)
            continue;

        // Input channel 0 or 1
        const int k = (vrt_packet.stream_id & (1u << channel_nums[0])) ? 0 : 1;

        if (vrt_packet.context) {

            input_context[k] = vrt_context;
            contexts_received |= (1u << k);

            if (vrt_context.gain != last_gain[k]) {
                // A gain change invalidates the level estimate
                levels_valid = false;
                last_gain[k] = vrt_context.gain;
            }

            // Forward as the context of output channel k
            if (k == 0 or not single)
                send_packet(responder, buffer, len, 1u << k, -1);
        }

        if (not start_rx and contexts_received == 0x3) {

            for (int c = 0; c < 2; c++) {
                if (input_context[c].data_item_size != 16) {
                    fprintf(stderr, "# Error: channel %u carries %u-bit samples, "
                                    "vrt_hybrid requires 16-bit. Insert "
                                    "'vrt_quantize --unpack' ahead of it.\n",
                            channel_nums[c], input_context[c].data_item_size);
                    exit(1);
                }
            }

            if (input_context[0].sample_rate != input_context[1].sample_rate)
                fprintf(stderr, "# WARNING: sample rates differ (%u and %u samples per second)\n",
                        input_context[0].sample_rate, input_context[1].sample_rate);

            vrt_print_context(&input_context[0]);
            vrt_print_context(&input_context[1]);

            printf("# Hybrid parameters:\n");
            printf("#    Channels: %u, %u\n", channel_nums[0], channel_nums[1]);
            printf("#    Angle [deg]: %.2f\n", angle);
            printf("#    Delay [samples]: %d\n", delay);
            if (fixed_gain)
                printf("#    Gain [dB]: %.2f (fixed)\n", gain_db);
            else
                printf("#    Equalizer time constant [sec]: %.2f\n", tau);
            printf("#    Output: %s\n", single ? "channel 0" : "channels 0 and 1");

            start_rx = true;
        }

        if (not (start_rx and vrt_packet.data))
            continue;

        if (vrt_packet.lost_frame)
            if (not continue_on_bad_packet)
                break;

        if (vrt_packet.num_rx_samps == 0 or vrt_packet.num_rx_samps > VRT_SAMPLES_PER_PACKET) {
            fprintf(stderr, "# Error: packet with %u samples, skipping\n", vrt_packet.num_rx_samps);
            continue;
        }

        if (pending[k].have and frame_count > 0) {
            // Its counterpart never arrived
            unpaired++;
        }

        pending[k].have                         = true;
        pending[k].integer_seconds_timestamp    = vrt_packet.integer_seconds_timestamp;
        pending[k].fractional_seconds_timestamp = vrt_packet.fractional_seconds_timestamp;
        pending[k].offset                       = vrt_packet.offset;
        pending[k].num_samps                    = vrt_packet.num_rx_samps;
        memcpy(pending[k].words.data(), buffer, len);

        if (not (pending[0].have and pending[1].have))
            continue;

        {
            const struct vrt_time_ps t0 = {(int64_t)pending[0].integer_seconds_timestamp,
                                           pending[0].fractional_seconds_timestamp};
            const struct vrt_time_ps t1 = {(int64_t)pending[1].integer_seconds_timestamp,
                                           pending[1].fractional_seconds_timestamp};
            const int64_t skew = vrt_time_diff_ps(t0, t1);

            if (skew != 0 or pending[0].num_samps != pending[1].num_samps) {
                // Drop the older one, the newer may still find its counterpart
                pending[skew < 0 ? 0 : 1].have = false;
                if (skew == 0)
                    pending[0].have = pending[1].have = false;
                // Before the first output frame this is just the stream
                // starting halfway a pair
                if (frame_count > 0)
                    unpaired++;
                continue;
            }
        }

        if (unpaired > 0 and not unpaired_reported) {
            fprintf(stderr, "# WARNING: dropping packets without a counterpart on the other "
                            "channel, the output will have gaps\n");
            unpaired_reported = true;
        }

        const uint32_t num_samps = pending[0].num_samps;

        // Equalize
        if (not fixed_gain) {

            level_estimate packet_level[2];
            for (int c = 0; c < 2; c++)
                accumulate_level(&pending[c].words[pending[c].offset], num_samps, &packet_level[c]);

            if (not levels_valid) {
                level[0] = packet_level[0];
                level[1] = packet_level[1];
                levels_valid = true;
            } else {
                // Track the level slowly, so that the estimate follows gain
                // and system temperature drift but not the signal
                double alpha = 1.0;
                if (input_context[0].sample_rate > 0 and tau > 0)
                    alpha = 1.0 - exp(-((double)num_samps / input_context[0].sample_rate) / tau);
                for (int c = 0; c < 2; c++) {
                    level[c].mean_re += alpha * (packet_level[c].mean_re - level[c].mean_re);
                    level[c].mean_im += alpha * (packet_level[c].mean_im - level[c].mean_im);
                    level[c].power   += alpha * (packet_level[c].power   - level[c].power);
                }
            }

            const double p0 = level_variance(&level[0]);
            const double p1 = level_variance(&level[1]);
            if (p0 > 0 and p1 > 0) {
                // Both to the geometric mean of the two powers
                const double ratio = p0 / p1;
                gain[0] = pow(ratio, -0.25);
                gain[1] = pow(ratio,  0.25);
            }
        }

        for (int c = 0; c < 2; c++)
            delay_samples(&pending[c].words[pending[c].offset], num_samps, hist[c], scratch,
                          samples[c].data());

        // Combine
        {
            const float               g0  = (float)(gain[0] * M_SQRT1_2);
            const std::complex<float> rot = std::complex<float>(rotation * gain[1] * M_SQRT1_2);
            const uint32_t            off = pending[0].offset;

            for (uint32_t i = 0; i < num_samps; i++) {
                const std::complex<float> x = g0  * samples[0][i];
                const std::complex<float> y = rot * samples[1][i];
                const std::complex<float> s = x + y;
                const std::complex<float> d = x - y;

                out_buffer[0][off + i] =
                    ((uint32_t)(uint16_t)saturate(s.imag(), &clipped) << 16)
                    | (uint16_t)saturate(s.real(), &clipped);
                out_buffer[1][off + i] =
                    ((uint32_t)(uint16_t)saturate(d.imag(), &clipped) << 16)
                    | (uint16_t)saturate(d.real(), &clipped);
            }

            // Header and fields of the first channel, timestamps are equal
            memcpy(out_buffer[0], pending[0].words.data(), off * sizeof(uint32_t));
            memcpy(out_buffer[1], pending[0].words.data(), off * sizeof(uint32_t));

            const int out_len = (int)(off + num_samps) * 4;
            send_packet(responder, out_buffer[0], out_len, 1u << 0, frame_count);
            if (not single)
                send_packet(responder, out_buffer[1], out_len, 1u << 1, frame_count);
            frame_count++;
        }

        pending[0].have = false;
        pending[1].have = false;

        num_total_samps += num_samps;
        stats_samps     += num_samps;

        if (progress)
            show_progress_stats(
                now,
                &last_update,
                &last_update_samps,
                &out_buffer[0][pending[0].offset],
                num_samps, 0
            );

        if (progress and stats_samps > 0
            and std::chrono::duration<double>(now - last_stats).count() > HYBRID_STATS_INTERVAL) {

            const double fs = 32768.0 * 32768.0;
            std::cout << "\t";
            if (not fixed_gain)
                std::cout << boost::format("CH%u: %6.1f dBFS, CH%u: %6.1f dBFS, ")
                                 % channel_nums[0] % (10*log10(level_variance(&level[0])/fs))
                                 % channel_nums[1] % (10*log10(level_variance(&level[1])/fs));
            std::cout << boost::format("gain %+.2f dB, %.2f%% clip.")
                             % (20*log10(gain[1]/gain[0]))
                             % (100.0*clipped/(4*stats_samps));
            std::cout << std::endl;

            clipped     = 0;
            stats_samps = 0;
            last_stats  = now;
        }
    }

    if (unpaired > 0)
        fprintf(stderr, "# %llu packets dropped without a counterpart on the other channel\n",
                (long long unsigned int)unpaired);

    zmq_close(subscriber);
    zmq_close(responder);
    zmq_ctx_destroy(context);

    return 0;

}
