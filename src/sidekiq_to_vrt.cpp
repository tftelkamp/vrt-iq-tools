//
// Copyright 2026 by Thomas Telkamp
//
// SPDX-License-Identifier: MIT
//
// Stream IQ samples from an Epiq Sidekiq SDR to VRT (Vita49) over ZeroMQ.
//
// Timestamps: the Sidekiq system timestamp (sys_freq ticks per second) is set
// to Unix time at startup, optionally on the external 1PPS edge. The first
// received block anchors the stream to absolute time; every packet timestamp
// after that is derived from the sample count (the per-block RF timestamp,
// which advances by one per sample) relative to that anchor. The stream is
// re-anchored after a discontinuity in the RF timestamp, which happens when
// samples are lost and after retuning (AD9361 timestamp slip).
//
// PPS (--pps): the clock is set on a PPS edge and verified on the next one.
// While streaming, the card latches the system and RF timestamps at every
// edge; the offset of the stream timestamps at the last edge is reported in
// the context packets (timestamp_adjustment). With --pps-track the timestamps
// themselves follow the PPS: every edge re-anchors them (sample at the edge =
// that second) and the sample rate is corrected by the clock error measured
// against the PPS, so they don't drift with the card's oscillator.
//

#include <dirent.h>
#include <pthread.h>
#include <sched.h>
#include <sys/resource.h>
#include <sys/syscall.h>
#include <unistd.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cinttypes>
#include <cmath>
#include <complex>
#include <csignal>
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

// Epiq
#include "sidekiq_api.h"

#include <zmq.h>

// VRT
#include <vrt/vrt_init.h>
#include <vrt/vrt_read.h>
#include <vrt/vrt_string.h>
#include <vrt/vrt_types.h>
#include <vrt/vrt_util.h>
#include <vrt/vrt_write.h>

#include <boost/format.hpp>
#include <boost/program_options.hpp>

// VRT tools functions
#include "vrt-tools.h"

#define DEFAULT_CARD_NUMBER 0
#define NUM_USEC_IN_MS (1000)

// TX block size in samples (words)
#define TX_BLOCK_SIZE 1020

// Sidekiq ADC samples are 12 bit, sign extended to 16 bit (not shifted by
// libsidekiq). They are scaled to 16 bit full scale in the VRT stream, like
// usrp_to_vrt, and the TX path scales back.
#define SIDEKIQ_FULL_SCALE 2048.0
#define SIDEKIQ_SAMPLE_SCALE 16

// Gaps in the sample stream up to this many samples are filled with zeros;
// larger ones re-anchor the timestamps.
#define MAX_ZERO_FILL VRT_SAMPLES_PER_PACKET

// After a retune the RF timestamp can slip for a few more blocks; within this
// many blocks a discontinuity re-anchors without being reported as sample loss.
#define RETUNE_SETTLE_BLOCKS 16

// A jump of more than this (s) between the time from the RF timestamp and the
// time from the system timestamp of a block is a timestamp slip.
#define MAX_TIMESTAMP_SLIP 1e-3L

namespace po = boost::program_options;

static std::atomic<bool> stop_signal_called{false};

// After a critical error libsidekiq must not be called anymore (by default it
// exits the process); this handler stops the tool cleanly instead.
static std::atomic<bool> critical_error{false};
static void critical_error_handler(int32_t status, void*)
{
    fprintf(stderr, "Error: libsidekiq critical error (status %d), stopping\n", (int)status);
    critical_error     = true;
    stop_signal_called = true;
}

static void release_card()
{
    if (not critical_error)
        skiq_exit();
}

// SIGINT/SIGTERM are blocked in all threads and taken by a dedicated thread,
// so they don't interrupt system calls in libsidekiq (which reports an
// interrupted USB poll as a receive error).
static void start_signal_thread()
{
    sigset_t stop_signals;
    sigemptyset(&stop_signals);
    sigaddset(&stop_signals, SIGINT);
    sigaddset(&stop_signals, SIGTERM);
    pthread_sigmask(SIG_BLOCK, &stop_signals, NULL);
    std::thread([stop_signals]() {
        int sig;
        sigwait(&stop_signals, &sig);
        stop_signal_called = true;
    }).detach();
}

// Convert a Sidekiq system timestamp (ticks of sys_freq Hz since the Unix
// epoch, as set at startup) to VRT time.
// 12 bit two's complement to int16
static inline int16_t sign_extend_12(uint32_t v)
{
    return (int16_t)((int16_t)(uint16_t)(v << 4) >> 4);
}

// Packed I/Q (--packed, skiq_write_iq_pack_mode): every 3 little-endian 32 bit
// words hold 4 samples as 8 fields of 12 bits, from the most significant bit
// down: word 0 a0 b0 a1[11:4], word 1 a1[3:0] b1 a2 b2[11:8], word 2 b2[7:0]
// a3 b3 (sidekiq_api.h calls a I and b Q). With skiq_iq_order_iq the card
// swaps the fields too, so a is Q and b is I (measured on the TX2i, as in
// ls-sdr-dev). Unpacked into I/Q int16 pairs, the same as the card sends them
// unpacked (12 bits, sign extended). num_samples is a multiple of 4.
static void unpack_iq(const void* packed, int16_t* iq, uint32_t num_samples)
{
    const uint8_t* p = (const uint8_t*)packed;
    for (uint32_t n = 0; n < num_samples; n += 4, p += 12, iq += 8) {
        uint32_t w0, w1, w2;
        memcpy(&w0, p, 4);
        memcpy(&w1, p + 4, 4);
        memcpy(&w2, p + 8, 4);
        iq[0] = sign_extend_12(w0 >> 8);                // I0 (b0)
        iq[1] = sign_extend_12(w0 >> 20);               // Q0 (a0)
        iq[2] = sign_extend_12(w1 >> 16);               // I1 (b1)
        iq[3] = sign_extend_12((w0 << 4) | (w1 >> 28)); // Q1 (a1)
        iq[4] = sign_extend_12((w1 << 8) | (w2 >> 24)); // I2 (b2)
        iq[5] = sign_extend_12(w1 >> 4);                // Q2 (a2)
        iq[6] = sign_extend_12(w2);                     // I3 (b3)
        iq[7] = sign_extend_12(w2 >> 12);               // Q3 (a3)
    }
}

static struct vrt_time_ps sys_ticks_to_time(uint64_t ticks, uint64_t sys_freq)
{
    struct vrt_time_ps t;
    t.seconds = (int64_t)(ticks / sys_freq);
    t.frac_ps = (uint64_t)(((unsigned __int128)(ticks % sys_freq) * VRT_PS_PER_SECOND) / sys_freq);
    return t;
}

// Convert VRT (Unix) time to Sidekiq system timestamp ticks.
static uint64_t time_to_sys_ticks(const struct vrt_time_ps& t, uint64_t sys_freq)
{
    return (uint64_t)t.seconds * sys_freq
           + (uint64_t)(((unsigned __int128)t.frac_ps * sys_freq) / VRT_PS_PER_SECOND);
}

// t + seconds (which may be negative), at picosecond resolution
static struct vrt_time_ps time_add_seconds(struct vrt_time_ps t, long double seconds)
{
    const long double ps_per_s = (long double)VRT_PS_PER_SECOND;
    const long double frac     = (long double)t.frac_ps + seconds * ps_per_s;
    const long double whole    = floorl(frac / ps_per_s);
    t.seconds += (int64_t)whole;
    t.frac_ps = (uint64_t)llroundl(frac - whole * ps_per_s);
    vrt_time_normalize(&t);
    return t;
}

// a - b in seconds
static long double time_diff(const struct vrt_time_ps& a, const struct vrt_time_ps& b)
{
    return (long double)(a.seconds - b.seconds)
           + ((long double)a.frac_ps - (long double)b.frac_ps) / (long double)VRT_PS_PER_SECOND;
}

static double seconds_since(std::chrono::steady_clock::time_point t)
{
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count();
}

/******************************** scheduling *********************************/

// Parse a CPU list like "1,2" or "1-2" (as used by taskset/isolcpus).
static bool parse_cpu_list(const std::string& list, std::vector<int>& cpus)
{
    std::stringstream ss(list);
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (item.empty())
            return false;
        char* end;
        long first = strtol(item.c_str(), &end, 10);
        long last  = first;
        if (*end == '-')
            last = strtol(end + 1, &end, 10);
        if (*end != '\0' or first < 0 or last < first or last >= CPU_SETSIZE)
            return false;
        for (long c = first; c <= last; c++)
            cpus.push_back((int)c);
    }
    return not cpus.empty();
}

static bool set_thread_affinity(pid_t tid, const std::vector<int>& cpus)
{
    cpu_set_t set;
    CPU_ZERO(&set);
    for (int c : cpus)
        CPU_SET(c, &set);
    return sched_setaffinity(tid, sizeof(set), &set) == 0;
}

// Realtime (SCHED_RR) scheduling for the calling thread only. It is set on the
// receive thread once streaming runs, so libsidekiq's own (USB) threads keep
// normal scheduling: on tx2i giving those threads realtime priority made sample
// loss worse, not better. Lowers the priority to RLIMIT_RTPRIO when that is
// smaller (non-root users). Returns the priority set, or 0 if realtime
// scheduling is not permitted.
static int set_realtime_priority(int priority)
{
    struct rlimit limit;
    if (geteuid() != 0 and getrlimit(RLIMIT_RTPRIO, &limit) == 0
        and limit.rlim_cur != RLIM_INFINITY and (rlim_t)priority > limit.rlim_cur) {
        if (limit.rlim_cur == 0)
            return 0;
        fprintf(stderr, "Warning: realtime priority limited to %d by RLIMIT_RTPRIO\n",
                (int)limit.rlim_cur);
        priority = (int)limit.rlim_cur;
    }
    struct sched_param param;
    param.sched_priority = priority;
    if (pthread_setschedparam(pthread_self(), SCHED_RR, &param) != 0)
        return 0;
    return priority;
}

// Pin the main (receive) thread to the first CPU and all other threads of the
// process to the remaining CPUs. On isolated cores (isolcpus) the kernel does
// not move threads between CPUs, so this spreads the load explicitly.
static void spread_threads(pid_t main_tid, const std::vector<int>& cpus)
{
    if (cpus.size() < 2)
        return;
    std::vector<int> main_cpu(cpus.begin(), cpus.begin() + 1);
    std::vector<int> other_cpus(cpus.begin() + 1, cpus.end());
    if (not set_thread_affinity(main_tid, main_cpu))
        fprintf(stderr, "Warning: unable to pin receive thread to CPU %d: %s\n", main_cpu[0],
                strerror(errno));
    int moved  = 0;
    DIR* tasks = opendir("/proc/self/task");
    if (tasks == NULL)
        return;
    while (struct dirent* entry = readdir(tasks)) {
        pid_t tid = (pid_t)atoi(entry->d_name);
        if (tid <= 0 or tid == main_tid)
            continue;
        if (set_thread_affinity(tid, other_cpus))
            moved++;
    }
    closedir(tasks);
    std::string others;
    for (int c : other_cpus)
        others += (others.empty() ? "" : ",") + std::to_string(c);
    printf("Info: receive thread on CPU %d, %d other threads on CPU %s\n", main_cpu[0], moved,
           others.c_str());
}

// Receive buffers; large enough for any message we accept (ZMQ_BUFFER_SIZE
// bytes). Static, so they don't live on the stack.
static uint32_t control_buffer[ZMQ_BUFFER_SIZE / sizeof(uint32_t)];
static uint32_t transmit_buffer[ZMQ_BUFFER_SIZE / sizeof(uint32_t)];

void transmit_worker(uint8_t card, skiq_tx_hdl_t curr_tx_hdl, void* zmq_transmit)
{
    skiq_tx_block_t* p_tx_block = skiq_tx_block_allocate(TX_BLOCK_SIZE);
    if (p_tx_block == NULL) {
        fprintf(stderr, "Error: unable to allocate Tx block\n");
        return;
    }

    if (skiq_start_tx_streaming(card, curr_tx_hdl) != 0) {
        fprintf(stderr, "Error: unable to start Tx streaming\n");
        skiq_tx_block_free(p_tx_block);
        return;
    }

    uint32_t counter = 0;

    // send data until the signal handler gets called; the socket has a receive
    // timeout so the stop flag is checked regularly
    while (not stop_signal_called) {

        int len = zmq_recv(zmq_transmit, transmit_buffer, sizeof(transmit_buffer), 0);
        if (len <= 0)
            continue;
        if ((size_t)len > sizeof(transmit_buffer)) {
            fprintf(stderr, "Warning: dropping oversized Tx message (%d bytes)\n", len);
            continue;
        }

        const int32_t size = len / sizeof(uint32_t);

        struct vrt_header h;
        struct vrt_fields f;

        int32_t offset = 0;
        int32_t rv = vrt_read_header(transmit_buffer + offset, size - offset, &h, true);
        if (rv < 0) {
            fprintf(stderr, "Failed to parse Tx header: %s\n", vrt_string_error(rv));
            continue;
        }
        offset += rv;

        if (h.packet_type != VRT_PT_IF_DATA_WITH_STREAM_ID)
            continue;

        rv = vrt_read_fields(&h, transmit_buffer + offset, size - offset, &f, true);
        if (rv < 0) {
            fprintf(stderr, "Failed to parse Tx fields section: %s\n", vrt_string_error(rv));
            continue;
        }
        offset += rv;

        int32_t num_tx_samps = (int32_t)h.packet_size - offset - (h.has.trailer ? 1 : 0);
        num_tx_samps = std::min(num_tx_samps, size - offset);
        if (num_tx_samps <= 0)
            continue;

        for (int32_t i = 0; i < num_tx_samps; i++) {
            int16_t re;
            memcpy(&re, (char*)&transmit_buffer[offset + i], 2);
            int16_t img;
            memcpy(&img, (char*)&transmit_buffer[offset + i] + 2, 2);

            // 16 bit VRT samples to the 12 bit DAC
            p_tx_block->data[counter * 2]     = re / SIDEKIQ_SAMPLE_SCALE;
            p_tx_block->data[counter * 2 + 1] = img / SIDEKIQ_SAMPLE_SCALE;

            counter++;

            if (counter == TX_BLOCK_SIZE) {
                int32_t status = skiq_transmit(card, curr_tx_hdl, p_tx_block, NULL);
                if (0 != status) {
                    fprintf(stderr, "Error: failed to transmit data (result code %" PRIi32 ")\n",
                            status);
                }
                counter = 0;
            }
        }
    }

    skiq_stop_tx_streaming(card, curr_tx_hdl);
    skiq_tx_block_free(p_tx_block);
}

int main(int argc, char* argv[])
{
    // variables to be set by po
    std::string ref, merge_address, cpu_list;
    int rt_priority;
    size_t total_num_samps;
    uint16_t port, instance, merge_port, card_arg;
    uint16_t tx_gain;
    uint16_t rx_gain_arg = 0;
    int hwm;
    double rate, bw, freq, total_time, setup_time, tx_freq, if_freq, pps_offset;
    bool merge;

    std::string args;

    // setup the program options
    po::options_description desc("Allowed options");
    // clang-format off
    desc.add_options()
        ("help", "help message")
        ("args", po::value<std::string>(&args)->default_value(""), "device args (unused)")
        ("card", po::value<uint16_t>(&card_arg)->default_value(DEFAULT_CARD_NUMBER), "Sidekiq card number")
        ("nsamps", po::value<size_t>(&total_num_samps)->default_value(0), "total number of samples to receive")
        ("duration", po::value<double>(&total_time)->default_value(0), "total number of seconds to receive")
        ("rate", po::value<double>(&rate)->default_value(1e6), "rate of incoming samples")
        ("freq", po::value<double>(&freq)->default_value(0.0), "RF center frequency in Hz")
        ("if-freq", po::value<double>(&if_freq)->default_value(0.0), "IF center frequency in Hz")
        ("gain", po::value<uint16_t>(&rx_gain_arg), "gain index for the RF chain (automatic gain if not set)")
        ("bw", po::value<double>(&bw),
            "channel bandwidth in Hz, RX and TX (RF IC FIR passband, analog filter set to match; default 0.8 x rate)")
        ("ref", po::value<std::string>(&ref)->default_value("internal"), "reference source (internal, external)")
        ("tx", "enable tx")
        ("tx-freq", po::value<double>(&tx_freq)->default_value(0.0), "TX RF center frequency in Hz")
        ("tx-gain", po::value<uint16_t>(&tx_gain)->default_value(0), "TX attenuation (in 0.25 dB steps)")
        ("setup", po::value<double>(&setup_time)->default_value(1.0), "seconds of setup time")
        ("progress", "periodically display short-term bandwidth")
        ("stats", "show average bandwidth on exit")
        ("pps", "synchronize the clock to the external PPS signal (exit if there is none)")
        ("pps-offset", po::value<double>(&pps_offset)->default_value(0),
            "time of the PPS edge relative to the integer second, in seconds")
        ("pps-track", "derive the timestamps from the PPS edges instead of the card clock (requires --pps)")
        ("temp", "read temperature sensor")
        ("int-second", "start streaming on a PPS edge (requires --pps)")
        ("null", "run without streaming")
        ("continue", "don't abort on a bad packet")
        ("no-dc-offset", "disable the FPGA RX DC offset correction (the RF IC's own DC tracking stays on)")
        ("packed", "transfer the 12 bit samples packed (4 samples in 3 words, 25% less data over USB)")
        ("port", po::value<uint16_t>(&port), "VRT ZMQ port")
        ("instance", po::value<uint16_t>(&instance)->default_value(0), "VRT ZMQ instance")
        ("merge", po::value<bool>(&merge)->default_value(true), "Merge another VRT ZMQ stream (SUB connect)")
        ("merge-port", po::value<uint16_t>(&merge_port)->default_value(50011), "VRT ZMQ merge port")
        ("merge-address", po::value<std::string>(&merge_address)->default_value("localhost"), "VRT ZMQ merge address")
        ("hwm", po::value<int>(&hwm)->default_value(10000), "VRT ZMQ HWM")
        ("priority", po::value<int>(&rt_priority)->implicit_value(50),
            "realtime (SCHED_RR) scheduling of the receive thread, optionally with priority 1-99 (default 50)")
        ("cpu", po::value<std::string>(&cpu_list),
            "pin to CPUs, e.g. 1,2 or 1-2: receive thread on the first, other threads on the rest")
    ;
    // clang-format on
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // print the help message
    if (vm.count("help")) {
        std::cout << boost::format("Epiq Sidekiq samples to Vita49 over ZMQ %s") % desc << std::endl;
        std::cout << std::endl
                  << "This application streams data from an Epiq Sidekiq "
                     "device to Vita49 over ZMQ.\n"
                  << std::endl;
        return ~0;
    }

    bool bw_summary  = vm.count("progress") > 0;
    bool stats       = vm.count("stats") > 0;
    bool null        = vm.count("null") > 0;
    bool enable_temp = vm.count("temp") > 0;
    bool enable_tx   = vm.count("tx") > 0;
    bool use_pps     = vm.count("pps") > 0;
    bool pps_track   = vm.count("pps-track") > 0;
    bool int_second  = vm.count("int-second") > 0;
    bool manual_gain = vm.count("gain") > 0;
    bool dc_offset   = vm.count("no-dc-offset") == 0;
    bool packed      = vm.count("packed") > 0;

    if (card_arg >= SKIQ_MAX_NUM_CARDS) {
        fprintf(stderr, "Error: card number %u out of range (0-%u)\n", card_arg,
                SKIQ_MAX_NUM_CARDS - 1);
        return EXIT_FAILURE;
    }

    // The pack mode applies to TX as well, and the TX path sends unpacked samples
    if (packed and enable_tx) {
        fprintf(stderr, "Error: --packed can't be combined with --tx\n");
        return EXIT_FAILURE;
    }

    if ((pps_track or int_second) and not use_pps) {
        fprintf(stderr, "Error: --pps-track and --int-second require --pps\n");
        return EXIT_FAILURE;
    }
    if (pps_offset <= -1.0 or pps_offset >= 1.0) {
        fprintf(stderr, "Error: --pps-offset must be within (-1, 1) seconds\n");
        return EXIT_FAILURE;
    }

    if (rate <= 0) {
        fprintf(stderr, "Error: invalid sample rate %f\n", rate);
        return EXIT_FAILURE;
    }
    if (ref != "internal" and ref != "external") {
        fprintf(stderr, "Error: unsupported reference source '%s' (internal or external)\n",
                ref.c_str());
        return EXIT_FAILURE;
    }
    if (manual_gain and rx_gain_arg > 255) {
        fprintf(stderr, "Error: gain index %u out of range\n", rx_gain_arg);
        return EXIT_FAILURE;
    }
    uint8_t rx_gain = (uint8_t)rx_gain_arg;

    std::vector<int> cpus;
    if (vm.count("cpu") and not parse_cpu_list(cpu_list, cpus)) {
        fprintf(stderr, "Error: invalid CPU list '%s'\n", cpu_list.c_str());
        return EXIT_FAILURE;
    }
    if (vm.count("priority") and (rt_priority < 1 or rt_priority > 99)) {
        fprintf(stderr, "Error: realtime priority %d out of range (1-99)\n", rt_priority);
        return EXIT_FAILURE;
    }

    start_signal_thread();

    /******************************** scheduling *********************************/

    // CPU affinity is set before the card is initialized, so libsidekiq's
    // threads inherit it; realtime priority is set once streaming runs.
    const pid_t main_tid = (pid_t)syscall(SYS_gettid);
    if (not cpus.empty()) {
        if (set_thread_affinity(main_tid, cpus)) {
            printf("Info: pinned to CPU %s\n", cpu_list.c_str());
        } else {
            fprintf(stderr, "Error: unable to pin to CPU %s: %s\n", cpu_list.c_str(),
                    strerror(errno));
            return EXIT_FAILURE;
        }
    }

    /******************************* Sidekiq setup *******************************/

    uint8_t card              = (uint8_t)card_arg;
    skiq_rx_hdl_t curr_rx_hdl = skiq_rx_hdl_A1;
    skiq_tx_hdl_t curr_tx_hdl = skiq_tx_hdl_A1;

    skiq_rx_stream_mode_t stream_mode = skiq_rx_stream_mode_high_tput;

    uint32_t read_sample_rate;
    double actual_sample_rate;
    uint32_t read_bandwidth;
    uint32_t actual_bandwidth;

    int32_t status = 0;
    pid_t owner    = 0;

    const auto startup_time = std::chrono::steady_clock::now();
    auto step_time          = startup_time;

    // Report how long each configuration step took (card initialization can
    // take many seconds).
    auto step_done = [&](const char* what) {
        printf("Info: %s (%.2f s)\n", what, seconds_since(step_time));
        step_time = std::chrono::steady_clock::now();
    };

    skiq_register_critical_error_callback(critical_error_handler, NULL);

    printf("Info: initializing card %" PRIu8 "...\n", card);
    status = skiq_init(skiq_xport_type_auto, skiq_xport_init_level_full, &card, 1);
    if (status != 0) {
        if ((EBUSY == status) and (0 != skiq_is_card_avail(card, &owner))) {
            fprintf(stderr,
                    "Error: card %" PRIu8 " is already in use (by process ID %u); cannot initialize card.\n",
                    card, (unsigned int)owner);
        } else if (-EINVAL == status) {
            fprintf(stderr,
                    "Error: unable to initialize libsidekiq; was a valid card specified? (result code %" PRIi32 ")\n",
                    status);
        } else {
            fprintf(stderr, "Error: unable to initialize libsidekiq with status %" PRIi32 "\n", status);
        }
        return EXIT_FAILURE;
    }
    step_done("card initialized");

    // From here on, release the card on errors
    auto fail = [&](const char* what, int32_t rc) {
        fprintf(stderr, "Error: %s (status %" PRIi32 ": %s)\n", what, rc, strerror(abs(rc)));
        release_card();
        return EXIT_FAILURE;
    };

    // Only write the reference clock when external is requested; not all
    // Sidekiq products support changing it (e.g. mPCIe).
    if (ref == "external") {
        printf("Info: setting reference clock to external\n");
        status = skiq_write_ref_clock_select(card, skiq_ref_clock_external);
        if (status != 0)
            return fail("failed to select external reference clock", status);
    }

    status = skiq_write_iq_order_mode(card, skiq_iq_order_iq);
    if (status != 0)
        return fail("failed to set I/Q order mode", status);

    status = skiq_write_iq_pack_mode(card, packed);
    if (status != 0)
        return fail(packed ? "failed to enable packed I/Q" : "failed to disable packed I/Q", status);
    printf("Info: I/Q transfer: %s\n", packed ? "packed 12 bit samples" : "16 bit samples");

    status = skiq_write_chan_mode(card, skiq_chan_mode_single);
    if (status != 0)
        return fail("failed to set Rx channel mode", status);

    uint32_t sample_rate = (uint32_t)rate;
    uint32_t bandwidth   = (vm.count("bw") > 0) ? (uint32_t)bw : (uint32_t)(0.8 * rate);

    status = skiq_write_rx_sample_rate_and_bandwidth(card, curr_rx_hdl, sample_rate, bandwidth);
    if (status != 0)
        return fail("failed to set Rx sample rate or bandwidth", status);

    status = skiq_read_rx_sample_rate_and_bandwidth(card, curr_rx_hdl, &read_sample_rate,
                                                    &actual_sample_rate, &read_bandwidth,
                                                    &actual_bandwidth);
    if (status != 0)
        return fail("failed to read Rx sample rate and bandwidth", status);

    printf("Info: requested sample rate is %u, requested bandwidth is %u\n", sample_rate, bandwidth);
    printf("Info: actual sample rate is %f, actual bandwidth is %u\n", actual_sample_rate,
           actual_bandwidth);

    if ((actual_sample_rate != sample_rate) || (actual_bandwidth < bandwidth)) {
        printf("Warning: sample rate or bandwidth does not match the requested value.\n");
    }
    step_done("sample rate configured");

    if (enable_tx) {
        double actual_tx_sample_rate;
        uint32_t actual_tx_bandwidth;
        status = skiq_write_tx_sample_rate_and_bandwidth(card, curr_tx_hdl, sample_rate, bandwidth);
        if (status != 0)
            return fail("failed to set Tx sample rate or bandwidth", status);

        status = skiq_read_tx_sample_rate_and_bandwidth(card, curr_tx_hdl, &read_sample_rate,
                                                        &actual_tx_sample_rate, &read_bandwidth,
                                                        &actual_tx_bandwidth);
        if (status == 0) {
            printf("Info: actual Tx sample rate is %f, actual Tx bandwidth is %u\n",
                   actual_tx_sample_rate, actual_tx_bandwidth);
            if ((actual_tx_sample_rate != sample_rate) || (actual_tx_bandwidth < bandwidth))
                printf("Warning: Tx sample rate or bandwidth does not match the requested value.\n");
        }
    }

    /* tune the Rx chain to the requested freq */
    status = skiq_write_rx_LO_freq(card, curr_rx_hdl, (uint64_t)freq);
    if (status != 0)
        return fail("failed to set Rx LO frequency", status);
    printf("Info: configured Rx LO freq to %" PRIu64 " Hz\n", (uint64_t)freq);

    if (enable_tx) {
        status = skiq_write_tx_LO_freq(card, curr_tx_hdl, (uint64_t)tx_freq);
        if (status != 0)
            return fail("failed to set Tx LO frequency", status);
        printf("Info: configured Tx LO freq to %" PRIu64 " Hz\n", (uint64_t)tx_freq);
    }
    step_done("LO configured");

    /* now that the Rx freq is set, set the gain mode and gain */
    status = skiq_write_rx_gain_mode(card, curr_rx_hdl,
                                     manual_gain ? skiq_rx_gain_manual : skiq_rx_gain_auto);
    if (status != 0)
        return fail("failed to set Rx gain mode", status);
    printf("Info: configured %s gain mode\n", manual_gain ? "manual" : "auto");

    if (manual_gain) {
        status = skiq_write_rx_gain(card, curr_rx_hdl, rx_gain);
        if (status != 0)
            return fail("failed to set Rx gain index", status);
        printf("Info: set gain index to %" PRIu8 "\n", rx_gain);
    }

    status = skiq_write_rx_dc_offset_corr(card, curr_rx_hdl, dc_offset);
    if (status != 0)
        return fail("failed to configure the Rx DC offset correction", status);
    if (not dc_offset)
        printf("Info: FPGA DC offset correction disabled\n");

    if (enable_tx) {
        if ((status = skiq_write_tx_attenuation(card, curr_tx_hdl, tx_gain)) != 0)
            return fail("unable to configure Tx attenuation", status);
        if ((status = skiq_write_tx_data_flow_mode(card, curr_tx_hdl,
                                                   skiq_tx_immediate_data_flow_mode)) != 0)
            return fail("unable to configure Tx data flow mode", status);
        if ((status = skiq_write_tx_block_size(card, curr_tx_hdl, TX_BLOCK_SIZE)) != 0)
            return fail("unable to configure Tx block size", status);
        if ((status = skiq_write_tx_transfer_mode(card, curr_tx_hdl, skiq_tx_transfer_mode_sync)) != 0)
            return fail("unable to configure Tx transfer mode", status);
    }

    int32_t block_size_in_words = skiq_read_rx_block_size(card, stream_mode) / 4;
    if (block_size_in_words <= 0)
        return fail("failed to read Rx block size", block_size_in_words);
    const int32_t payload_words = block_size_in_words - (int32_t)SKIQ_RX_HEADER_SIZE_IN_WORDS;
    const uint32_t block_samples =
        packed ? (uint32_t)SKIQ_NUM_PACKED_SAMPLES_IN_BLOCK(payload_words) : (uint32_t)payload_words;
    printf("Info: acquiring at %" PRIu32 " samples per block\n", block_samples);

    // packed blocks are unpacked into this buffer (I/Q int16 pairs)
    std::vector<int16_t> unpacked_iq(packed ? 2 * (size_t)block_samples : 0);

    // skiq_receive() waits at most this long (us) for a block, so the main
    // loop keeps servicing control/merge sockets without spinning.
    status = skiq_set_rx_transfer_timeout(card, 10000);
    if (status != 0)
        return fail("unable to set Rx transfer timeout", status);

    uint64_t sys_freq = 0;
    status            = skiq_read_sys_timestamp_freq(card, &sys_freq);
    if (status != 0 or sys_freq == 0)
        return fail("unable to read system timestamp frequency", status);
    printf("Info: system timestamp frequency %" PRIu64 " Hz\n", sys_freq);

    /***************************** time synchronization *****************************/

    // Set the system timestamp to the host clock (Unix time)
    status = skiq_update_timestamps(card, time_to_sys_ticks(vrt_time_now(), sys_freq));
    if (status != 0)
        return fail("unable to set timestamps", status);

    // PPS edges: the card latches the system and RF timestamps at every edge.
    uint32_t pps_calibration_second = 0; // second at which the clock was set on PPS
    uint64_t sync_edge_rf = 0, sync_edge_sys = 0;

    // Wait (up to timeout_s) for a PPS edge other than the one at prev_sys.
    auto wait_for_pps_edge = [&](uint64_t prev_sys, double timeout_s, uint64_t& edge_rf,
                                 uint64_t& edge_sys) {
        const auto t0 = std::chrono::steady_clock::now();
        while (not stop_signal_called and seconds_since(t0) < timeout_s) {
            if (skiq_read_last_1pps_timestamp(card, &edge_rf, &edge_sys) == 0
                and edge_sys != prev_sys)
                return true;
            std::this_thread::sleep_for(std::chrono::milliseconds(5));
        }
        return false;
    };

    if (use_pps) {
        status = skiq_write_1pps_source(card, skiq_1pps_source_external);
        if (status != 0)
            return fail("unable to select external 1PPS source", status);

        uint64_t edge_rf = 0, edge_sys = 0;
        skiq_read_last_1pps_timestamp(card, &edge_rf, &edge_sys);
        printf("Info: waiting for PPS...\n");
        if (not wait_for_pps_edge(edge_sys, 2.5, edge_rf, edge_sys)) {
            fprintf(stderr, "Error: no PPS signal on the external 1PPS input\n");
            release_card();
            return EXIT_FAILURE;
        }

        // Just after an edge, arm the card to load the time of the next edge
        // (taken from the host clock, which only needs to be within 0.5 s),
        // then verify on the edge after that.
        const int max_attempts = 5;
        bool synced            = false;
        for (int attempt = 1; attempt <= max_attempts and not synced and not stop_signal_called;
             attempt++) {
            uint64_t s0, s1, s2, rf;
            if (not wait_for_pps_edge(edge_sys, 1.5, rf, s0))
                break;
            const struct vrt_time_ps host = vrt_time_now();
            const int64_t this_second =
                llroundl(time_diff(host, {0, 0}) - (long double)pps_offset);
            const int64_t next_second = this_second + 1;
            const uint64_t next_ticks =
                (uint64_t)llroundl(((long double)next_second + pps_offset) * (long double)sys_freq);

            status = skiq_write_timestamp_update_on_1pps(card, s0 + sys_freq / 2, next_ticks);
            if (status != 0)
                return fail("unable to set timestamp on 1PPS", status);

            if (not wait_for_pps_edge(s0, 1.5, rf, s1) or not wait_for_pps_edge(s1, 1.5, rf, s2))
                break;
            edge_sys = s2;

            // the verification edge must be exactly one second after the
            // update edge (within 20 ppm), and the card clock within 0.5 s
            // of the host clock
            const int64_t err_ticks = (int64_t)(s2 - (next_ticks + sys_freq));
            uint64_t card_now;
            skiq_read_curr_sys_timestamp(card, &card_now);
            const long double host_diff =
                time_diff(sys_ticks_to_time(card_now, sys_freq), vrt_time_now());
            if (std::llabs(err_ticks) < (int64_t)(sys_freq / 50000) and fabsl(host_diff) < 0.5L) {
                synced                 = true;
                pps_calibration_second = (uint32_t)next_second;
                sync_edge_rf           = rf;
                sync_edge_sys          = s2;
                printf("Info: clock set on PPS edge of second %" PRIi64
                       ", verified on the next edge (%+" PRIi64 " ticks, card - host %+.3f ms)\n",
                       next_second, err_ticks, (double)host_diff * 1e3);
            } else {
                fprintf(stderr,
                        "Warning: PPS sync attempt %d failed (%+" PRIi64
                        " ticks on verification edge, card - host %+.3f s), retrying\n",
                        attempt, err_ticks, (double)host_diff);
            }
        }
        if (not synced) {
            fprintf(stderr, "Error: unable to synchronize to the PPS signal\n");
            release_card();
            return EXIT_FAILURE;
        }
    }

    uint64_t epiq_timestamp;
    skiq_read_curr_sys_timestamp(card, &epiq_timestamp);
    struct vrt_time_ps host_now = vrt_time_now();
    std::cout << boost::format("Epiq Device time right now: %.6f seconds\n")
                     % ((double)epiq_timestamp / (double)sys_freq);
    std::cout << boost::format("PC Clock time:              %.6f seconds\n")
                     % ((double)host_now.seconds + (double)host_now.frac_ps / 1e12);
    step_done("card configured");

    /*********************************** ZeroMQ ***********************************/

    // Sockets are set up only now, once the card is ready, so clients don't
    // connect to a stream that has not started yet.
    void* context = zmq_ctx_new();

    auto zmq_fail = [&](const char* what, const std::string& endpoint) {
        fprintf(stderr, "Error: unable to %s %s: %s\n", what, endpoint.c_str(), zmq_strerror(errno));
        release_card();
        return EXIT_FAILURE;
    };

    // Like the other tools: the data port follows --port or --instance, the
    // control and TX ports are at fixed offsets from it.
    const uint16_t main_port =
        vm.count("port") > 0 ? port : (uint16_t)(DEFAULT_MAIN_PORT + MAX_CHANNELS * instance);
    printf("Info: VRT ZMQ data port %u, control port %u%s\n", main_port, main_port + 200,
           enable_tx ? (", TX port " + std::to_string(main_port + 400)).c_str() : "");

    void* zmq_server   = zmq_socket(context, ZMQ_PUB);
    std::string endpoint = "tcp://*:" + std::to_string(main_port);
    zmq_setsockopt(zmq_server, ZMQ_SNDHWM, &hwm, sizeof hwm);
    if (zmq_bind(zmq_server, endpoint.c_str()) != 0)
        return zmq_fail("bind", endpoint);

    void* zmq_control = zmq_socket(context, ZMQ_SUB);
    endpoint          = "tcp://*:" + std::to_string(main_port + 200);
    if (zmq_bind(zmq_control, endpoint.c_str()) != 0)
        return zmq_fail("bind", endpoint);
    zmq_setsockopt(zmq_control, ZMQ_SUBSCRIBE, "", 0);

    void* zmq_transmit = NULL;
    if (enable_tx) {
        zmq_transmit = zmq_socket(context, ZMQ_SUB);
        endpoint     = "tcp://*:" + std::to_string(main_port + 400);
        if (zmq_bind(zmq_transmit, endpoint.c_str()) != 0)
            return zmq_fail("bind", endpoint);
        zmq_setsockopt(zmq_transmit, ZMQ_SUBSCRIBE, "", 0);
        int timeout_ms = 100;
        zmq_setsockopt(zmq_transmit, ZMQ_RCVTIMEO, &timeout_ms, sizeof timeout_ms);
    }

    void* merge_zmq = NULL;
    if (merge) {
        merge_zmq = zmq_socket(context, ZMQ_SUB);
        endpoint  = "tcp://" + merge_address + ":" + std::to_string(merge_port);
        if (zmq_connect(merge_zmq, endpoint.c_str()) != 0)
            return zmq_fail("connect to", endpoint);
        zmq_setsockopt(merge_zmq, ZMQ_SUBSCRIBE, "", 0);
    }

    /*********************************** VRT ***********************************/

    struct vrt_packet p;
    vrt_init_packet(&p);
    vrt_init_data_packet(&p);
    // Only 1 channel
    p.fields.stream_id = 1;

    // packet assembly buffer
    static int16_t packet_samples[2 * VRT_SAMPLES_PER_PACKET];
    uint32_t packet_fill     = 0; // samples in packet_samples
    uint64_t packet_first_rf = 0; // RF timestamp of the first sample in packet_samples

    // timestamp anchor: absolute time of the sample with RF timestamp anchor_rf
    bool anchored                  = false;
    bool resync_pending            = false; // re-anchor after retune (expected slip)
    uint32_t resync_blocks_left    = 0;     // blocks in which the slip may still settle
    uint64_t anchor_rf             = 0;
    struct vrt_time_ps anchor_time = {0, 0};
    uint64_t expected_rf           = 0;
    uint64_t last_block_rf         = 0;
    uint64_t last_block_sys        = 0;
    // samples per second of stream time: the nominal rate, or with --pps-track
    // the rate measured against the PPS
    double stream_rate = actual_sample_rate;

    // PPS monitoring (--pps)
    const int PPS_HISTORY = 16;          // edges used to measure the clock rate
    struct pps_edge {
        uint64_t sys;
        int64_t second;
    } pps_history[PPS_HISTORY];
    int pps_history_len              = 0;
    bool pps_have_edge               = false; // last accepted edge below is valid
    uint64_t pps_edge_sys            = 0;
    uint64_t pps_edge_rf             = 0;
    uint64_t pps_stale_rf            = 0;     // edges with an unusable RF timestamp
    int64_t pps_edge_second          = 0;
    uint64_t pps_last_read_sys       = 0;     // last edge read (accepted or not)
    bool pps_present                 = use_pps;
    uint64_t pps_slip_sys            = 0;     // RF timestamps latched before this are stale
    double clock_ppm                 = 0;     // card clock rate error against the PPS
    bool clock_ppm_valid             = false;
    bool pps_tracking                = false; // timestamps follow the PPS (--pps-track)
    bool anchor_from_sys             = true;  // anchor taken from a block system timestamp
    long double sys_anchor_bias      = 0;     // its measured time error (s)
    const long double PPS_TRACK_TAU       = 30.0;  // s, to steer out a residual
    const long double PPS_TRACK_MAX_SLEW  = 10e-6; // s, larger residuals are stepped
    bool have_timestamp_adjustment   = false;
    int64_t timestamp_adjustment_ps  = 0;     // stream time - PPS time at the last edge
    uint64_t pps_edges = 0, pps_glitches = 0, pps_losses = 0;
    int64_t max_abs_adjustment_ps    = 0;

    uint32_t frame_count = 0;

    // statistics
    unsigned long long num_total_samps = 0;
    uint64_t num_gaps                  = 0;
    uint64_t num_lost_samps            = 0;
    uint64_t num_overruns              = 0;
    uint64_t num_overload_blocks       = 0;
    bool overload_since_context        = false;
    bool sample_loss_since_context     = false;
    bool context_changed               = true;

    // progress statistics
    uint64_t progress_samps   = 0;
    int16_t progress_peak     = 0;
    double progress_sum_sq    = 0;
    uint64_t progress_clipped = 0;

    unsigned long long num_requested_samples = total_num_samps;
    if (total_time > 0)
        num_requested_samples = (unsigned long long)(total_time * rate);

    // flush merge queue
    if (merge) {
        zmq_msg_t msg;
        zmq_msg_init(&msg);
        while (zmq_msg_recv(&msg, merge_zmq, ZMQ_DONTWAIT) >= 0) {
        }
        zmq_msg_close(&msg);
    }

    // TX
    std::thread transmit_thread;
    if (enable_tx) {
        transmit_thread = std::thread([&]() { transmit_worker(card, curr_tx_hdl, zmq_transmit); });
    }

    // RX
    uint64_t start_edge_sys = 0;
    if (int_second) {
        // start on the first PPS edge at least half a second after the last
        // one; this call blocks until that edge
        uint64_t edge_rf = 0;
        skiq_read_last_1pps_timestamp(card, &edge_rf, &start_edge_sys);
        printf("Info: start streaming on next PPS\n");
        status = skiq_start_rx_streaming_on_1pps(card, curr_rx_hdl, start_edge_sys + sys_freq / 2);
    } else {
        status = skiq_start_rx_streaming(card, curr_rx_hdl);
    }
    if (status != 0) {
        fprintf(stderr, "Error: receive streaming failed to start with status code %" PRIi32 "\n",
                status);
        stop_signal_called = true;
    } else {
        printf("Info: streaming started %.2f s after startup\n", seconds_since(startup_time));
        // libsidekiq starts its receive threads with streaming
        spread_threads(main_tid, cpus);
        if (vm.count("priority")) {
            int prio = set_realtime_priority(rt_priority);
            if (prio > 0)
                printf("Info: realtime scheduling enabled (SCHED_RR, priority %d)\n", prio);
            else
                fprintf(stderr, "Warning: realtime scheduling not permitted (check RLIMIT_RTPRIO, "
                                "/etc/security/limits.d); continuing without\n");
        }
        std::cout << "Press Ctrl + C to stop streaming..." << std::endl;
    }

    const auto start_time = std::chrono::steady_clock::now();
    auto last_update      = start_time;
    auto last_context     = start_time - std::chrono::milliseconds(2 * VRT_CONTEXT_INTERVAL);

    auto packet_time = [&](uint64_t rf) {
        const int64_t delta = (int64_t)(rf - anchor_rf);
        if (stream_rate == actual_sample_rate and delta >= 0)
            return vrt_time_add_samples(anchor_time, (uint64_t)delta, actual_sample_rate);
        return time_add_seconds(anchor_time, (long double)delta / stream_rate);
    };

    // Time of the first sample of a block with system timestamp sys: the card
    // clock as is, or with --pps-track mapped through the last PPS edge and
    // the measured clock rate, corrected for the measured offset between the
    // block system timestamps and the RF timestamps latched at PPS edges.
    auto sys_time = [&](uint64_t sys) {
        if (pps_tracking and pps_have_edge) {
            const long double ticks_per_s = (long double)sys_freq * (1.0L + clock_ppm * 1e-6L);
            return time_add_seconds(
                time_add_seconds({pps_edge_second, 0}, (long double)pps_offset),
                (long double)(int64_t)(sys - pps_edge_sys) / ticks_per_s - sys_anchor_bias);
        }
        return sys_ticks_to_time(sys, sys_freq);
    };

    // sample rate of the stream (in samples per PPS second) without steering
    auto pps_base_rate = [&]() {
        return actual_sample_rate * (clock_ppm_valid ? 1.0 + clock_ppm * 1e-6 : 1.0);
    };

    // Process a PPS edge latched at system timestamp sys and RF timestamp rf.
    auto handle_pps_edge = [&](uint64_t rf, uint64_t sys) {
        pps_last_read_sys = sys;

        int64_t second;
        int64_t k = 0; // seconds since the previous edge
        if (pps_have_edge) {
            // must be a whole number of seconds after the previous edge
            const int64_t dsys = (int64_t)(sys - pps_edge_sys);
            k                  = llroundl((long double)dsys / (long double)sys_freq);
            const int64_t err  = dsys - k * (int64_t)sys_freq;
            if (k < 1 or std::llabs(err) > k * (int64_t)(sys_freq / 50000)) {
                pps_glitches++;
                fprintf(stderr, "Warning: ignoring PPS edge %.6f s after the previous one\n",
                        (double)dsys / (double)sys_freq);
                return;
            }
            second = pps_edge_second + k;
        } else {
            second = llroundl((long double)sys / (long double)sys_freq - (long double)pps_offset);
        }
        pps_edges++;

        if (not pps_present) {
            printf("Info: PPS signal restored\n");
            pps_present = true;
        }

        // clock rate over the last PPS_HISTORY edges
        if (pps_history_len == PPS_HISTORY) {
            std::move(pps_history + 1, pps_history + PPS_HISTORY, pps_history);
            pps_history_len--;
        }
        pps_history[pps_history_len++] = {sys, second};
        if (pps_history_len >= 2) {
            const pps_edge& first = pps_history[0];
            const long double ticks_per_s =
                (long double)(sys - first.sys) / (long double)(second - first.second);
            clock_ppm       = (double)((ticks_per_s / (long double)sys_freq - 1.0L) * 1e6L);
            clock_ppm_valid = true;
        }

        const struct vrt_time_ps edge_time =
            time_add_seconds({second, 0}, (long double)pps_offset);

        // The RF timestamp is latched in the sample clock domain: a short or
        // noisy PPS pulse can be missed there while the system timestamp
        // (40 MHz) still catches it, leaving a stale RF timestamp. Only use
        // it if it advanced by the expected number of samples since the
        // previous edge (which must be after any timestamp slip).
        bool rf_latch_ok = pps_have_edge and pps_edge_sys > pps_slip_sys;
        if (rf_latch_ok) {
            const long double expected = (long double)k * (long double)actual_sample_rate
                                         * (1.0L + (clock_ppm_valid ? clock_ppm * 1e-6L : 0.0L));
            const long double drf      = (long double)(int64_t)(rf - pps_edge_rf);
            if (fabsl(drf - expected) > 2.0L + expected * 20e-6L) {
                rf_latch_ok = false;
                if (pps_stale_rf++ == 0)
                    fprintf(stderr,
                            "Warning: RF timestamp latched at a PPS edge advanced by %.0f samples "
                            "instead of %.0f; not using it (PPS pulse too short or noisy?)\n",
                            (double)drf, (double)expected);
            }
        }

        // The RF timestamp latched at the edge is only usable if it is from
        // the current RF timestamp epoch (no retune since).
        if (rf_latch_ok and anchored and not resync_pending and resync_blocks_left == 0
            and sys > pps_slip_sys) {
            // stream time - PPS time at the edge (positive: ahead)
            // The latched RF timestamp is the sample period the edge fell in;
            // on average the edge is in the middle of it.
            const struct vrt_time_ps model_time =
                time_add_seconds(packet_time(rf), 0.5L / (long double)stream_rate);
            const long double residual = time_diff(model_time, edge_time);
            timestamp_adjustment_ps   = (int64_t)llroundl(residual * (long double)VRT_PS_PER_SECOND);
            have_timestamp_adjustment = true;
            max_abs_adjustment_ps =
                std::max(max_abs_adjustment_ps, (int64_t)std::llabs(timestamp_adjustment_ps));

            if (pps_track) {
                // the time error of an anchor taken from a block system timestamp
                if (anchor_from_sys) {
                    sys_anchor_bias += residual;
                    anchor_from_sys = false;
                }
                // The RF timestamp latched at an edge has a resolution of one
                // sample, so instead of re-anchoring on every edge (steps of up
                // to a sample) the residual is removed by steering the rate
                // over PPS_TRACK_TAU seconds; only a large error is stepped.
                const long double max_residual =
                    std::max(PPS_TRACK_MAX_SLEW, 4.0L / (long double)actual_sample_rate);
                long double remaining = residual;
                if (not pps_tracking or fabsl(residual) > max_residual) {
                    printf("Info: timestamps %s the PPS (step of %+.3f us)\n",
                           pps_tracking ? "re-aligned to" : "now follow",
                           (double)(-residual * 1e6L));
                    anchor_time = time_add_seconds(edge_time, -0.5L / (long double)stream_rate);
                    remaining   = 0;
                } else {
                    anchor_time = packet_time(rf); // continuous
                }
                anchor_rf    = rf;
                stream_rate  = pps_base_rate() * (double)(1.0L + remaining / PPS_TRACK_TAU);
                pps_tracking = true;
            }
        }

        pps_edge_sys    = sys;
        pps_edge_rf     = rf;
        pps_edge_second = second;
        pps_have_edge   = true;
    };

    // the edge the clock was verified on
    if (use_pps)
        handle_pps_edge(sync_edge_rf, sync_edge_sys);
    auto last_pps_poll = start_time;
    auto last_pps_seen = start_time;

    auto send_data_packet = [&]() {
        struct vrt_time_ps t = packet_time(packet_first_rf);

        p.body                                = packet_samples;
        p.header.packet_count                 = (uint8_t)(frame_count % 16);
        p.fields.integer_seconds_timestamp    = (uint32_t)t.seconds;
        p.fields.fractional_seconds_timestamp = t.frac_ps;

        if (not null) {
            zmq_msg_t msg;
            zmq_msg_init_size(&msg, VRT_DATA_PACKET_SIZE * 4);
            int32_t rv = vrt_write_packet(&p, zmq_msg_data(&msg), VRT_DATA_PACKET_SIZE, true);
            if (rv < 0) {
                fprintf(stderr, "Failed to write packet: %s\n", vrt_string_error(rv));
                zmq_msg_close(&msg);
            } else {
                zmq_msg_send(&msg, zmq_server, 0);
            }
        }

        frame_count++;
        num_total_samps += VRT_SAMPLES_PER_PACKET;
    };

    // Add samples (or zeros if src is NULL) starting at RF timestamp rf to
    // the packet being assembled, sending packets as they fill up.
    auto append_samples = [&](const int16_t* src, uint32_t count, uint64_t rf) {
        uint32_t done = 0;
        while (done < count) {
            if (packet_fill == 0)
                packet_first_rf = rf + done;
            const uint32_t n =
                std::min(count - done, (uint32_t)VRT_SAMPLES_PER_PACKET - packet_fill);
            if (src) {
                // 12 bit ADC samples (sign extended) to 16 bit full scale
                int16_t* dst       = &packet_samples[2 * packet_fill];
                const int16_t* in  = &src[2 * done];
                for (uint32_t i = 0; i < 2 * n; i++)
                    dst[i] = (int16_t)(in[i] * SIDEKIQ_SAMPLE_SCALE);
            } else
                memset(&packet_samples[2 * packet_fill], 0, n * 2 * sizeof(int16_t));
            packet_fill += n;
            done += n;
            if (packet_fill == VRT_SAMPLES_PER_PACKET) {
                send_data_packet();
                packet_fill = 0;
            }
        }
    };

    // sample loss is reported at most once per second
    uint64_t gaps_since_report = 0;
    uint64_t lost_since_report = 0;
    auto last_loss_report      = start_time;

    bool first_frame = true;

    while (not stop_signal_called
           and (num_requested_samples > num_total_samps or num_requested_samples == 0)) {

        skiq_rx_block_t* p_rx_block;
        uint32_t len;

        skiq_rx_status_t rx_status = skiq_receive(card, &curr_rx_hdl, &p_rx_block, &len);

        if (rx_status == skiq_rx_status_success) {

            // len is in bytes; unpacked one word per sample, packed 4 samples in 3 words
            uint32_t num_rx_samps        = len / 4 - SKIQ_RX_HEADER_SIZE_IN_WORDS;
            const uint64_t rf_timestamp  = p_rx_block->rf_timestamp;
            const uint64_t sys_timestamp = p_rx_block->sys_timestamp;
            const int16_t* samples       = (const int16_t*)(const void*)p_rx_block->data;
            if (packed) {
                num_rx_samps = std::min((uint32_t)SKIQ_NUM_PACKED_SAMPLES_IN_BLOCK(num_rx_samps),
                                        block_samples);
                unpack_iq((const void*)p_rx_block->data, unpacked_iq.data(), num_rx_samps);
                samples = unpacked_iq.data();
            }

            // Over-range: the RF overload flag (only delivered over PCIe; the
            // USB transport zeroes it) or samples at ADC full scale.
            bool over_range = p_rx_block->overload;
            for (uint32_t i = 0; i < 2 * num_rx_samps and not over_range; i++)
                over_range = samples[i] >= SIDEKIQ_FULL_SCALE - 1 or samples[i] <= -SIDEKIQ_FULL_SCALE;
            if (over_range) {
                num_overload_blocks++;
                overload_since_context = true;
            }

            bool reanchor = not anchored or resync_pending;

            if (anchored and not resync_pending and rf_timestamp != expected_rf
                and resync_blocks_left > 0) {
                // the RF timestamp is still settling after a retune
                packet_fill = 0;
                reanchor    = true;
            } else if (anchored and not resync_pending and rf_timestamp != expected_rf) {
                // unexpected discontinuity: samples were lost
                const int64_t missing = (int64_t)(rf_timestamp - expected_rf);
                num_gaps++;
                gaps_since_report++;
                sample_loss_since_context = true;

                if (missing > 0 and missing <= MAX_ZERO_FILL) {
                    // Small gap: fill with zeros, so the stream (and its
                    // timestamps) stays continuous.
                    num_lost_samps += missing;
                    lost_since_report += missing;
                    append_samples(NULL, (uint32_t)missing, expected_rf);
                } else {
                    // Large gap or timestamp jump: drop the partial packet,
                    // skip a packet count so receivers see a lost frame, and
                    // re-anchor.
                    fprintf(stderr,
                            "Warning: RF timestamp jump of %" PRIi64 " samples (RF timestamp %" PRIu64
                            ", expected %" PRIu64 "), re-anchoring timestamps\n",
                            missing, rf_timestamp, expected_rf);
                    if (missing > 0) {
                        num_lost_samps += missing;
                        lost_since_report += missing;
                    }
                    if (packet_fill > 0) {
                        packet_fill = 0;
                        frame_count++;
                    }
                    reanchor = true;
                }
            }

            // A retune pauses the RF timestamp counter together with the
            // samples, so it does not always show as a jump in the RF
            // timestamp. It does show as a jump between the time derived from
            // the RF timestamp and the time of the block system timestamp,
            // which runs from the same oscillator but never pauses.
            if (anchored and not reanchor) {
                const long double slip =
                    time_diff(sys_time(sys_timestamp), packet_time(rf_timestamp));
                if (fabsl(slip) > MAX_TIMESTAMP_SLIP) {
                    printf("Info: RF timestamp slipped by %+.3f ms, re-anchoring timestamps\n",
                           (double)(slip * 1e3L));
                    if (packet_fill > 0) {
                        packet_fill = 0;
                        frame_count++;
                    }
                    pps_slip_sys = sys_timestamp; // RF timestamps latched before are stale
                    reanchor     = true;
                }
            }

            if (resync_blocks_left > 0)
                resync_blocks_left--;

            if (reanchor) {
                if (resync_pending) {
                    // don't mix samples from before and after a retune
                    packet_fill = 0;
                    printf("Info: re-anchored timestamps after retune\n");
                    resync_blocks_left = RETUNE_SETTLE_BLOCKS;
                }

                anchor_rf      = rf_timestamp;
                anchor_time    = sys_time(sys_timestamp);
                anchor_from_sys = true;
                anchored       = true;
                resync_pending = false;
            }

            expected_rf    = rf_timestamp + num_rx_samps;
            last_block_rf  = rf_timestamp;
            last_block_sys = sys_timestamp;

            if (first_frame) {
                struct vrt_time_ps t = packet_time(rf_timestamp);
                std::cout << boost::format("First frame: %u samples, %u full secs, %.09f frac secs")
                                 % num_rx_samps % t.seconds % ((double)t.frac_ps / 1e12)
                          << std::endl;
                if (int_second) {
                    // the first sample should be the one at the start edge
                    uint64_t edge_rf, edge_sys;
                    if (skiq_read_last_1pps_timestamp(card, &edge_rf, &edge_sys) == 0)
                        printf("Info: first sample is %" PRIi64 " samples after the PPS edge"
                               " (%s)\n",
                               (int64_t)(rf_timestamp - edge_rf),
                               edge_sys > start_edge_sys ? "start edge" : "edge not yet seen");
                }
                first_frame = false;
                last_update = std::chrono::steady_clock::now();
            }

            if (bw_summary) {
                progress_samps += num_rx_samps;
                for (uint32_t i = 0; i < 2 * num_rx_samps; i++) {
                    const int16_t v = samples[i] < 0 ? -samples[i] : samples[i];
                    progress_peak   = std::max(progress_peak, v);
                    progress_sum_sq += (double)samples[i] * samples[i];
                    if (v >= SIDEKIQ_FULL_SCALE - 1)
                        progress_clipped++;
                }
            }

            append_samples(samples, num_rx_samps, rf_timestamp);

        } else if (rx_status == skiq_rx_status_error_overrun) {
            // data was overwritten before it was read; the gap shows up in
            // the RF timestamp of the next block
            num_overruns++;
            fprintf(stderr, "Warning: Rx overrun, host does not consume data fast enough\n");
        } else if (rx_status != skiq_rx_status_no_data) {
            fprintf(stderr, "Error: skiq_receive failed with status %d\n", (int)rx_status);
        }

        const auto now = std::chrono::steady_clock::now();

        if (gaps_since_report > 0 and now - last_loss_report > std::chrono::seconds(1)) {
            fprintf(stderr,
                    "Warning: %" PRIu64 " sample discontinuities, %" PRIu64
                    " samples lost (zero-filled) in the last %.1f s\n",
                    gaps_since_report, lost_since_report,
                    std::chrono::duration<double>(now - last_loss_report).count());
            gaps_since_report = 0;
            lost_since_report = 0;
            last_loss_report  = now;
        }

        // PPS: poll for new edges (the card latches the timestamps at each
        // one). Reading them takes milliseconds over USB, so only poll once
        // the next edge is due (by the card clock), then every 100 ms.
        const bool pps_edge_due = not pps_have_edge or last_block_sys == 0
                                  or last_block_sys > pps_edge_sys + sys_freq + sys_freq / 100;
        if (use_pps and pps_edge_due and now - last_pps_poll > std::chrono::milliseconds(100)) {
            last_pps_poll = now;
            uint64_t edge_rf, edge_sys;
            if (skiq_read_last_1pps_timestamp(card, &edge_rf, &edge_sys) == 0
                and edge_sys != pps_last_read_sys) {
                last_pps_seen = now;
                handle_pps_edge(edge_rf, edge_sys);
            }
            if (pps_present and now - last_pps_seen > std::chrono::milliseconds(1500)) {
                pps_present = false;
                pps_losses++;
                if (pps_tracking and anchored) {
                    // stop steering; continue at the measured clock rate
                    anchor_time = packet_time(last_block_rf);
                    anchor_rf   = last_block_rf;
                    stream_rate = pps_base_rate();
                }
                fprintf(stderr, "Warning: PPS signal lost%s\n",
                        pps_tracking ? ", timestamps continue at the last measured rate" : "");
            }
        }

        // Context
        if (anchored and now - last_context > std::chrono::milliseconds(VRT_CONTEXT_INTERVAL)) {

            last_context = now;

            struct vrt_packet pc;
            vrt_init_packet(&pc);
            vrt_init_context_packet(&pc);

            // time of the most recent block
            struct vrt_time_ps t                   = packet_time(last_block_rf);
            pc.fields.integer_seconds_timestamp    = (uint32_t)t.seconds;
            pc.fields.fractional_seconds_timestamp = t.frac_ps;

            if (enable_temp) {
                int8_t temp_in_deg_C;
                if (skiq_read_temp(card, &temp_in_deg_C) == 0) {
                    pc.if_context.has.temperature = true;
                    pc.if_context.temperature     = (float)temp_in_deg_C;
                }
            }

            uint8_t gain_index = rx_gain;
            if (not manual_gain)
                skiq_read_rx_gain(card, curr_rx_hdl, &gain_index);

            pc.if_context.context_field_change_indicator = context_changed;

            pc.fields.stream_id                         = p.fields.stream_id;
            pc.if_context.bandwidth                     = actual_bandwidth;
            pc.if_context.sample_rate                   = actual_sample_rate;
            pc.if_context.rf_reference_frequency        = freq;
            pc.if_context.rf_reference_frequency_offset = 0;
            pc.if_context.if_reference_frequency        = if_freq; // 0 for Zero-IF
            pc.if_context.if_band_offset                = 0;
            pc.if_context.gain.stage1                   = gain_index;
            pc.if_context.gain.stage2                   = 0;

            pc.if_context.state_and_event_indicators.has.reference_lock = true;
            pc.if_context.state_and_event_indicators.reference_lock     = (ref == "external");

            // time is calibrated while the PPS sync is verified and edges
            // keep arriving
            pc.if_context.state_and_event_indicators.has.calibrated_time = true;
            pc.if_context.state_and_event_indicators.calibrated_time     = use_pps and pps_present;

            if (pps_calibration_second != 0) {
                pc.if_context.has.timestamp_calibration_time = true;
                pc.if_context.timestamp_calibration_time     = pps_calibration_second;
            }
            // how far the stream timestamps were off at the last PPS edge
            // (positive: ahead), like usrp_to_vrt
            if (have_timestamp_adjustment and pps_present) {
                pc.if_context.has.timestamp_adjustment = true;
                pc.if_context.timestamp_adjustment     = timestamp_adjustment_ps;
            }

            pc.if_context.state_and_event_indicators.has.agc_or_mgc = true;
            pc.if_context.state_and_event_indicators.agc_or_mgc =
                manual_gain ? VRT_AOM_MGC : VRT_AOM_AGC;

            pc.if_context.state_and_event_indicators.has.over_range = true;
            pc.if_context.state_and_event_indicators.over_range     = overload_since_context;

            pc.if_context.state_and_event_indicators.has.sample_loss = true;
            pc.if_context.state_and_event_indicators.sample_loss     = sample_loss_since_context;

            static uint32_t context_buffer[VRT_DATA_PACKET_SIZE];
            int32_t rv = vrt_write_packet(&pc, context_buffer, VRT_DATA_PACKET_SIZE, true);
            if (rv < 0) {
                fprintf(stderr, "Failed to write context packet: %s\n", vrt_string_error(rv));
            } else {
                zmq_send(zmq_server, context_buffer, rv * 4, 0);
            }

            context_changed           = false;
            overload_since_context    = false;
            sample_loss_since_context = false;
        }

        // Merge: forward messages as-is
        if (merge) {
            zmq_msg_t msg;
            zmq_msg_init(&msg);
            while (zmq_msg_recv(&msg, merge_zmq, ZMQ_DONTWAIT) >= 0) {
                if (zmq_msg_send(&msg, zmq_server, 0) < 0)
                    zmq_msg_close(&msg);
                zmq_msg_init(&msg);
            }
            zmq_msg_close(&msg);
        }

        // Control
        int clen = zmq_recv(zmq_control, control_buffer, sizeof(control_buffer), ZMQ_DONTWAIT);
        if (clen > 0) {
            if ((size_t)clen > sizeof(control_buffer)) {
                fprintf(stderr, "Warning: ignoring oversized control message (%d bytes)\n", clen);
            } else {
                printf("-> Control context received\n");

                const int32_t size = clen / sizeof(uint32_t);
                struct vrt_header h;
                struct vrt_fields f;
                struct vrt_if_context c;

                int32_t offset = 0;
                int32_t rv     = vrt_read_header(control_buffer, size, &h, true);
                if (rv < 0) {
                    fprintf(stderr, "Failed to parse control header: %s\n", vrt_string_error(rv));
                } else if (h.packet_type == VRT_PT_IF_CONTEXT) {
                    offset += rv;
                    rv = vrt_read_fields(&h, control_buffer + offset, size - offset, &f, true);
                    if (rv < 0) {
                        fprintf(stderr, "Failed to parse control fields section: %s\n",
                                vrt_string_error(rv));
                    } else {
                        offset += rv;
                        rv = vrt_read_if_context(control_buffer + offset, size - offset, &c, true);
                        if (rv < 0) {
                            fprintf(stderr, "Failed to parse control IF context section: %s\n",
                                    vrt_string_error(rv));
                        } else {
                            if (c.has.rf_reference_frequency) {
                                double new_freq = (double)round(c.rf_reference_frequency);
                                std::cout << boost::format("    Setting RX Freq: %f MHz...")
                                                 % (new_freq / 1e6)
                                          << std::endl;
                                status = skiq_write_rx_LO_freq(card, curr_rx_hdl, (uint64_t)new_freq);
                                if (status != 0) {
                                    fprintf(stderr, "Error: failed to set Rx LO freq (status %" PRIi32 ")\n",
                                            status);
                                } else {
                                    freq = new_freq;
                                }
                                // retuning slips the RF timestamp, also the
                                // one latched at PPS edges before now
                                skiq_read_curr_sys_timestamp(card, &pps_slip_sys);
                                resync_pending = true;
                                context_changed = true;
                            }
                            if (c.has.gain) {
                                if (c.gain.stage1 < 0 or c.gain.stage1 > 255) {
                                    fprintf(stderr, "Error: gain index %f out of range\n", c.gain.stage1);
                                } else {
                                    uint8_t new_gain = (uint8_t)c.gain.stage1;
                                    if (not manual_gain)
                                        skiq_write_rx_gain_mode(card, curr_rx_hdl, skiq_rx_gain_manual);
                                    status = skiq_write_rx_gain(card, curr_rx_hdl, new_gain);
                                    if (status != 0) {
                                        fprintf(stderr, "Error: failed to set gain index to %" PRIu8 "\n",
                                                new_gain);
                                    } else {
                                        manual_gain = true;
                                        rx_gain     = new_gain;
                                        printf("Info: set gain index to %" PRIu8 "\n", rx_gain);
                                    }
                                }
                                context_changed = true;
                            }
                            // Trigger context update (next)
                            last_context = start_time - std::chrono::milliseconds(2 * VRT_CONTEXT_INTERVAL);
                        }
                    }
                }
            }
        }

        if (bw_summary and now - last_update > std::chrono::seconds(1)) {
            const double dt = std::chrono::duration<double>(now - last_update).count();
            const double n  = 2.0 * (double)std::max<uint64_t>(progress_samps, 1);
            std::cout << boost::format("\t%.6f Msps, peak %5.1f dBFS, rms %5.1f dBFS, "
                                       "%.2f%% clip, over-range blocks %u, gaps %u")
                             % ((double)progress_samps / dt / 1e6)
                             % (20 * log10(std::max<double>(progress_peak, 1) / SIDEKIQ_FULL_SCALE))
                             % (10 * log10(std::max(progress_sum_sq / n, 1e-3)
                                           / (SIDEKIQ_FULL_SCALE * SIDEKIQ_FULL_SCALE)))
                             % (100.0 * (double)progress_clipped / n) % num_overload_blocks
                             % num_gaps;
            if (use_pps) {
                std::cout << (pps_present ? ", PPS ok" : ", PPS lost");
                if (clock_ppm_valid)
                    std::cout << boost::format(", clock %+.3f ppm") % clock_ppm;
                if (have_timestamp_adjustment)
                    std::cout << boost::format(", ts adj %+.0f ns") % ((double)timestamp_adjustment_ps / 1e3);
            }
            std::cout << std::endl;
            last_update      = now;
            progress_samps   = 0;
            progress_peak    = 0;
            progress_sum_sq  = 0;
            progress_clipped = 0;
        }
    }

    const auto actual_stop_time = std::chrono::steady_clock::now();

    // clean up transmit worker
    stop_signal_called = true;

    if (enable_tx and transmit_thread.joinable()) {
        transmit_thread.join();
    }

    if (not critical_error)
        skiq_stop_rx_streaming(card, curr_rx_hdl);
    release_card();

    if (stats) {
        std::cout << std::endl;
        const double actual_duration_seconds =
            std::chrono::duration<double>(actual_stop_time - start_time).count();

        std::cout << boost::format("Received %d samples in %f seconds") % num_total_samps
                         % actual_duration_seconds
                  << std::endl;
        const double avg_rate = (double)num_total_samps / actual_duration_seconds;
        std::cout << (avg_rate / 1e6) << " Msps" << std::endl;
        std::cout << boost::format("Discontinuities: %u (%u samples lost), overruns: %u, "
                                   "over-range blocks: %u")
                         % num_gaps % num_lost_samps % num_overruns % num_overload_blocks
                  << std::endl;
        if (use_pps) {
            std::cout << boost::format("PPS: %u edges, %u ignored, %u with stale RF timestamp, "
                                       "lost %u times, max |timestamp adjustment| %.0f ns")
                             % pps_edges % pps_glitches % pps_stale_rf % pps_losses
                             % ((double)max_abs_adjustment_ps / 1e3);
            if (clock_ppm_valid)
                std::cout << boost::format(", clock %+.3f ppm") % clock_ppm;
            std::cout << std::endl;
        }
    }

    // wait for ZMQ a bit
    std::this_thread::sleep_for(std::chrono::milliseconds(int64_t(1000 * setup_time)));

    // don't block on exit for messages a slow subscriber has not taken yet
    int linger = 0;
    for (void* s : {zmq_server, zmq_control, zmq_transmit, merge_zmq}) {
        if (s) {
            zmq_setsockopt(s, ZMQ_LINGER, &linger, sizeof linger);
            zmq_close(s);
        }
    }
    zmq_ctx_destroy(context);

    // finished
    std::cout << std::endl << "Done!" << std::endl << std::endl;

    return critical_error ? EXIT_FAILURE : EXIT_SUCCESS;
}
