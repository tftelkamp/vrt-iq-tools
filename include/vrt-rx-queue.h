/* Receive path shared by the SDR tools whose driver delivers samples in a
 * callback (airspy_to_vrt, airspyhf_to_vrt).
 *
 * Samples travel from the driver's callback thread to the main loop in whole
 * transfer blocks. Every block carries the index of its first sample on the
 * stream timeline. That index also counts the samples that were lost on the
 * way, so packet timestamps follow the sample clock across a gap.
 *
 * The wall clock is read once, for the first block. All later timestamps are
 * counted in samples from there, with vrt_time_add_samples().
 */

#ifndef _VRTRXQUEUE_H
#define _VRTRXQUEUE_H

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <deque>
#include <mutex>
#include <vector>

#include "vrt-tools.h"

struct vrt_rx_block {
    std::vector<int16_t> data; // interleaved I/Q
    uint64_t first_sample = 0;
};

struct vrt_rx_queue {
    double rate       = 0;
    bool int_second   = false;
    size_t max_blocks = 64;

    // Blocks are discarded until the main loop sets this.
    std::atomic<bool> enabled{false};

    // Used by the callback thread only.
    bool started         = false;
    uint64_t next_sample = 0;

    // Protected by mutex.
    std::mutex mutex;
    std::condition_variable cv;
    std::deque<vrt_rx_block> queue;
    std::vector<std::vector<int16_t>> pool;
    struct vrt_time_ps time_first_sample = {0, 0};

    // Samples lost before the driver called us, and samples we had to throw
    // away because the main loop did not keep up.
    std::atomic<uint64_t> dropped_driver{0};
    std::atomic<uint64_t> dropped_queue{0};
};

/* Time of the sample that lies `samples` before t. */
inline struct vrt_time_ps vrt_time_before_samples(struct vrt_time_ps t, uint64_t samples, double rate) {
    uint64_t back_ps = (uint64_t)std::llround((double)samples * (double)VRT_PS_PER_SECOND / rate);

    t.seconds -= (int64_t)(back_ps / VRT_PS_PER_SECOND);
    back_ps %= VRT_PS_PER_SECOND;
    if (t.frac_ps < back_ps) {
        t.seconds -= 1;
        t.frac_ps += VRT_PS_PER_SECOND;
    }
    t.frac_ps -= back_ps;
    return t;
}

/* Queue one block of samples. Call this from the driver's callback.
 *
 * arrival is the wall clock at the moment the callback was entered, count the
 * number of complex samples, and dropped_samples the number of samples the
 * driver reports as lost ahead of this block. */
inline void vrt_rx_push(vrt_rx_queue* rx, const struct vrt_time_ps& arrival, const int16_t* samples,
                        uint64_t count, uint64_t dropped_samples) {
    if (not rx->enabled.load() or samples == nullptr or count == 0)
        return;

    uint64_t skip = 0;
    struct vrt_time_ps time_first_sample = {0, 0};

    if (not rx->started) {
        // The block has just arrived, so its first sample was taken one block
        // length ago.
        time_first_sample = vrt_time_before_samples(arrival, count, rx->rate);

        if (rx->int_second) {
            // Start at the sample that falls on the next integer second.
            const struct vrt_time_ps boundary = {
                time_first_sample.seconds + (time_first_sample.frac_ps > 0 ? 1 : 0), 0};
            skip = (uint64_t)std::llround(
                (double)vrt_time_diff_ps(boundary, time_first_sample) * 1e-12 * rx->rate);
            if (skip >= count)
                return; // the second does not start in this block
            time_first_sample = boundary;
        }
    } else if (dropped_samples > 0) {
        // The driver dropped whole transfers ahead of this one.
        rx->next_sample += dropped_samples;
        rx->dropped_driver += dropped_samples;
    }

    const uint64_t num_samps = count - skip;
    const int16_t* first     = samples + 2 * skip;
    bool queued              = false;

    {
        std::lock_guard<std::mutex> lock(rx->mutex);

        if (not rx->started)
            rx->time_first_sample = time_first_sample;

        if (rx->queue.size() < rx->max_blocks) {
            vrt_rx_block block;
            if (not rx->pool.empty()) {
                block.data = std::move(rx->pool.back());
                rx->pool.pop_back();
            }
            block.data.assign(first, first + 2 * num_samps);
            block.first_sample = rx->next_sample;
            rx->queue.push_back(std::move(block));
            queued = true;
        }
    }

    if (queued)
        rx->cv.notify_one();
    else
        rx->dropped_queue += num_samps;

    rx->started = true;
    rx->next_sample += num_samps;
}

/* Take the next block off the queue, waiting up to timeout for one. Also
 * gives the time of sample 0 of the stream. Returns false if there is none. */
inline bool vrt_rx_pop(vrt_rx_queue* rx, vrt_rx_block* block, struct vrt_time_ps* time_first_sample,
                       std::chrono::milliseconds timeout) {
    std::unique_lock<std::mutex> lock(rx->mutex);
    rx->cv.wait_for(lock, timeout, [rx] { return not rx->queue.empty(); });
    if (rx->queue.empty())
        return false;

    *block = std::move(rx->queue.front());
    rx->queue.pop_front();
    *time_first_sample = rx->time_first_sample;
    return true;
}

/* Hand a block's buffer back for reuse. */
inline void vrt_rx_recycle(vrt_rx_queue* rx, vrt_rx_block* block) {
    std::lock_guard<std::mutex> lock(rx->mutex);
    rx->pool.push_back(std::move(block->data));
}

/* Cuts the blocks into packets of a fixed number of samples.
 *
 * A packet never spans a gap in the stream: what was collected before the gap
 * is dropped, and the first packet after it is flagged, so that the caller can
 * skip a packet count and consumers see a lost frame. */
struct vrt_rx_packetizer {
    size_t samps_per_packet;
    std::vector<int16_t> body; // interleaved I/Q of the packet being handed out

    size_t fill                  = 0; // complex samples collected for the next packet
    uint64_t packet_first_sample = 0; // timeline index of the first of them
    uint64_t expected_sample     = 0; // timeline index that continues the stream
    uint64_t dropped_partial     = 0; // samples of packets cut short by a gap
    bool gap_pending             = false;

    explicit vrt_rx_packetizer(size_t samps) : samps_per_packet(samps), body(2 * samps) {}

    /* Calls on_packet(first_sample, after_gap) for every packet that the block
     * completes; body holds the packet's samples during the call. */
    template <typename F> void feed(const vrt_rx_block& block, F&& on_packet) {
        const size_t block_samps = block.data.size() / 2;

        if (block.first_sample != expected_sample) {
            dropped_partial += fill;
            fill        = 0;
            gap_pending = true;
        }

        size_t pos = 0;

        while (pos < block_samps) {
            if (fill == 0)
                packet_first_sample = block.first_sample + pos;

            const size_t take = std::min(samps_per_packet - fill, block_samps - pos);
            std::memcpy(&body[2 * fill], &block.data[2 * pos], take * 2 * sizeof(int16_t));
            fill += take;
            pos += take;

            if (fill < samps_per_packet)
                break;

            fill = 0;
            const bool after_gap = gap_pending;
            gap_pending          = false;
            on_packet(packet_first_sample, after_gap);
        }

        expected_sample = block.first_sample + block_samps;
    }
};

/* Print one line of progress, in the format of usrp_to_vrt: the measured
 * sample rate, and the peak level and clipping of one packet of samples. */
inline void vrt_rx_print_progress(double rate, unsigned channel, const int16_t* iq, size_t num_samps) {
    const double datatype_max = 32767.;

    double max_iq    = 0;
    uint32_t clip_iq = 0;

    for (size_t i = 0; i < num_samps; i++) {
        const double sample_i = std::fabs((double)iq[2 * i]);
        const double sample_q = std::fabs((double)iq[2 * i + 1]);
        max_iq = std::fmax(max_iq, std::fmax(sample_i, sample_q));
        if (sample_i > datatype_max * 0.99 || sample_q > datatype_max * 0.99)
            clip_iq++;
    }

    printf("\t%.6f Msps, CH%u: %3.0f dBFS (%2.0f/%d bits), %2.0f%% clip. \n", rate / 1e6, channel,
           20 * std::log10(max_iq / datatype_max), std::ceil(std::log2(max_iq) + 1),
           (int)std::ceil(std::log2(datatype_max) + 1), 100.0 * clip_iq / num_samps);
    fflush(stdout);
}

#endif
