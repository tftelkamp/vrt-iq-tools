#ifndef VRT_COMMON_H
#define VRT_COMMON_H

#include <chrono>
#include <complex>
#include <cmath>
#include <csignal>
#include <cstdint>
#include <iostream>

#include <boost/format.hpp>

// Shared signal handler
static bool stop_signal_called = false;
inline void sig_int_handler(int) { stop_signal_called = true; }

// Sample absolute value helpers (used by progress reporting)
template <typename samp_type>
inline float get_abs_val(samp_type t) { return std::fabs(t); }

inline float get_abs_val(std::complex<int16_t> t) { return std::fabs(t.real()); }
inline float get_abs_val(std::complex<int8_t> t) { return std::fabs(t.real()); }

// Progress reporting for consumer tools that receive ci16_le VRT data
struct progress_state {
    std::chrono::time_point<std::chrono::steady_clock> last_update;
    unsigned long long last_update_samps = 0;
};

inline void report_progress(
    progress_state& state,
    std::chrono::time_point<std::chrono::steady_clock> now,
    const uint32_t* buffer,
    uint32_t num_rx_samps)
{
    state.last_update_samps += num_rx_samps;

    const auto time_since_last_update = now - state.last_update;
    if (time_since_last_update > std::chrono::seconds(1)) {
        const double time_since_last_update_s =
            std::chrono::duration<double>(time_since_last_update).count();
        const double rate = double(state.last_update_samps) / time_since_last_update_s;
        std::cout << "\t" << (rate / 1e6) << " Msps, ";

        state.last_update_samps = 0;
        state.last_update = now;

        float sum_i = 0;
        uint32_t clip_i = 0;
        double datatype_max = 32768.;

        for (uint32_t i = 0; i < num_rx_samps; i++) {
            auto sample_i = get_abs_val((std::complex<int16_t>)buffer[i]);
            sum_i += sample_i;
            if (sample_i > datatype_max * 0.99)
                clip_i++;
        }
        sum_i = sum_i / num_rx_samps;
        std::cout << boost::format("%.0f") % (100.0 * log2(sum_i) / log2(datatype_max)) << "% I (";
        std::cout << boost::format("%.0f") % ceil(log2(sum_i) + 1) << " of ";
        std::cout << (int)ceil(log2(datatype_max) + 1) << " bits), ";
        std::cout << boost::format("%.0f") % (100.0 * clip_i / num_rx_samps) << "% I clip, ";
        std::cout << std::endl;
    }
}

#endif
