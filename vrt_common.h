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

#endif
