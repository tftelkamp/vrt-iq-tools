//
// Copyright 2025 by Thomas Telkamp 
//
// SPDX-License-Identifier: MIT
//

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include <chrono>
#include <cmath>
#include <complex>
#include <csignal>
#include <cstring>
#include <fstream>
#include <iostream>
#include <thread>
#include <vector>

#include <zmq.h>
#include <assert.h>

// VRT
#include <vrt/vrt_init.h>
#include <vrt/vrt_string.h>
#include <vrt/vrt_types.h>
#include <vrt/vrt_util.h>
#include <vrt/vrt_write.h>
#include <vrt/vrt_read.h>

#include <inttypes.h>

#include <libairspy/airspy.h>

// VRT tools functions
#include "vrt-tools.h"
#include "vrt-rx-queue.h"

#define DEFAULT_VGA_IF_GAIN (5)
#define DEFAULT_LNA_GAIN (1)
#define DEFAULT_MIXER_GAIN (5)

#define VGA_GAIN_MAX (15)
#define MIXER_GAIN_MAX (15)
#define LNA_GAIN_MAX (14)
#define LINEARITY_GAIN_MAX (21)
#define SENSITIVITY_GAIN_MAX (21)

#define FREQ_HZ_MIN (24000000.0)
#define FREQ_HZ_MAX (1900000000.0) // as airspy_rx; officially 1750 MHz
#define FREQ_HZ_SPEC_MAX (1750000000.0)

unsigned long long num_total_samps = 0;

// R820T gain per index step, in tenths of a dB. Measured by Steve Markgraf on
// an RTL-SDR dongle at 928 MHz (librtlsdr, tuner_r82xx.c), so the gain reported
// for an Airspy is approximate, and relative rather than absolute.
#define R820T_VGA_BASE_GAIN (-47)
static const int r820t_vga_gain_steps[] = {
    0, 26, 26, 30, 42, 35, 24, 13, 14, 32, 36, 34, 35, 37, 35, 36
};
static const int r820t_lna_gain_steps[] = {
    0, 9, 13, 40, 38, 13, 31, 22, 26, 31, 26, 14, 19, 5, 35, 13
};
static const int r820t_mixer_gain_steps[] = {
    0, 5, 10, 10, 19, 9, 10, 25, 17, 10, 8, 16, 13, 6, 3, -8
};

// The VGA, mixer and LNA index that libairspy sets for each linearity and
// sensitivity gain value (libairspy, airspy.c). Entry 0 is the highest gain:
// a gain value v selects entry GAIN_COUNT - 1 - v.
#define GAIN_COUNT (22)
static const uint8_t linearity_vga_gains[GAIN_COUNT] = { 13, 12, 11, 11, 11, 11, 11, 10, 10, 10, 10, 10, 10, 10, 10, 10, 9, 8, 7, 6, 5, 4 };
static const uint8_t linearity_mixer_gains[GAIN_COUNT] = { 12, 12, 11, 9, 8, 7, 6, 6, 5, 0, 0, 1, 0, 0, 2, 2, 1, 1, 1, 1, 0, 0 };
static const uint8_t linearity_lna_gains[GAIN_COUNT] = { 14, 14, 14, 13, 12, 10, 9, 9, 8, 9, 8, 6, 5, 3, 1, 0, 0, 0, 0, 0, 0, 0 };
static const uint8_t sensitivity_vga_gains[GAIN_COUNT] = { 13, 12, 11, 10, 9, 8, 7, 6, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4 };
static const uint8_t sensitivity_mixer_gains[GAIN_COUNT] = { 12, 12, 12, 12, 11, 10, 10, 9, 9, 8, 7, 4, 4, 4, 3, 2, 2, 1, 0, 0, 0, 0 };
static const uint8_t sensitivity_lna_gains[GAIN_COUNT] = { 14, 14, 14, 14, 14, 14, 14, 14, 14, 13, 12, 12, 9, 9, 8, 7, 6, 5, 3, 2, 1, 0 };

// Approximate gain in dB for a set of R820T gain indexes.
static double r820t_gain_db(uint32_t vga, uint32_t mixer, uint32_t lna)
{
    int tenth_db = R820T_VGA_BASE_GAIN;

    for (uint32_t i = 1; i <= vga and i < 16; i++)
        tenth_db += r820t_vga_gain_steps[i];
    for (uint32_t i = 1; i <= mixer and i < 16; i++)
        tenth_db += r820t_mixer_gain_steps[i];
    for (uint32_t i = 1; i <= lna and i < 16; i++)
        tenth_db += r820t_lna_gain_steps[i];

    return tenth_db / 10.0;
}

namespace po = boost::program_options;

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

// Frequency the R820T is really tuned to, for a requested centre frequency.
// Returns 0 if the tuner PLL cannot reach the frequency.
//
// After Phil Karn (ka9q). The PLL words are computed exactly as
// r820t_set_pll() in the Airspy firmware does, including its truncating
// divider: rounding to the nearest step instead gives a result that is one
// PLL step off for about 0.5% of all frequencies.
//
// The offset was measured against a signal generator on the same 10 MHz
// reference (firmware v1.0.0-rc10-6, 2.5 and 10 Msps, all six divider
// ranges). With a non-zero fraction the VCO runs a quarter of a step above the
// programmed value. With a zero fraction the firmware switches the sigma-delta
// modulator off, and there is no offset. Model and measurement agree to 0.01 Hz.
double true_freq(uint64_t freq_hz, uint32_t rate, int8_t fw_version) {
  const uint32_t vco_min = 1770000000u; // 1.77 GHz
  const uint32_t vco_max = 3900000000u; // 3.9 GHz

  // Tuner reference: 20 MHz for v1.0.0-rc10-0-g946184a 2016-09-19, 25 MHz for
  // v1.0.0-rc10-6-g4008185 2020-05-08, and assumed unchanged in newer firmware.
  const uint32_t xtal_freq = (fw_version == 0) ? 20000000u : 25000000u;

  const uint32_t if_freq = rate/2;
  const uint32_t lo_freq = (uint32_t)freq_hz + if_freq;

  // The firmware works in 32 bits throughout.
  uint32_t ref = xtal_freq >> 1;
  uint32_t div_num;
  uint32_t vco = 0;
  bool div_found = false;

  // Find divider to put VCO = f*2^(d+1) in range vco_min to vco_max
  //          MHz             step, Hz (25 MHz reference)
  // 0: 885.0     1950.0      190.735
  // 1: 442.50     885.00      95.367
  // 2: 221.25     442.50      47.684
  // 3: 110.625    221.25      23.842
  // 4:  55.3125   110.625     11.921
  // 5:  27.65625   55.312      5.960
  for (div_num = 0; div_num <= 5; div_num++) {
    vco = lo_freq << (div_num + 1);
    if (vco >= vco_min && vco <= vco_max) {
      div_found = true;
      break;
    }
  }
  if (!div_found)
    return 0; // Frequency out of range

  // PLL programming bits: Nint in bits 16 and up, Nfract (sdm) in the lower 16.
  // Long division of the VCO frequency by the step, with a rounding term.
  vco += ref >> 16;
  ref <<= 8;
  uint32_t mask = 1u << 23;
  uint32_t rem = 0;
  while (mask > 0 && vco > 0) {
    if (vco >= ref) {
      rem |= mask;
      vco -= ref;
    }
    ref >>= 1;
    mask >>= 1;
  }

  const uint32_t sdm = rem & 0xffff;
  const double offset = (sdm == 0) ? 0.0 : 0.25;

  // Freq steps are (xtal_freq/2) / 2^(16 + div_num) Hz
  return ((double)rem + offset) * (xtal_freq / 2.0) / (double)(1u << (div_num + 16)) - if_freq;
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
        return AIRSPY_SUCCESS;
    } else {
        return AIRSPY_ERROR_INVALID_PARAM;
    }
}

int rx_callback(airspy_transfer_t* transfer)
{
    const struct vrt_time_ps arrival = vrt_time_now();

    if (transfer->sample_count > 0)
        vrt_rx_push((vrt_rx_queue*)transfer->ctx, arrival, (const int16_t*)transfer->samples,
            (uint64_t)transfer->sample_count, transfer->dropped_samples);

    return 0;
}

int main(int argc, char* argv[])
{
    // variables to be set by po
    std::string merge_address, dev_given;
    size_t total_num_samps = 0;
    uint16_t instance, port, merge_port;
    uint32_t stream_id;
    int hwm;
    double gain = 0; // approximate, in dB
    double rate, freq, total_time = 0, setup_time, if_freq;
    bool merge;
    double real_freq;

    // Airspy
    struct airspy_device* device = 0;
    char airspy_version[255 + 1] = "";
    uint64_t serial_number_val;
    uint32_t *supported_samplerates;
    uint32_t sample_rate_u32 = 0;
    uint32_t nsrates;
    uint8_t register_value;
    int8_t fw_version = -1;

    uint32_t vga_gain = DEFAULT_VGA_IF_GAIN;
    uint32_t lna_gain = DEFAULT_LNA_GAIN;
    uint32_t mixer_gain = DEFAULT_MIXER_GAIN;

    uint32_t linearity_gain_val = 0;
    uint32_t sensitivity_gain_val = 0;
   
    // setup the program options
    po::options_description desc("Allowed options");
    // clang-format off
    desc.add_options()
        ("help", "help message")
        ("rate", po::value<double>(&rate)->default_value(2500000), "rate of incoming samples")
        ("freq", po::value<double>(&freq)->default_value(100e6), "RF center frequency in Hz")
        ("serial", po::value<std::string>(&dev_given), "device serial")
        ("if-freq", po::value<double>(&if_freq)->default_value(0.0), "IF center frequency in Hz")
        ("setup", po::value<double>(&setup_time)->default_value(1.0), "seconds of setup time")
        ("vga", po::value<uint32_t>(&vga_gain)->default_value(DEFAULT_VGA_IF_GAIN), "VGA gain for the RF chain")
        ("lna", po::value<uint32_t>(&lna_gain)->default_value(DEFAULT_LNA_GAIN), "LNA gain for the RF chain")
        ("mixer", po::value<uint32_t>(&mixer_gain)->default_value(DEFAULT_MIXER_GAIN), "MIXER gain for the RF chain")
        ("sensitivity", po::value<uint32_t>(&sensitivity_gain_val), "Sensitivity gain for the RF chain")
        ("linearity", po::value<uint32_t>(&linearity_gain_val), "Linearity gain for the RF chain")
        ("packing", "set packing for samples")
        ("progress", "periodically display short-term bandwidth")
        ("stats", "show average bandwidth on exit")
        ("int-second", "align start of reception to integer second")
        ("null", "run without streaming")
        ("continue", "don't abort on a bad packet")
        ("bias-tee", "Enable Bias Tee power")
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
        std::cout << boost::format("Airspy RX samples to VRT. %s") % desc << std::endl;
        std::cout << std::endl
                  << "This application streams data from a single channel of an Airspy "
                     "device to VRT.\n"
                  << std::endl;
        return ~0;
    }

    bool bw_summary             = vm.count("progress") > 0;
    bool stats                  = vm.count("stats") > 0;
    bool null                   = vm.count("null") > 0;
    bool continue_on_bad_packet = vm.count("continue") > 0;
    bool bias_tee               = vm.count("bias-tee") > 0;
    bool sensitivity_gain       = vm.count("sensitivity") > 0;
    bool linearity_gain         = vm.count("linearity") > 0;
    bool serial_number          = vm.count("serial") > 0;
    bool packing                = vm.count("packing") > 0;

    if (serial_number) {
        parse_u64(dev_given.c_str(), &serial_number_val);
    }

    if(vga_gain > VGA_GAIN_MAX) {
        fprintf(stderr, "argument error: vga_gain out of range\n");
        return EXIT_FAILURE;
    }

    if(mixer_gain > MIXER_GAIN_MAX) {
        fprintf(stderr, "argument error: mixer_gain out of range\n");
        return EXIT_FAILURE;
    }

    if(lna_gain > LNA_GAIN_MAX) {
        fprintf(stderr, "argument error: lna_gain out of range\n");
        return EXIT_FAILURE;
    }

    if(linearity_gain_val > LINEARITY_GAIN_MAX) {
        fprintf(stderr, "argument error: linearity_gain out of range\n");
        return EXIT_FAILURE;
    }

    if(sensitivity_gain_val > SENSITIVITY_GAIN_MAX) {
        fprintf(stderr, "argument error: sensitivity_gain out of range\n");
        return EXIT_FAILURE;
    }

    if( (linearity_gain == true) && (sensitivity_gain == true) )
    {
        fprintf(stderr, "argument error: linearity_gain and sensitivity_gain are both set (choose only one option)\n");
        return EXIT_FAILURE;
    }

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

    // create an airspy device

    int result = airspy_init();
    if( result != AIRSPY_SUCCESS ) {
        fprintf(stderr, "airspy_init() failed: %s (%d)\n", airspy_error_name((airspy_error)result), result);
        return EXIT_FAILURE;
    }

    if(serial_number == true)
    {
        result = airspy_open_sn(&device, serial_number_val);
        if( result != AIRSPY_SUCCESS ) {
            fprintf(stderr, "airspy_open_sn() failed: %s (%d)\n", airspy_error_name((airspy_error)result), result);
            airspy_exit();
            return EXIT_FAILURE;
        }
    } else {
        result = airspy_open(&device);
        if( result != AIRSPY_SUCCESS ) {
            fprintf(stderr, "airspy_open() failed: %s (%d)\n", airspy_error_name((airspy_error)result), result);
            airspy_exit();
            return EXIT_FAILURE;
        }
    }

    enum airspy_sample_type sample_type_val = AIRSPY_SAMPLE_INT16_IQ;

    result = airspy_set_sample_type(device, sample_type_val);
    if (result != AIRSPY_SUCCESS) {
        fprintf(stderr, "airspy_set_sample_type() failed: %s (%d)\n", airspy_error_name((airspy_error)result), result);
        airspy_close(device);
        airspy_exit();
        return EXIT_FAILURE;
    }

    // check sample rates
    airspy_get_samplerates(device, &nsrates, 0);
    supported_samplerates = (uint32_t *) malloc(nsrates * sizeof(uint32_t));
    airspy_get_samplerates(device, supported_samplerates, nsrates);

    int s;
    fprintf(stderr, "Supported sample rates:");
    for (s = 0; s < nsrates; s++) {
        fprintf(stderr, " %d kS/s", supported_samplerates[s]/1000);
        if (rate == supported_samplerates[s])
            sample_rate_u32 = supported_samplerates[s];
    }
    fprintf(stderr, "\n");

    free(supported_samplerates);

    if (sample_rate_u32==0) {
        fprintf(stderr, "unsupported sample rate: %f\n", rate);
        exit(1);
    }

    // set sample rate
    if (airspy_set_samplerate(device, sample_rate_u32) != AIRSPY_SUCCESS) {
        fprintf(stderr, "airspy_set_samplerate() failed: %d\n", sample_rate_u32);
        exit(1);
    }

    rate = sample_rate_u32;

    fprintf(stderr, "Sample rate: %f MS/s %s\n", sample_rate_u32 * 0.000001f, "IQ");

    airspy_read_partid_serialno_t read_partid_serialno;

    if (airspy_board_partid_serialno_read(device, &read_partid_serialno) != AIRSPY_SUCCESS) {
        fprintf(stderr, "airspy_board_partid_serialno_read() failed\n");
        exit(1);
    } else
        fprintf(stderr, "Device Serial Number: 0x%08X%08X\n",
                read_partid_serialno.serial_no[2],
                read_partid_serialno.serial_no[3]);

    result = airspy_version_string_read(device, &airspy_version[0], 255);
    if (result != AIRSPY_SUCCESS) {
        fprintf(stderr, "airspy_version_string_read() failed: %s (%d)\n",
            airspy_error_name((airspy_error)result), result);
    }

    if (strcmp(airspy_version,"AirSpy NOS v1.0.0-rc10-6-g4008185 2020-05-08") == 0)
        fw_version = 1;
    else if (strcmp(airspy_version,"AirSpy NOS v1.0.0-rc10-0-g946184a 2016-09-19") == 0)
        fw_version = 0;

    fprintf(stderr, "Firmware Version: %s (airspy_to_vrt: %i)\n", airspy_version, fw_version);

    result = airspy_set_packing(device, (uint8_t)0);
    if (packing) {
        result = airspy_set_packing(device, (uint8_t)packing);
        if( result != AIRSPY_SUCCESS ) {
            fprintf(stderr, "airspy_set_packing() failed: %s (%d)\n", airspy_error_name((airspy_error)result), result);
            airspy_close(device);
            airspy_exit();
            return EXIT_FAILURE;
        }
        fprintf(stderr, "Packing enabled\n");
    }

    result = airspy_set_rf_bias(device, (int)bias_tee);
    if( result != AIRSPY_SUCCESS ) {
        fprintf(stderr, "airspy_set_rf_bias() failed: %s (%d)\n", airspy_error_name((airspy_error)result), result);
        airspy_close(device);
        airspy_exit();
        return EXIT_FAILURE;
    }

    // check external ref
    bool ref = false;
    bool ref_signal = false;

    result = airspy_si5351c_read(device, 0, &register_value); 
    if (result == AIRSPY_SUCCESS)
        if ((register_value & 16) == 0)
            ref_signal = true;

    fprintf(stderr, "External reference signal: %s\n", ref_signal ? "on" : "off"); 

    result = airspy_si5351c_read(device, 15, &register_value); 
    if (result == AIRSPY_SUCCESS)
        if ((register_value & 0x0C) != 0)
            ref = true;

    fprintf(stderr, "External reference used: %s\n", ref ? "yes" : "no"); 

    // Frequency

    // The tuner takes whole Hz.
    freq = std::round(freq - if_freq);

    if (freq < FREQ_HZ_MIN or freq >= FREQ_HZ_MAX) {
        fprintf(stderr, "argument error: frequency %.0f Hz is outside the tuner range of %.0f to %.0f MHz\n",
            freq, FREQ_HZ_MIN/1e6, FREQ_HZ_MAX/1e6);
        airspy_close(device);
        airspy_exit();
        return EXIT_FAILURE;
    }
    if (freq > FREQ_HZ_SPEC_MAX)
        fprintf(stderr, "Warning: frequency %.0f Hz is above the specified maximum of %.0f MHz.\n",
            freq, FREQ_HZ_SPEC_MAX/1e6);

    // The firmware does not report it when the tuner PLL cannot reach a
    // frequency, it leaves the tuner where it was.
    const double pll_freq = true_freq(freq, sample_rate_u32, fw_version);
    if (pll_freq == 0) {
        fprintf(stderr, "argument error: frequency %.0f Hz is outside the range of the tuner PLL at this sample rate (minimum %.0f Hz)\n",
            freq, 27656250.0 - sample_rate_u32/2);
        airspy_close(device);
        airspy_exit();
        return EXIT_FAILURE;
    }

    if (ref) {
        if (fw_version < 0)
            fprintf(stderr, "Warning: unknown firmware version, the true frequency assumes a 25 MHz tuner reference.\n");
        real_freq = pll_freq;
        printf("True frequency (R820T): %f Hz\n",real_freq);
    } else
        real_freq = freq;

    if( (linearity_gain == false) && (sensitivity_gain == false) )
    {
        // Manual gain: switch the AGCs off first, as libairspy does for its
        // linearity and sensitivity settings. The tuner keeps whatever a
        // previous application left behind.
        result = airspy_set_mixer_agc(device, 0);
        if( result != AIRSPY_SUCCESS ) {
            fprintf(stderr, "airspy_set_mixer_agc() failed: %s (%d)\n", airspy_error_name((airspy_error)result), result);
        }

        result = airspy_set_lna_agc(device, 0);
        if( result != AIRSPY_SUCCESS ) {
            fprintf(stderr, "airspy_set_lna_agc() failed: %s (%d)\n", airspy_error_name((airspy_error)result), result);
        }

        result = airspy_set_vga_gain(device, vga_gain);
        if( result != AIRSPY_SUCCESS ) {
            fprintf(stderr, "airspy_set_vga_gain() failed: %s (%d)\n", airspy_error_name((airspy_error)result), result);
        }

        result = airspy_set_mixer_gain(device, mixer_gain);
        if( result != AIRSPY_SUCCESS ) {
            fprintf(stderr, "airspy_set_mixer_gain() failed: %s (%d)\n", airspy_error_name((airspy_error)result), result);
        }

        result = airspy_set_lna_gain(device, lna_gain);
        if( result != AIRSPY_SUCCESS ) {
            fprintf(stderr, "airspy_set_lna_gain() failed: %s (%d)\n", airspy_error_name((airspy_error)result), result);
        }
        gain = r820t_gain_db(vga_gain, mixer_gain, lna_gain);
    } else
    {
        if( linearity_gain == true )
        {
            result =  airspy_set_linearity_gain(device, linearity_gain_val);
            if( result != AIRSPY_SUCCESS ) {
                fprintf(stderr, "airspy_set_linearity_gain() failed: %s (%d)\n", airspy_error_name((airspy_error)result), result);
            }
            const uint32_t entry = GAIN_COUNT - 1 - linearity_gain_val;
            gain = r820t_gain_db(linearity_vga_gains[entry], linearity_mixer_gains[entry], linearity_lna_gains[entry]);
        }

        if( sensitivity_gain == true )
        {
            result =  airspy_set_sensitivity_gain(device, sensitivity_gain_val);
            if( result != AIRSPY_SUCCESS ) {
                fprintf(stderr, "airspy_set_sensitivity_gain() failed: %s (%d)\n", airspy_error_name((airspy_error)result), result);
            }
            const uint32_t entry = GAIN_COUNT - 1 - sensitivity_gain_val;
            gain = r820t_gain_db(sensitivity_vga_gains[entry], sensitivity_mixer_gains[entry], sensitivity_lna_gains[entry]);
        }
    }

    // Samples received during setup are discarded: the callback only queues
    // blocks once rx.enabled is set, just before the main loop.
    vrt_rx_queue rx;
    rx.rate = rate;

    fprintf(stderr, "Gain: %.1f dB (approximate)\n", gain);

    if( airspy_start_rx(device, rx_callback, &rx) != AIRSPY_SUCCESS ) {
        fprintf(stderr, "airspy_start() failed.\n");
        exit(1);
    }

    if( airspy_set_freq(device, (uint32_t)freq) != AIRSPY_SUCCESS ) {
        fprintf(stderr, "airspy_set_freq() failed.\n");
        exit(1);
    }

    // Sleep setup time
    std::this_thread::sleep_for(std::chrono::milliseconds(int64_t(1000 * setup_time)));

    struct vrt_time_ps time_now = vrt_time_now();

    // seed random generator with the current time
    srand((unsigned)(time_now.frac_ps/1000000 + time_now.seconds));

 	// Receive

    if (total_num_samps == 0) {
        std::signal(SIGINT, &sig_int_handler);
        std::cout << "Press Ctrl + C to stop streaming..." << std::endl;
    }

	size_t samps_per_buff = VRT_SAMPLES_PER_PACKET;

	unsigned long long num_requested_samples = total_num_samps;
    double time_requested = total_time;
    bool int_second             = (bool)vm.count("int-second");

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
    auto stop_time =
        start_time + std::chrono::milliseconds(int64_t(1000 * time_requested));
    
    // Track time and samps between updating the BW summary
    auto last_update                     = start_time;
    auto last_context                    = start_time;
    unsigned long long last_update_samps = 0;


    // Run this loop until either time expired (if a duration was given), until
    // the requested number of samples were collected (if such a number was
    // given), or until Ctrl-C was pressed.

    uint32_t frame_count = 0;

    uint32_t num_words_read=0;

    vrt_rx_packetizer packetizer(samps_per_buff);
    std::vector<int16_t>& bodydata = packetizer.body;

    struct vrt_time_ps time_first_sample = {0, 0};

    uint64_t dropped_reported = 0;
    auto last_drop_report = start_time;
    bool stream_failed = false;

    // flush merge queue
    if (merge)
        while ( zmq_recv(merge_zmq, merge_buffer.data(), merge_buffer.size(), ZMQ_NOBLOCK) > 0 ) { }

    // Start accepting samples. With --int-second the callback starts at the
    // sample that falls on the next integer second.
    rx.int_second = int_second;
    rx.enabled = true;

    while (not stop_signal_called) {

        const auto now = std::chrono::steady_clock::now();

        // libairspy stops for good on the first failed USB transfer.
        if (airspy_is_streaming(device) != AIRSPY_TRUE) {
            fprintf(stderr, "Error: Airspy stopped streaming (USB transfer failed).\n");
            stream_failed = true;
            break;
        }

        vrt_rx_block block;

        if (vrt_rx_pop(&rx, &block, &time_first_sample, std::chrono::milliseconds(1))) {

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

                num_words_read = samps_per_buff;

                num_total_samps += num_words_read;

                p.body = bodydata.data();
                p.header.packet_count = (uint8_t)frame_count%16;
                p.fields.integer_seconds_timestamp = time_now.seconds;
                p.fields.fractional_seconds_timestamp = time_now.frac_ps;

                zmq_msg_t msg;
                int rc = zmq_msg_init_size (&msg, SIZE*4);

                int32_t rv = vrt_write_packet(&p, zmq_msg_data(&msg), SIZE, true);

                frame_count++;

                // VRT
                zmq_msg_send(&msg, zmq_server, 0);
                zmq_msg_close(&msg);

                const auto time_since_last_context = now - last_context;
                if (time_since_last_context > std::chrono::milliseconds(200)) {

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
                    pc.if_context.rf_reference_frequency            = real_freq+if_freq;
                    pc.if_context.rf_reference_frequency_offset     = 0;
                    pc.if_context.if_reference_frequency            = if_freq; // 0 for Zero-IF
                    pc.if_context.gain.stage1                       = gain;
                    pc.if_context.gain.stage2                       = 0;

                    pc.if_context.state_and_event_indicators.reference_lock = ref;

                    pc.if_context.state_and_event_indicators.calibrated_time = false;

                    if (context_changed)
                        pc.if_context.context_field_change_indicator = true;
                    else
                        pc.if_context.context_field_change_indicator = false;

                    int32_t rv = vrt_write_packet(&pc, buffer, SIZE, true);
                    if (rv < 0) {
                        fprintf(stderr, "Failed to write packet: %s\n", vrt_string_error(rv));
                    }

                    // ZMQ
                    zmq_send (zmq_server, buffer, rv*4, 0);

                    context_changed = false;
                }

                if (bw_summary) {
                    last_update_samps += num_words_read;
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

            vrt_rx_recycle(&rx, &block);
        }

        // Report lost samples, at most once a second.
        const uint64_t dropped_total = rx.dropped_driver + rx.dropped_queue + packetizer.dropped_partial;
        if (dropped_total != dropped_reported and now - last_drop_report > std::chrono::seconds(1)) {
            fprintf(stderr, "Warning: lost %llu samples (%llu in total).\n",
                (unsigned long long)(dropped_total - dropped_reported),
                (unsigned long long)dropped_total);
            dropped_reported = dropped_total;
            last_drop_report = now;
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

        const uint64_t dropped_driver = rx.dropped_driver;
        const uint64_t dropped_queue  = rx.dropped_queue;
        std::cout << boost::format("Lost %d samples (%d in the driver, %d in the queue, %d in cut packets).")
                         % (dropped_driver + dropped_queue + packetizer.dropped_partial)
                         % dropped_driver % dropped_queue % packetizer.dropped_partial
                  << std::endl;
    }
  
    /* switch bias-tee off */
    airspy_set_rf_bias(device, 0);

    /* clean up */
    if( airspy_stop_rx(device) != AIRSPY_SUCCESS ) {
        fprintf(stderr, "airspy_stop() failed\n");
    }

    if( airspy_close(device) != AIRSPY_SUCCESS ) {
        fprintf(stderr,"airspy_close() failed\n");
    }

    // finished
    std::cout << std::endl << "Done!" << std::endl << std::endl;

    return stream_failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
