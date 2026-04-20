//
// Copyright 2026 by Tammo Jan Dijkema
//
// SPDX-License-Identifier: MIT
//

#define _FILE_OFFSET_BITS 64

#include <hackrf.h>

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <math.h>
#include <chrono>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <inttypes.h>

#include <zmq.h>
#include <assert.h>

// VRT
#include <iostream>
#include <vrt/vrt_init.h>
#include <vrt/vrt_string.h>
#include <vrt/vrt_types.h>
#include <vrt/vrt_util.h>
#include <vrt/vrt_write.h>
#include <vrt/vrt_read.h>

#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <boost/circular_buffer.hpp>
#include <boost/thread/thread.hpp>

#include "vrt-tools.h"

#if defined(__GNUC__)
	#include <unistd.h>
	#include <sys/time.h>
#endif

#include <signal.h>

#define FD_BUFFER_SIZE (8 * 1024)

#define FREQ_ONE_MHZ (1000000.0)

#define DEFAULT_FREQ_HZ (900000000.0)  /* 900MHz */
#define FREQ_ABS_MIN_HZ (0.0)          /* 0 Hz */
#define FREQ_MIN_HZ     (1000000.0)    /* 1MHz */
#define FREQ_MAX_HZ     (6000000000.0) /* 6000MHz */
#define FREQ_ABS_MAX_HZ (7250000000.0) /* 7250MHz */
#define IF_ABS_MIN_HZ   (2000000000.0)
#define IF_MIN_HZ       (2170000000.0)
#define IF_MAX_HZ       (2740000000.0)
#define IF_ABS_MAX_HZ   (3000000000.0)
#define LO_MIN_HZ       (84375000.0)
#define LO_MAX_HZ       (5400000000.0)
#define DEFAULT_LO_HZ   (1000000000.0)

#define SAMPLE_RATE_MIN_HZ     (2000000.0)  /* 2MHz min sample rate */
#define SAMPLE_RATE_MAX_HZ     (20000000.0) /* 20MHz max sample rate */
#define DEFAULT_SAMPLE_RATE_HZ (10000000.0) /* 10MHz default sample rate */

#define DEFAULT_BASEBAND_FILTER_BANDWIDTH (5000000) /* 5MHz default */

#define SAMPLES_TO_XFER_MAX (0x8000000000000000ull) /* Max value */

#define BASEBAND_FILTER_BW_MIN (1750000)  /* 1.75 MHz min value */
#define BASEBAND_FILTER_BW_MAX (28000000) /* 28 MHz max value */

namespace po = boost::program_options;

boost::circular_buffer<int16_t> cb(65536*2*4*2);
boost::shared_mutex _access;

typedef enum {
	TRANSCEIVER_MODE_OFF = 0,
	TRANSCEIVER_MODE_RX = 1,
	TRANSCEIVER_MODE_TX = 2,
	TRANSCEIVER_MODE_SS = 3,
} transceiver_mode_t;

typedef enum {
	HW_SYNC_MODE_OFF = 0,
	HW_SYNC_MODE_ON = 1,
} hw_sync_mode_t;

typedef struct {
	uint64_t m0_total;
	uint64_t m4_total;
} stats_t;

static transceiver_mode_t transceiver_mode = TRANSCEIVER_MODE_RX;

static float TimevalDiff(const struct timeval* a, const struct timeval* b)
{
	return (a->tv_sec - b->tv_sec) + 1e-6f * (a->tv_usec - b->tv_usec);
}

static volatile bool do_exit = false;
static volatile bool interrupted = false;
static volatile bool flush_complete = false;

volatile uint32_t byte_count = 0;

/* sum of power of all samples, reset on the periodic report */
volatile uint64_t stream_power = 0;

struct timeval time_start;
struct timeval t_start;

bool automatic_tuning = false;
double freq_hz;

bool if_freq = false;
double if_freq_hz;

bool lo_freq = false;
double lo_freq_hz = DEFAULT_LO_HZ;

bool image_reject = false;
uint32_t image_reject_selection;

bool amp = false;

bool bias_tee = false;

bool sample_rate = false;
double sample_rate_hz;

bool force_ranges = false;

bool limit_num_samples = false;
uint64_t samples_to_xfer = 0;
size_t bytes_to_xfer = 0;

bool display_stats = false;

bool baseband_filter_bw = false;
double baseband_filter_bw_hz = 0;

bool crystal_correct = false;
uint32_t crystal_correct_ppm;

void stop_main_loop(void)
{
	do_exit = true;
}

int rx_callback(hackrf_transfer* transfer)
{
	size_t bytes_to_write;
	unsigned int i;

	/* Accumulate power (magnitude squared). */
	bytes_to_write = transfer->valid_length;
	uint64_t sum = 0;
	for (i = 0; i < bytes_to_write; i++) {
		int8_t value = transfer->buffer[i];
		sum += value * value;
	}

	/* Update both running totals at approximately the same time. */
	byte_count += transfer->valid_length;
	stream_power += sum;

	if (limit_num_samples) {
		if (bytes_to_write >= bytes_to_xfer) {
			bytes_to_write = bytes_to_xfer;
		}
		bytes_to_xfer -= bytes_to_write;
		if (bytes_to_xfer == 0) {
			stop_main_loop();
		}
	}

	/* Push samples (as int16_t) into circular buffer. */
	{
		boost::unique_lock<boost::shared_mutex> lock(_access);
		int8_t* buf = (int8_t*)transfer->buffer;
		for (size_t j = 0; j < bytes_to_write; j++) {
			cb.push_back((int16_t)buf[j]);
		}
	}

	return 0;
}


static int update_stats(hackrf_device* device, hackrf_m0_state* state, stats_t* stats)
{
	int result = hackrf_get_m0_state(device, state);

	if (result == HACKRF_SUCCESS) {
		if (state->m0_count < (stats->m0_total & 0xFFFFFFFF))
			stats->m0_total += 0x100000000;
		if (state->m4_count < (stats->m4_total & 0xFFFFFFFF))
			stats->m4_total += 0x100000000;
		stats->m0_total =
			(stats->m0_total & 0xFFFFFFFF00000000) | state->m0_count;
		stats->m4_total =
			(stats->m4_total & 0xFFFFFFFF00000000) | state->m4_count;
	}

	return result;
}

static hackrf_device* device = NULL;

void sigint_callback_handler(int signum)
{
	interrupted = true;
	fprintf(stderr, "Caught signal %d\n", signum);
	do_exit = true;
}


#define PATH_FILE_MAX_LEN (FILENAME_MAX)
#define DATE_TIME_MAX_LEN (32)

int main(int argc, char** argv)
{
	int opt;
	char path_file[PATH_FILE_MAX_LEN];
	char date_time[DATE_TIME_MAX_LEN];
	std::string path;
	std::string serial_number;
	char* endptr = NULL;
	int result;
	time_t rawtime;
	struct tm* timeinfo;
	long int file_pos;
	int exit_code = EXIT_SUCCESS;
	struct timeval t_end;
	float time_diff;
	unsigned int lna_gain = 8, vga_gain = 20;
	hackrf_m0_state state;
	stats_t stats = {0, 0};

	uint16_t instance, port;
	int hwm;

	po::options_description desc("hackrf_to_vrt options");
	desc.add_options()
	("help,h", "help message")
	("hw-sync", "enable hardware sync")
	("serial", po::value<std::string>(&serial_number), "serial number of desired HackRF")
	("freq", po::value<double>(&freq_hz),
		(boost::format(
		"RF center frequency in Hz [%.0fMHz to %.0fMHz supported, %.0fMHz to %.0fMHz forceable]")
		% (FREQ_MIN_HZ    / FREQ_ONE_MHZ)
		% (FREQ_MAX_HZ    / FREQ_ONE_MHZ)
		% (FREQ_ABS_MIN_HZ / FREQ_ONE_MHZ)
		% (FREQ_ABS_MAX_HZ / FREQ_ONE_MHZ)
		).str().c_str())
	("if-freq", po::value<double>(&if_freq_hz),
		(boost::format(
		"Intermediate frequency in Hz [%.0fMHz to %.0fMHz supported, %.0fMHz to %.0fMHz forceable]")
		% (IF_MIN_HZ / FREQ_ONE_MHZ)
		% (IF_MAX_HZ / FREQ_ONE_MHZ)
		% (IF_ABS_MIN_HZ / FREQ_ONE_MHZ)
		% (IF_ABS_MAX_HZ / FREQ_ONE_MHZ)
		).str().c_str())
	("lo-freq", po::value<double>(&lo_freq_hz),
		(boost::format(
		"LO frequency in Hz [%.0fMHz to %.0fMHz]")
		% (LO_MIN_HZ / FREQ_ONE_MHZ)
		% (LO_MAX_HZ / FREQ_ONE_MHZ)
		).str().c_str())
	("image-reject", po::value<uint32_t>(&image_reject_selection), "image rejection filter selection (0=auto, 1=low, 2=high)")
	("amp", "enable RF amplifier enable")
	("bias-tee", "enable bias tee power")
	("lna-gain,l", po::value<uint32_t>(&lna_gain), "LNA gain [0-40dB, 8dB steps]")
	("vga-gain,g", po::value<uint32_t>(&vga_gain), "VGA (baseband) RX gain [0-62dB, 2dB steps]")
	("rate", po::value<double>(&sample_rate_hz),
		(boost::format(
		"sample rate in Hz (%.0f-%.0fMHz supported, default %.0fMHz)")
		% (SAMPLE_RATE_MIN_HZ / FREQ_ONE_MHZ)
		% (SAMPLE_RATE_MAX_HZ / FREQ_ONE_MHZ)
		% (DEFAULT_SAMPLE_RATE_HZ / FREQ_ONE_MHZ)
		).str().c_str())
	("force-ranges,F", "force setting ranges")
	("nsamps", po::value<uint64_t>(&samples_to_xfer), "number of samples to transfer")
	("progress", "periodically display short-term bandwidth (not implemented)")
	("stats", "display stats")
	("bw", po::value<double>(&baseband_filter_bw_hz), "baseband filter bandwidth in Hz")
	("crystal-correct", po::value<uint32_t>(&crystal_correct_ppm), "crystal correction in ppm")
	("port", po::value<uint16_t>(&port), "VRT ZMQ port")
        ("hwm", po::value<int>(&hwm)->default_value(10000), "VRT ZMQ HWM")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("help")) {
	    std::cout << "HackRF RX samples to VRT" << std::endl << desc << std::endl;
	    return EXIT_SUCCESS;
	}

	// Boolean flags
	bool hw_sync        = vm.count("hw-sync") > 0;
	bool force_ranges   = vm.count("force-ranges") > 0;
	bool bw_summary     = vm.count("progress") > 0;
	bool display_stats  = vm.count("stats") > 0;

	// Mode flags
	bool automatic_tuning = vm.count("freq") > 0;
	bool if_freq        = vm.count("if-freq") > 0;
	bool lo_freq        = vm.count("lo-freq") > 0;
	bool image_reject   = vm.count("image-reject") > 0;
	bool amp            = vm.count("amp") > 0;
	bool bias_tee       = vm.count("bias-tee") > 0;
	bool sample_rate    = vm.count("rate") > 0;
	bool limit_num_samples = vm.count("nsamps") > 0;
	bool baseband_filter_bw = vm.count("bw") > 0;
	bool crystal_correct = vm.count("crystal-correct") > 0;

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

	void *context = zmq_ctx_new();
	void *responder = zmq_socket(context, ZMQ_PUB);
	int rc = zmq_setsockopt (responder, ZMQ_SNDHWM, &hwm, sizeof hwm);
	assert(rc == 0);

	std::string connect_string = "tcp://*:" + std::to_string(main_port);
	rc = zmq_bind(responder, connect_string.c_str());
	assert (rc == 0);
	zmq_server = responder;

	// Continue HackRF
	if (limit_num_samples) {
	    bytes_to_xfer = samples_to_xfer * 2ull;
	}

	if (lna_gain % 8)
		fprintf(stderr, "warning: lna_gain (-l) must be a multiple of 8\n");

	if (vga_gain % 2)
		fprintf(stderr, "warning: vga_gain (-g) must be a multiple of 2\n");

	if (samples_to_xfer >= SAMPLES_TO_XFER_MAX) {
		std::cerr << boost::format(
			"argument error: num_samples must be less than %" PRIu64 "/%" PRIu64 "Mio\n")
			% SAMPLES_TO_XFER_MAX
			% (SAMPLES_TO_XFER_MAX / (uint64_t)FREQ_ONE_MHZ);
		return EXIT_FAILURE;
	}

	if (if_freq || lo_freq || image_reject) {
		/* explicit tuning selected */
		if (!if_freq) {
			fprintf(stderr,
				"argument error: if_freq_hz must be specified for explicit tuning.\n");
			return EXIT_FAILURE;
		}
		if (!image_reject) {
			fprintf(stderr,
				"argument error: image_reject must be specified for explicit tuning.\n");
			return EXIT_FAILURE;
		}
		if (!lo_freq && (image_reject_selection != RF_PATH_FILTER_BYPASS)) {
			fprintf(stderr,
				"argument error: lo_freq_hz must be specified for explicit tuning unless image_reject is set to bypass.\n");
			return EXIT_FAILURE;
		}
		if (((if_freq_hz > IF_MAX_HZ) || (if_freq_hz < IF_MIN_HZ)) &&
		    !force_ranges) {
			std::cerr << boost::format(
				"argument error: if_freq_hz should be between %.0f and %.0f.\n")
				% IF_MIN_HZ % IF_MAX_HZ;
			return EXIT_FAILURE;
		}
		if ((if_freq_hz > IF_ABS_MAX_HZ) || (if_freq_hz < IF_ABS_MIN_HZ)) {
			std::cerr << boost::format(
				"argument error: if_freq_hz must be between %.0f and %.0f.\n")
				% IF_ABS_MIN_HZ % IF_ABS_MAX_HZ;
			return EXIT_FAILURE;
		}
		if ((lo_freq_hz > LO_MAX_HZ) || (lo_freq_hz < LO_MIN_HZ)) {
			std::cerr << boost::format(
				"argument error: lo_freq_hz shall be between %.0f and %.0f.\n")
				% LO_MIN_HZ % LO_MAX_HZ;
			return EXIT_FAILURE;
		}
		if (image_reject_selection > 2) {
			fprintf(stderr,
				"argument error: image_reject must be 0, 1, or 2 .\n");
			return EXIT_FAILURE;
		}
		if (automatic_tuning) {
			fprintf(stderr,
				"warning: freq_hz ignored by explicit tuning selection.\n");
			automatic_tuning = false;
		}
		switch (image_reject_selection) {
		case RF_PATH_FILTER_BYPASS:
			freq_hz = if_freq_hz;
			break;
		case RF_PATH_FILTER_LOW_PASS:
			freq_hz = std::fabs(if_freq_hz - lo_freq_hz);
			break;
		case RF_PATH_FILTER_HIGH_PASS:
			freq_hz = if_freq_hz + lo_freq_hz;
			break;
		default:
			freq_hz = DEFAULT_FREQ_HZ;
			break;
		}
		std::cerr << boost::format("explicit tuning specified for %.0f Hz.\n") % freq_hz;

	} else if (automatic_tuning) {
		if (((freq_hz > FREQ_MAX_HZ) || (freq_hz < FREQ_MIN_HZ)) &&
		    !force_ranges) {
			std::cerr << boost::format(
				"argument error: freq_hz should be between %.0f and %.0f.\n")
				% FREQ_MIN_HZ % FREQ_MAX_HZ;
			return EXIT_FAILURE;
		}
		if (freq_hz > FREQ_ABS_MAX_HZ) {
			std::cerr << boost::format(
				"argument error: freq_hz must be between %.0f and %.0f.\n")
				% FREQ_ABS_MIN_HZ % FREQ_ABS_MAX_HZ;
			return EXIT_FAILURE;
		}
	} else {
		/* Use default freq */
		freq_hz = DEFAULT_FREQ_HZ;
		automatic_tuning = true;
	}

	if (sample_rate) {
		if (sample_rate_hz > SAMPLE_RATE_MAX_HZ && !force_ranges) {
			std::cerr << boost::format(
				"argument error: sample_rate_hz should be <= %.0f Hz / %.03f MHz\n")
				% SAMPLE_RATE_MAX_HZ % (SAMPLE_RATE_MAX_HZ / FREQ_ONE_MHZ);
			return EXIT_FAILURE;
		}
		if (sample_rate_hz < SAMPLE_RATE_MIN_HZ && !force_ranges) {
			std::cerr << boost::format(
				"argument error: sample_rate_hz should be >= %.0f Hz / %.03f MHz\n")
				% SAMPLE_RATE_MIN_HZ % (SAMPLE_RATE_MIN_HZ / FREQ_ONE_MHZ);
			return EXIT_FAILURE;
		}
	} else {
		std::cerr << "Setting default sample rate " << std::endl;
		sample_rate_hz = DEFAULT_SAMPLE_RATE_HZ;
	}

	if (baseband_filter_bw) {
		if (baseband_filter_bw_hz > BASEBAND_FILTER_BW_MAX) {
			std::cerr << boost::format(
				"argument error: baseband_filter_bw_hz must be <= %.0f Hz / %.03f MHz\n")
				% BASEBAND_FILTER_BW_MAX % (BASEBAND_FILTER_BW_MAX / FREQ_ONE_MHZ);
			return EXIT_FAILURE;
		}

		if (baseband_filter_bw_hz < BASEBAND_FILTER_BW_MIN) {
			std::cerr << boost::format(
				"argument error: baseband_filter_bw_hz must be >= %.0f Hz / %.03f MHz\n")
				% BASEBAND_FILTER_BW_MIN % (BASEBAND_FILTER_BW_MIN / FREQ_ONE_MHZ);
			return EXIT_FAILURE;
		}

		/* Compute nearest freq for bw filter */
		baseband_filter_bw_hz =
			hackrf_compute_baseband_filter_bw((uint32_t)baseband_filter_bw_hz);
	}

	transceiver_mode = TRANSCEIVER_MODE_RX;

	// Change the freq and sample rate to correct the crystal clock error.
	if (crystal_correct) {
		double ppm_factor = (1000000.0 - crystal_correct_ppm) / 1000000.0;
		sample_rate_hz = sample_rate_hz * ppm_factor;
		freq_hz        = freq_hz        * ppm_factor;
	}

	result = hackrf_init();
	if (result != HACKRF_SUCCESS) {
		fprintf(stderr,
			"hackrf_init() failed: %s (%d)\n",
			hackrf_error_name((hackrf_error)result),
			result);
		return EXIT_FAILURE;
	}

	result = hackrf_open_by_serial(serial_number.c_str(), &device);
	if (result != HACKRF_SUCCESS) {
		fprintf(stderr,
			"hackrf_open() failed: %s (%d)\n",
			hackrf_error_name((hackrf_error)result),
			result);
		return EXIT_FAILURE;
	}

	signal(SIGINT, &sigint_callback_handler);
	signal(SIGILL, &sigint_callback_handler);
	signal(SIGFPE, &sigint_callback_handler);
	signal(SIGSEGV, &sigint_callback_handler);
	signal(SIGTERM, &sigint_callback_handler);
	signal(SIGABRT, &sigint_callback_handler);
	signal(SIGHUP, &sigint_callback_handler);

	std::cerr << boost::format("call hackrf_set_sample_rate(%.0f Hz / %.03f MHz)\n")
		% sample_rate_hz % (sample_rate_hz / FREQ_ONE_MHZ);
	result = hackrf_set_sample_rate(device, (uint32_t)sample_rate_hz);
	if (result != HACKRF_SUCCESS) {
		fprintf(stderr,
			"hackrf_set_sample_rate() failed: %s (%d)\n",
			hackrf_error_name((hackrf_error)result),
			result);
		return EXIT_FAILURE;
	}

	if (baseband_filter_bw) {
		std::cerr << boost::format("call hackrf_set_baseband_filter_bandwidth(%.0f Hz / %.03f MHz)\n")
			% baseband_filter_bw_hz % (baseband_filter_bw_hz / FREQ_ONE_MHZ);
		result = hackrf_set_baseband_filter_bandwidth(device, (uint32_t)baseband_filter_bw_hz);
		if (result != HACKRF_SUCCESS) {
			fprintf(stderr,
				"hackrf_set_baseband_filter_bandwidth() failed: %s (%d)\n",
				hackrf_error_name((hackrf_error)result),
				result);
			return EXIT_FAILURE;
		}
	}

	fprintf(stderr, "call hackrf_set_hw_sync_mode(%d)\n", hw_sync ? 1 : 0);
	result = hackrf_set_hw_sync_mode(
		device,
		hw_sync ? HW_SYNC_MODE_ON : HW_SYNC_MODE_OFF);
	if (result != HACKRF_SUCCESS) {
		fprintf(stderr,
			"hackrf_set_hw_sync_mode() failed: %s (%d)\n",
			hackrf_error_name((hackrf_error)result),
			result);
		return EXIT_FAILURE;
	}

	if (result != HACKRF_SUCCESS) {
		fprintf(stderr,
			"hackrf_start_?x() failed: %s (%d)\n",
			hackrf_error_name((hackrf_error)result),
			result);
		return EXIT_FAILURE;
	}

	if (automatic_tuning) {
		std::cerr << boost::format("call hackrf_set_freq(%.0f Hz / %.03f MHz)\n")
			% freq_hz % (freq_hz / FREQ_ONE_MHZ);
		result = hackrf_set_freq(device, (uint64_t)freq_hz);
		if (result != HACKRF_SUCCESS) {
			fprintf(stderr,
				"hackrf_set_freq() failed: %s (%d)\n",
				hackrf_error_name((hackrf_error)result),
				result);
			return EXIT_FAILURE;
		}
	} else {
		std::cerr << boost::format("call hackrf_set_freq_explicit() with %.0f Hz IF, %.0f Hz LO, %s\n")
			% if_freq_hz % lo_freq_hz
			% hackrf_filter_path_name((rf_path_filter)image_reject_selection);
		result = hackrf_set_freq_explicit(
			device,
			(uint64_t)if_freq_hz,
			(uint64_t)lo_freq_hz,
			(rf_path_filter)image_reject_selection);
		if (result != HACKRF_SUCCESS) {
			fprintf(stderr,
				"hackrf_set_freq_explicit() failed: %s (%d)\n",
				hackrf_error_name((hackrf_error)result),
				result);
			return EXIT_FAILURE;
		}
	}

	if (amp) {
		fprintf(stderr, "call hackrf_set_amp_enable(%u)\n", 1);
		result = hackrf_set_amp_enable(device, 1);
		if (result != HACKRF_SUCCESS) {
			fprintf(stderr,
				"hackrf_set_amp_enable() failed: %s (%d)\n",
				hackrf_error_name((hackrf_error)result),
				result);
			return EXIT_FAILURE;
		}
	}

	if (bias_tee) {
		fprintf(stderr, "call hackrf_set_antenna_enable(%u)\n", 1);
		result = hackrf_set_antenna_enable(device, 1);
		if (result != HACKRF_SUCCESS) {
			fprintf(stderr,
				"hackrf_set_antenna_enable() failed: %s (%d)\n",
				hackrf_error_name((hackrf_error)result),
				result);
			return EXIT_FAILURE;
		}
	}

	result = hackrf_set_vga_gain(device, vga_gain);
	result |= hackrf_set_lna_gain(device, lna_gain);
	result |= hackrf_start_rx(device, rx_callback, NULL);

	if (limit_num_samples) {
		std::cerr << boost::format("samples_to_xfer %" PRIu64 " / %" PRIu64 " Mio\n")
			% samples_to_xfer
			% (samples_to_xfer / (uint64_t)FREQ_ONE_MHZ);
	}

	gettimeofday(&t_start, NULL);
	gettimeofday(&time_start, NULL);

	fprintf(stderr, "Stop with Ctrl-C\n");

	size_t samps_per_buff = VRT_SAMPLES_PER_PACKET;

	uint32_t buffer[SIZE];
	int16_t bodydata[samps_per_buff * 2];

	uint32_t frame_count = 0;
	bool first_frame = true;
	bool context_changed = true;

	/* Warn if not standards compliant */
	if (vrt_is_platform_little_endian()) {
		printf("Warning: little endian support is work in progress.\n");
	}

	auto start_time = std::chrono::steady_clock::now();
	auto last_context = start_time;

	while (!do_exit) {
		struct timeval time_now;

		{
			boost::shared_lock<boost::shared_mutex> lock(_access);
			while (cb.size() > 2 * samps_per_buff) {

				gettimeofday(&time_now, NULL);

				if (first_frame) {
					fprintf(stderr, "First frame: %lu full secs, %.09f frac secs\n",
						time_now.tv_sec, time_now.tv_usec / 1e6);
					first_frame = false;
				}

				for (uint32_t i = 0; i < 2 * samps_per_buff; i++) {
					bodydata[i] = cb.front();
					cb.pop_front();
				}

				p.body = bodydata;
				p.header.packet_count = (uint8_t)(frame_count % 16);
				p.fields.integer_seconds_timestamp = time_now.tv_sec;
				p.fields.fractional_seconds_timestamp = 1e6 * time_now.tv_usec;

				zmq_msg_t msg;
				zmq_msg_init_size(&msg, SIZE * 4);
				int32_t rv = vrt_write_packet(&p, zmq_msg_data(&msg), SIZE, true);
				zmq_msg_send(&msg, zmq_server, 0);
				zmq_msg_close(&msg);

				frame_count++;

				const auto now = std::chrono::steady_clock::now();
				const auto time_since_last_context = now - last_context;
				if (time_since_last_context > std::chrono::milliseconds(200)) {

					last_context = now;

					struct vrt_packet pc;
					vrt_init_packet(&pc);
					vrt_init_context_packet(&pc);

					gettimeofday(&time_now, NULL);
					pc.fields.integer_seconds_timestamp = time_now.tv_sec;
					pc.fields.fractional_seconds_timestamp = 1e3 * time_now.tv_usec;

					pc.fields.stream_id = p.fields.stream_id;

					pc.if_context.bandwidth             = sample_rate_hz;
					pc.if_context.sample_rate           = sample_rate_hz;
					pc.if_context.rf_reference_frequency = (double)freq_hz;
					pc.if_context.rf_reference_frequency_offset = 0;
					pc.if_context.if_reference_frequency = 0;
					pc.if_context.gain.stage1            = (float)vga_gain;
					pc.if_context.gain.stage2            = (float)lna_gain;

					pc.if_context.state_and_event_indicators.reference_lock = false;
					pc.if_context.state_and_event_indicators.calibrated_time = false;

					pc.if_context.context_field_change_indicator = context_changed;

					int32_t rv = vrt_write_packet(&pc, buffer, SIZE, true);
					if (rv < 0) {
						fprintf(stderr, "Failed to write context packet: %s\n", vrt_string_error(rv));
					} else {
						zmq_send(zmq_server, buffer, rv * 4, 0);
					}

					context_changed = false;
				}
			}
		}
		usleep(10);
	}

	result = hackrf_is_streaming(device);
	if (do_exit) {
		fprintf(stderr, "\nExiting...\n");
	} else {
		fprintf(stderr,
			"\nExiting... hackrf_is_streaming() result: %s (%d)\n",
			hackrf_error_name((hackrf_error)result),
			result);
	}

	gettimeofday(&t_end, NULL);
	time_diff = TimevalDiff(&t_end, &t_start);
	fprintf(stderr, "Total time: %5.5f s\n", time_diff);

	if (device != NULL) {
		result = hackrf_stop_rx(device);
		if (result != HACKRF_SUCCESS) {
			fprintf(stderr,
				"hackrf_stop_rx() failed: %s (%d)\n",
				hackrf_error_name((hackrf_error)result),
				result);
		} else {
			fprintf(stderr, "hackrf_stop_rx() done\n");
		}

		if (display_stats) {
			result = update_stats(device, &state, &stats);
			if (result != HACKRF_SUCCESS) {
				fprintf(stderr,
					"hackrf_get_m0_state() failed: %s (%d)\n",
					hackrf_error_name((hackrf_error)result),
					result);
			} else {
				fprintf(stderr,
					"Transfer statistics:\n"
					"%" PRIu64 " bytes transferred by M0\n"
					"%" PRIu64 " bytes transferred by M4\n"
					"%u %s, longest %u bytes\n",
					stats.m0_total,
					stats.m4_total,
					state.num_shortfalls,
					"overruns",
					state.longest_shortfall);
			}
		}

		result = hackrf_close(device);
		if (result != HACKRF_SUCCESS) {
			fprintf(stderr,
				"hackrf_close() failed: %s (%d)\n",
				hackrf_error_name((hackrf_error)result),
				result);
		} else {
			fprintf(stderr, "hackrf_close() done\n");
		}

		hackrf_exit();
		fprintf(stderr, "hackrf_exit() done\n");
	}

	fprintf(stderr, "exit\n");
	return exit_code;
}
