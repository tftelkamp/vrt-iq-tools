//
// Copyright 2021/2022 by Thomas Telkamp 
//
// SPDX-License-Identifier: GPL-3.0-or-later
//

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/circular_buffer.hpp>

#include <chrono>
#include <complex>
#include <csignal>
#include <fstream>
#include <iostream>
#include <thread>

#include <zmq.h>
#include <assert.h>

// VRT
#include <vrt/vrt_init.h>
#include <vrt/vrt_string.h>
#include <vrt/vrt_types.h>
#include <vrt/vrt_util.h>
#include <vrt/vrt_write.h>
#include <vrt/vrt_read.h>

#include <boost/thread/thread.hpp>

#include <inttypes.h>

#include "rtl-sdr.h"
#include "convenience.h"

// VRT tools functions
#include "vrt-tools.h"

#define DEFAULT_SAMPLE_RATE     1000000
// #define DEFAULT_BUF_LENGTH      (16 * 16384)
// #define DEFAULT_BUF_LENGTH      (8 * 16384)
#define DEFAULT_BUF_LENGTH      (20480)

#define MINIMAL_BUF_LENGTH      512
#define MAXIMAL_BUF_LENGTH      (256 * 16384)

static int do_exit = 0;
static uint32_t bytes_to_read = 0;
static rtlsdr_dev_t *dev = NULL;

int n_read;
uint8_t *rtl_buffer;
uint32_t out_block_size = DEFAULT_BUF_LENGTH;

unsigned long long num_total_samps = 0;

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

int main(int argc, char* argv[])
{
    // variables to be set by po
    std::string merge_address, dev_given;
    size_t total_num_samps = 0;
    uint16_t instance, port, merge_port;
    uint32_t stream_id;
    int hwm, gain;
    double rate, freq, total_time, setup_time, if_freq;
    bool merge;

    // setup the program options
    po::options_description desc("Allowed options");
    // clang-format off
    desc.add_options()
        ("help", "help message")
        ("rate", po::value<double>(&rate)->default_value(1e6), "rate of incoming samples")
        ("freq", po::value<double>(&freq)->default_value(0.0), "RF center frequency in Hz")
        ("device", po::value<std::string>(&dev_given), "device name or index")
        ("if-freq", po::value<double>(&if_freq)->default_value(0.0), "IF center frequency in Hz")
        ("gain", po::value<int>(&gain)->default_value(0), "gain for the RF chain (default AGC)")
        ("setup", po::value<double>(&setup_time)->default_value(1.0), "seconds of setup time")
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
        std::cout << boost::format("RTL-SDR RX samples to VRT. %s") % desc << std::endl;
        std::cout << std::endl
                  << "This application streams data from a single channel of an RTL-SDR "
                     "device to VRT.\n"
                  << std::endl;
        return ~0;
    }

    bool bw_summary             = vm.count("progress") > 0;
    bool stats                  = vm.count("stats") > 0;
    bool null                   = vm.count("null") > 0;
    bool continue_on_bad_packet = vm.count("continue") > 0;
    bool bias_tee               = vm.count("bias-tee") > 0;

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

    // Merge
    void *merge_zmq = zmq_socket(context, ZMQ_SUB);
    if (merge) {
        connect_string = "tcp://" + merge_address + ":" + std::to_string(merge_port);
        rc = zmq_connect(merge_zmq, connect_string.c_str());
        assert(rc == 0);
        zmq_setsockopt(merge_zmq, ZMQ_SUBSCRIBE, "", 0);
    }


    // create an rtlsdr device

    int r, opt;
    int ppm_error = 0;
    int sync_mode = 0;
    int dev_index = 0;

    if(out_block_size < MINIMAL_BUF_LENGTH ||
       out_block_size > MAXIMAL_BUF_LENGTH ){
        fprintf(stderr,
            "Output block size wrong value, falling back to default\n");
        fprintf(stderr,
            "Minimal length: %u\n", MINIMAL_BUF_LENGTH);
        fprintf(stderr,
            "Maximal length: %u\n", MAXIMAL_BUF_LENGTH);
        out_block_size = DEFAULT_BUF_LENGTH;
    }

    rtl_buffer = (uint8_t*)malloc(out_block_size * sizeof(uint8_t));

    if (vm.count("device") > 0) {
        dev_index = verbose_device_search(dev_given.c_str());
    }

    if (dev_index < 0) {
        exit(1);
    }

    r = rtlsdr_open(&dev, (uint32_t)dev_index);
    if (r < 0) {
        fprintf(stderr, "Failed to open rtlsdr device #%d.\n", dev_index);
        exit(1);
    }

    /* Set the sample rate */
    verbose_set_sample_rate(dev, (uint32_t)rate);

    rate = (double)rtlsdr_get_sample_rate(dev);

    /* Set the frequency */
    verbose_set_frequency(dev, (uint32_t)freq);

    freq =  (double)rtlsdr_get_center_freq(dev);

    if (0 == gain) {
         /* Enable automatic gain */
        verbose_auto_gain(dev);
    } else {
        /* Enable manual gain */
        gain = nearest_gain(dev, gain*10);
        verbose_gain_set(dev, gain);
        gain = gain/10;
    }

    verbose_ppm_set(dev, ppm_error);

    /* Bias Tee */
    rtlsdr_set_bias_tee(dev, bias_tee);

    /* Reset endpoint before we start reading from it (mandatory) */
    verbose_reset_buffer(dev);

    // Sleep setup time
    std::this_thread::sleep_for(std::chrono::milliseconds(int64_t(1000 * setup_time)));

    struct timeval time_now{};
    gettimeofday(&time_now, nullptr);

    // seed random generator with seconds and microseconds
    srand(time_now.tv_usec + time_now.tv_sec);

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
   
    bool first_frame = true;

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


    if (int_second) {
    	gettimeofday(&time_now, nullptr);
    	struct timeval new_time{};
    	gettimeofday(&new_time, nullptr);
        while (time_now.tv_sec==new_time.tv_sec)
        	gettimeofday(&new_time, nullptr);
    }

    // Run this loop until either time expired (if a duration was given), until
    // the requested number of samples were collected (if such a number was
    // given), or until Ctrl-C was pressed.

    uint32_t frame_count = 0;

    uint32_t len;
    int32_t status=0;
    uint32_t num_words_read=0;

    int16_t bodydata[samps_per_buff*2];

    // Create a circular buffer with a capacity for xxx.
	boost::circular_buffer<int8_t> cb(samps_per_buff*3*2);

    // flush merge queue
    if (merge)
        while ( zmq_recv(merge_zmq, buffer, 100000, ZMQ_NOBLOCK) > 0 ) { }

    while (not stop_signal_called) {
 
        const auto now = std::chrono::steady_clock::now();

        r = rtlsdr_read_sync(dev, rtl_buffer, out_block_size, &n_read);
        if (r < 0) {
            fprintf(stderr, "WARNING: sync read failed.\n");
            break;
        } else {
        	// stuff
        	num_words_read = n_read/2; 
        	// printf("buffer bytes read: %u\n", n_read);
        	for (uint32_t i = 0; i < n_read; i++) {
        		cb.push_back( (int16_t)rtl_buffer[i]-127);
        	}

        	while (cb.size() > 2*samps_per_buff ) {
	        	// printf("circular: %li\n", cb.size());	

	        	gettimeofday(&time_now, nullptr);

		        if (first_frame) {
		            std::cout << boost::format(
		                             "First frame: %u samples, %u full secs, %.09f frac secs")
		                             % n_read % time_now.tv_sec
		                             % (time_now.tv_usec/1e6)
		                      << std::endl;
		            first_frame = false;
		        }

	        	int8_t this_sample;
		        for (uint32_t i = 0; i < 2*samps_per_buff; i++) {
		        	this_sample = cb.front();
        			cb.pop_front();
        			bodydata[i] = this_sample;
        		}	
	   
		        num_total_samps += num_words_read;

		        p.body = bodydata;
		        p.header.packet_count = (uint8_t)frame_count%16;
		        p.fields.integer_seconds_timestamp = time_now.tv_sec;
		        p.fields.fractional_seconds_timestamp = 1e6*time_now.tv_usec;
		
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

		            gettimeofday(&time_now, nullptr);
		            pc.fields.integer_seconds_timestamp = time_now.tv_sec;
		            pc.fields.fractional_seconds_timestamp = 1e3*time_now.tv_usec;

		            pc.fields.stream_id = p.fields.stream_id;

		            pc.if_context.bandwidth                         = rate;
		            pc.if_context.sample_rate                       = rate;
		            pc.if_context.rf_reference_frequency            = freq+if_freq;
		            pc.if_context.rf_reference_frequency_offset     = 0;
		            pc.if_context.if_reference_frequency            = if_freq; // 0 for Zero-IF
		            pc.if_context.gain.stage1                       = gain;
		            pc.if_context.gain.stage2                       = 0;

		            pc.if_context.state_and_event_indicators.reference_lock = false;

		            pc.if_context.state_and_event_indicators.calibrated_time = false;

		            int32_t rv = vrt_write_packet(&pc, buffer, SIZE, true);
		            if (rv < 0) {
		                fprintf(stderr, "Failed to write packet: %s\n", vrt_string_error(rv));
		            }

		            // ZMQ
                    zmq_send (zmq_server, buffer, rv*4, 0);

		        }
	        
		        if (bw_summary) {
		            last_update_samps += num_words_read;
		            const auto time_since_last_update = now - last_update;
		            if (time_since_last_update > std::chrono::seconds(1)) {

		                const double time_since_last_update_s =
		                    std::chrono::duration<double>(time_since_last_update).count();
		                const double rate = double(last_update_samps) / time_since_last_update_s;
		                std::cout << "\t" << (rate / 1e6) << " Msps, ";
		                
		                last_update_samps = 0;
		                last_update       = now;

		                float sum_i = 0;
		                uint32_t clip_i = 0;

		                double datatype_max = 128.;

		                for (int i=0; i<10000; i++ ) {
		                    auto sample_i = get_abs_val(bodydata[2*i]);
		                    sum_i += sample_i;
		                    if (sample_i > datatype_max*0.99)
		                        clip_i++;
		                }
		                sum_i = sum_i/10000;
		                std::cout << boost::format("%.0f") % (100.0*log2(sum_i)/log2(datatype_max)) << "% I (";
		                std::cout << boost::format("%.0f") % ceil(log2(sum_i)+1) << " of ";
		                std::cout << (int)ceil(log2(datatype_max)+1) << " bits), ";
		                std::cout << "" << boost::format("%.0f") % (100.0*clip_i/10000) << "% I clip.";
		                std::cout << std::endl;

		            }
		        }
		    }
		}

        // Merge
        if (merge) {
            int mergelen;
            while ( (mergelen = zmq_recv(merge_zmq, buffer, 100000, ZMQ_NOBLOCK)) > 0  ) {
                // zmq_send (zmq_server, buffer, mergelen, 0);
                zmq_msg_t msg;
                zmq_msg_init_size (&msg, mergelen);
                memcpy (zmq_msg_data(&msg), buffer, mergelen);
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
    }
  
     /* clean up */
    rtlsdr_set_bias_tee(dev, false);
    rtlsdr_close(dev);

    // finished
    std::cout << std::endl << "Done!" << std::endl << std::endl;

    return EXIT_SUCCESS;
}
