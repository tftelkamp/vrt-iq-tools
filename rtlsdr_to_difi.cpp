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

// UDP
#include <sys/socket.h> 
#include <arpa/inet.h> 
#include <netinet/in.h> 

#include <boost/thread/thread.hpp>

#include <inttypes.h>

#include "rtl-sdr.h"
#include "convenience.h"

// DIFI tools functions
#include "difi-tools.h"

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
    std::string gain_list, udp_forward;
    size_t total_num_samps;
    uint16_t port;
    uint32_t stream_id;
    int hwm;
    double rate, freq, bw, total_time, setup_time, lo_offset;

    // recv_frame_size=1024, num_recv_frames=1024, recv_buff_size
    std::string stdargs = "num_recv_frames=1024";
    std::string args;

    // setup the program options
    po::options_description desc("Allowed options");
    // clang-format off
    desc.add_options()
        ("help", "help message")
        ("args", po::value<std::string>(&args)->default_value(""), "multi uhd device address args")
        ("rate", po::value<double>(&rate)->default_value(1e6), "rate of incoming samples")
        ("freq", po::value<double>(&freq)->default_value(0.0), "RF center frequency in Hz")
        // ("lo-offset", po::value<double>(&lo_offset)->default_value(0.0),
        //     "Offset for frontend LO in Hz (optional)")
        ("gain", po::value<std::string>(&gain_list), "gain(s) for the RF chain")
        ("bw", po::value<double>(&bw), "analog frontend filter bandwidth in Hz")
        ("setup", po::value<double>(&setup_time)->default_value(1.0), "seconds of setup time")
        ("udp", po::value<std::string>(&udp_forward), "DIFI UDP forward address")
        ("progress", "periodically display short-term bandwidth")
        ("stats", "show average bandwidth on exit")
        ("vrt", "publish IQ using VRT over ZeroMQ (PUB on port 50100")
        ("int-second", "align start of reception to integer second")
        ("null", "run without streaming")
        ("continue", "don't abort on a bad packet")
        ("skip-lo", "skip checking LO lock status")
        // ("stream-id", po::value<uint32_t>(&stream_id), "DIFI Stream ID")
        ("port", po::value<uint16_t>(&port)->default_value(50100), "DIFI ZMQ port")
        ("hwm", po::value<int>(&hwm)->default_value(10000), "DIFI ZMQ HWM")
    ;
    // clang-format on
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // print the help message
    if (vm.count("help")) {
        std::cout << boost::format("RTL-SDR RX samples to Vita49 %s") % desc << std::endl;
        std::cout << std::endl
                  << "This application streams data from a single channel of an RTL-SDR "
                     "device to Vita49.\n"
                  << std::endl;
        return ~0;
    }

    bool bw_summary             = vm.count("progress") > 0;
    bool stats                  = vm.count("stats") > 0;
    bool null                   = vm.count("null") > 0;
    bool continue_on_bad_packet = vm.count("continue") > 0;
    bool vrt                    = vm.count("vrt") > 0;
    bool zmq                    = vm.count("zmq") > 0;
    bool enable_udp             = vm.count("udp") > 0;

    vrt = true;

    // create an rtlsdr device

    int r, opt;
    int gain = 0;
    int ppm_error = 0;
    int sync_mode = 0;
    int dev_index = 0;
    int dev_given = 0;
    uint32_t frequency = freq;
    uint32_t samp_rate = rate;

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

    if (!dev_given) {
        dev_index = verbose_device_search((char*)"0");
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
    verbose_set_sample_rate(dev, samp_rate);

    /* Set the frequency */
    verbose_set_frequency(dev, frequency);

    if (0 == gain) {
         /* Enable automatic gain */
        verbose_auto_gain(dev);
    } else {
        /* Enable manual gain */
        gain = nearest_gain(dev, gain);
        verbose_gain_set(dev, gain);
    }

    verbose_ppm_set(dev, ppm_error);

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

	size_t samps_per_buff = 10000; // spb

	unsigned long long num_requested_samples = total_num_samps;
    double time_requested = total_time;
    bool int_second             = (bool)vm.count("int-second");

    uint32_t buffer[SIZE];
   
    bool first_frame = true;

    struct vrt_packet p;
    vrt_init_packet(&p);

    /* Warn if not standards compliant */
    if (vrt_is_platform_little_endian()) {
        printf("Warning: little endian support is work in progress.\n");
    }

    /* DIFI init */
    difi_init_data_packet(&p);
    
    p.fields.stream_id = 1;
    
    // ZMQ
    void *zmq_server;
    if (vrt) {
        void *context = zmq_ctx_new();
        void *responder = zmq_socket(context, ZMQ_PUB);
        int rc = zmq_setsockopt (responder, ZMQ_SNDHWM, &hwm, sizeof hwm);
        assert(rc == 0);

        std::string connect_string = "tcp://*:" + std::to_string(port);
        rc = zmq_bind(responder, connect_string.c_str());
        assert (rc == 0);
        zmq_server = responder;
    }

    // UDP DI-FI

    int sockfd; 
    struct sockaddr_in servaddr, cliaddr; 
    if (enable_udp) {

        printf("Enable UDP\n");
            
        // Creating socket file descriptor 
        if ( (sockfd = socket(AF_INET, SOCK_DGRAM, 0)) < 0 ) { 
            perror("socket creation failed"); 
            exit(EXIT_FAILURE); 
        } 
            
        memset(&servaddr, 0, sizeof(servaddr)); 
        memset(&cliaddr, 0, sizeof(cliaddr)); 
            
        // Filling server information 
        servaddr.sin_family    = AF_INET; // IPv4 
        servaddr.sin_addr.s_addr = inet_addr(udp_forward.c_str()); /* Server's Address   */
        servaddr.sin_port = htons(50000);  // 4991?
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
        		cb.push_back( (int16_t)rtl_buffer[i]-127 );
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

		        // UDP
		        if (enable_udp) {
		            if (sendto(sockfd, zmq_msg_data(&msg), SIZE*4, 0,
		                         (struct sockaddr *)&servaddr, sizeof(servaddr)) < 0)
		            {
		               printf("UDP fail\n");
		            }
		        }

		        // VRT
		        if (vrt) {
		            zmq_msg_send(&msg, zmq_server, 0);
		        }

		        zmq_msg_close(&msg);

		        const auto time_since_last_context = now - last_context;
		        if (time_since_last_context > std::chrono::milliseconds(200)) {

		            last_context = now;

		            // VITA 49.2
		            /* Initialize to reasonable values */
		            struct vrt_packet pc;
		            vrt_init_packet(&pc);

		            /* DIFI Configure. Note that context packets cannot have a trailer word. */
		            difi_init_context_packet(&pc);

		            gettimeofday(&time_now, nullptr);
		            pc.fields.integer_seconds_timestamp = time_now.tv_sec;
		            pc.fields.fractional_seconds_timestamp = 1e3*time_now.tv_usec;

		            pc.fields.stream_id = p.fields.stream_id;

		            pc.if_context.bandwidth                         = rate;
		            pc.if_context.sample_rate                       = rate;
		            pc.if_context.rf_reference_frequency            = freq;
		            pc.if_context.rf_reference_frequency_offset     = 0;
		            pc.if_context.if_reference_frequency            = 0; // Zero-IF
		            pc.if_context.gain.stage1                       = gain;
		            pc.if_context.gain.stage2                       = 0;

		            pc.if_context.state_and_event_indicators.reference_lock = false;

		            pc.if_context.state_and_event_indicators.calibrated_time = false;

		            int32_t rv = vrt_write_packet(&pc, buffer, SIZE, true);
		            if (rv < 0) {
		                fprintf(stderr, "Failed to write packet: %s\n", vrt_string_error(rv));
		            }

		            // ZMQ
		            if (vrt) {
		                  zmq_send (zmq_server, buffer, rv*4, 0);
		            }

		            if (enable_udp) {
		                if (sendto(sockfd, buffer, rv*4, 0,
		                     (struct sockaddr *)&servaddr, sizeof(servaddr)) < 0)
		                {
		                   printf("UDP fail\n");
		                }
		            }
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
    rtlsdr_close(dev);

    // finished
    std::cout << std::endl << "Done!" << std::endl << std::endl;

    return EXIT_SUCCESS;
}
