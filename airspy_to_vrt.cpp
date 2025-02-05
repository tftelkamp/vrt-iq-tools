//
// Copyright 2025 by Thomas Telkamp 
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

// #include <airspyhf.h>
#include <airspy.h>

// VRT tools functions
#include "vrt-tools.h"

#define DEFAULT_VGA_IF_GAIN (5)
#define DEFAULT_LNA_GAIN (1)
#define DEFAULT_MIXER_GAIN (5)

#define VGA_GAIN_MAX (15)
#define MIXER_GAIN_MAX (15)
#define LNA_GAIN_MAX (14)
#define LINEARITY_GAIN_MAX (21)
#define SENSITIVITY_GAIN_MAX (21)

unsigned long long num_total_samps = 0;

namespace po = boost::program_options;

// Create a circular buffer with a capacity for xxx.
boost::circular_buffer<int16_t> cb(65536*2*4);

boost::shared_mutex _access;

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

// From Phil Karn ka9q
double true_freq(uint64_t freq_hz){
  const int VCO_MIN=1770000000u; // 1.77 GHz
  const int VCO_MAX=3900000000u; // 3.54 GHz
  const int MAX_DIV = 5;

  // Clock divider set to 2 for the best resolution
  // const uint32_t pll_ref = 25000000u/2; // 12.5 MHz

  // For external ref (10 MHz) this works only if I set pll_ref to 10e6, but the R820T clock is 25 MHz. Why?
  const uint32_t pll_ref = 10000000u; // 10 MHz ???? 
 
  // Find divider to put VCO = f*2^(d+1) in range VCO_MIN to VCO_MAX
  //          MHz             step, Hz
  // 0: 885.0     1770.0      190.735
  // 1: 442.50     885.00      95.367
  // 2: 221.25     442.50      47.684
  // 3: 110.625    221.25      23.842
  // 4:  55.3125   110.625     11.921
  // 5:  27.65625   55.312      5.960
  int8_t div_num;
  for (div_num = 0; div_num <= MAX_DIV; div_num++){
    uint32_t vco = freq_hz << (div_num + 1);
    if (VCO_MIN <= vco && vco <= VCO_MAX)
      break;
  }
  if(div_num > MAX_DIV)
    return 0; // Frequency out of range

  // PLL programming bits: Nint in upper 16 bits, Nfract in lower 16 bits
  // Freq steps are pll_ref / 2^(16 + div_num) Hz
  // Note the '+ (pll_ref >> 1)' term simply rounds the division to the nearest integer
  uint32_t r = (((uint64_t) freq_hz << (div_num + 16)) + (pll_ref >> 1)) / pll_ref;
 
  // This is a puzzle; is it related to spur suppression?
  double offset = 0;
  switch(div_num){
  default: // 2, 1, 0
   case 3:
    offset = 0.25;
    break;
  case 4:
    offset = 0.5;
    break;
  case 5:
    offset = 1.0;
    break;
  }

  // Compute true frequency
  return ((double)(r + offset) * pll_ref) / (double)(1 << (div_num + 16));
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
    uint32_t ints_to_write;
    int16_t* pt_rx_buffer;
    
    ints_to_write = transfer->sample_count * 2;
    pt_rx_buffer = (int16_t*)transfer->samples;

    {
        boost::unique_lock< boost::shared_mutex > lock(_access);

        if(pt_rx_buffer) {
            for (uint32_t i = 0; i < ints_to_write; i++) {
                cb.push_back(pt_rx_buffer[i]);
            }
        }
    }
    return 0;
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
    double real_freq;

    // Airspy
    struct airspy_device* device = 0;
    uint64_t serial_number_val;
    uint32_t *supported_samplerates;
    uint32_t sample_rate_u32 = 0;
    uint32_t nsrates;

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
        ("gain", po::value<int>(&gain)->default_value(0), "(sensitivity) gain for the RF chain")
        ("setup", po::value<double>(&setup_time)->default_value(1.0), "seconds of setup time")
        ("vga", po::value<uint32_t>(&vga_gain)->default_value(DEFAULT_VGA_IF_GAIN), "VGA gain for the RF chain")
        ("lna", po::value<uint32_t>(&lna_gain)->default_value(DEFAULT_LNA_GAIN), "VGA gain for the RF chain")
        ("mixer", po::value<uint32_t>(&mixer_gain)->default_value(DEFAULT_MIXER_GAIN), "VGA gain for the RF chain")
        ("sensitivity", po::value<uint32_t>(&sensitivity_gain_val), "Sensitivity gain for the RF chain")
        ("linearity", po::value<uint32_t>(&linearity_gain_val), "Linearity gain for the RF chain")
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

    fprintf(stderr, "%f MS/s %s\n", sample_rate_u32 * 0.000001f, "IQ");

    airspy_read_partid_serialno_t read_partid_serialno;

    if (airspy_board_partid_serialno_read(device, &read_partid_serialno) != AIRSPY_SUCCESS) {
        fprintf(stderr, "airspy_board_partid_serialno_read() failed\n");
        exit(1);
    } else
        fprintf(stderr, "Device Serial Number: 0x%08X%08X\n",
                read_partid_serialno.serial_no[2],
                read_partid_serialno.serial_no[3]);


    result = airspy_set_packing(device, 0);
    if( result != AIRSPY_SUCCESS ) {
        fprintf(stderr, "airspy_set_packing() failed: %s (%d)\n", airspy_error_name((airspy_error)result), result);
        airspy_close(device);
        airspy_exit();
        return EXIT_FAILURE;
    }

    result = airspy_set_rf_bias(device, (int)bias_tee);
    if( result != AIRSPY_SUCCESS ) {
        fprintf(stderr, "airspy_set_rf_bias() failed: %s (%d)\n", airspy_error_name((airspy_error)result), result);
        airspy_close(device);
        airspy_exit();
        return EXIT_FAILURE;
    }

    // check external ref
    uint8_t register_value;
    bool ref = false;
    bool ref_signal = false;

    result = airspy_si5351c_read(device, 0, &register_value); 
    if (result == AIRSPY_SUCCESS)
        if ((register_value & 16) == 0)
            ref_signal = true;

    printf("External reference signal: %s\n", ref_signal ? "on" : "off"); 

    result = airspy_si5351c_read(device, 15, &register_value); 
    if (result == AIRSPY_SUCCESS)
        if ((register_value & 0x0C) != 0)
            ref = true;

    printf("External reference used: %s\n", ref ? "yes" : "no"); 

    // Frequency
    if (ref) {
        real_freq = true_freq(freq);
        printf("True frequency (R820T): %f Hz\n",real_freq);
    } else
        real_freq = freq;

    if( (linearity_gain == false) && (sensitivity_gain == false) )
    {
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
        gain = vga_gain + mixer_gain + lna_gain;
    } else
    {
        if( linearity_gain == true )
        {
            result =  airspy_set_linearity_gain(device, linearity_gain_val);
            if( result != AIRSPY_SUCCESS ) {
                fprintf(stderr, "airspy_set_linearity_gain() failed: %s (%d)\n", airspy_error_name((airspy_error)result), result);
            }
            gain = linearity_gain_val;
        }

        if( sensitivity_gain == true )
        {
            result =  airspy_set_sensitivity_gain(device, sensitivity_gain_val);
            if( result != AIRSPY_SUCCESS ) {
                fprintf(stderr, "airspy_set_sensitivity_gain() failed: %s (%d)\n", airspy_error_name((airspy_error)result), result);
            }
            gain = sensitivity_gain_val;
        }
    }

    if( airspy_start_rx(device, rx_callback, NULL) != AIRSPY_SUCCESS ) {
        fprintf(stderr, "airspy_start() failed.\n");
        exit(1);
    }

    if( airspy_set_freq(device, (uint32_t)freq) != AIRSPY_SUCCESS ) {
        fprintf(stderr, "airspy_set_freq() failed.\n");
        exit(1);
    }

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


    // flush merge queue
    if (merge)
        while ( zmq_recv(merge_zmq, buffer, 100000, ZMQ_NOBLOCK) > 0 ) { }

    while (not stop_signal_called) {
 
        const auto now = std::chrono::steady_clock::now();

        {
        boost::shared_lock< boost::shared_mutex > lock(_access);
    	while (cb.size() > 2*samps_per_buff ) {

        	gettimeofday(&time_now, nullptr);

	        if (first_frame) {
	            std::cout << boost::format(
	                             "First frame: %u samples, %u full secs, %.09f frac secs")
	                             % 0 % time_now.tv_sec
	                             % (time_now.tv_usec/1e6)
	                      << std::endl;
	            first_frame = false;
	        }

            num_words_read = samps_per_buff;

        	int16_t this_sample;
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
	            pc.if_context.rf_reference_frequency            = real_freq+if_freq;
	            pc.if_context.rf_reference_frequency_offset     = 0;
	            pc.if_context.if_reference_frequency            = if_freq; // 0 for Zero-IF
	            pc.if_context.gain.stage1                       = gain;
	            pc.if_context.gain.stage2                       = 0;

	            pc.if_context.state_and_event_indicators.reference_lock = ref;

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

	                double datatype_max = 32767.;

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

    return EXIT_SUCCESS;
}
