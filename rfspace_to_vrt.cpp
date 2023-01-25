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

#include <sys/time.h>

#include <zmq.h>
#include <assert.h>

// VRT
#include <vrt/vrt_init.h>
#include <vrt/vrt_string.h>
#include <vrt/vrt_types.h>
#include <vrt/vrt_util.h>
#include <vrt/vrt_write.h>
#include <vrt/vrt_read.h>

// TCP/UDP
#include <sys/socket.h> 
#include <arpa/inet.h> 
#include <netinet/in.h> 
#include <netinet/tcp.h>
#include <netinet/udp.h>
#include <netdb.h>

#include <boost/thread/thread.hpp>

// DIFI tools functions
#include "difi-tools.h"

unsigned long long num_total_samps = 0;

int sockfd, connfd;

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

bool transaction(const unsigned char *cmd, size_t size, std::vector<unsigned char> &response)
{
    size_t rx_bytes = 0;
    unsigned char data[1024*2];

    response.clear();

    if ( write(sockfd, cmd, size) != (int)size )
      return false;

    int nbytes = read(sockfd, data, 2); /* read header */
    if ( nbytes != 2 )
      return false;

    int length = (data[1] & 0x1f) | data[0];

    if ( (length < 2) || (length > (int)sizeof(data)) )
      return false;

    length -= 2; /* subtract header size */

    nbytes = read(sockfd, &data[2], length); /* read payload */
    if ( nbytes != length )
      return false;

    rx_bytes = 2 + length; /* header + payload */

    response.resize( rx_bytes );
    memcpy( response.data(), data, rx_bytes );

#if 0
  printf("> ");
  for (size_t i = 0; i < rx_bytes; i++)
    printf("%02x ", (unsigned char) data[i]);
  printf("\n");
#endif

  return true;
}

int main(int argc, char* argv[])
{
    // variables to be set by po
    std::string udp_forward, ref, sdrhost;
    size_t total_num_samps;
    uint16_t port;
    uint32_t stream_id;
    int hwm;
    int16_t gain;
    double rate, freq, bw, total_time, setup_time, lo_offset;

    bool context_changed = true;

    // setup the program options
    po::options_description desc("Allowed options");
    // clang-format off
    desc.add_options()
        ("help", "help message")
        ("freq", po::value<double>(&freq)->default_value(0.0), "RF center frequency in Hz")
        // ("lo-offset", po::value<double>(&lo_offset)->default_value(0.0),
        //     "Offset for frontend LO in Hz (optional)")
        ("rate", po::value<double>(&rate)->default_value(1e6), "rate of incoming samples")
        ("gain", po::value<int16_t>(&gain), "gain for the RF chain")
        ("bw", po::value<double>(&bw), "analog frontend filter bandwidth in Hz")
        ("pps", "use external pps signal")
        ("ref", po::value<std::string>(&ref)->default_value("internal"), "reference source (internal, external)")
        ("setup", po::value<double>(&setup_time)->default_value(1.0), "seconds of setup time")
        ("udp", po::value<std::string>(&udp_forward), "DIFI UDP forward address")
        ("progress", "periodically display short-term bandwidth")
        ("stats", "show average bandwidth on exit")
        // ("vrt", "publish IQ using VRT over ZeroMQ (PUB on port 50100")
        ("int-second", "align start of reception to integer second")
        ("null", "run without streaming")
        ("continue", "don't abort on a bad packet")
        ("skip-lo", "skip checking LO lock status")
        ("sdr", po::value<std::string>(&sdrhost)->default_value("cloudsdr"), "RFSPACE SDR hostname")
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
        std::cout << boost::format("RFSpace SDR RX samples to Vita49 %s") % desc << std::endl;
        std::cout << std::endl
                  << "This application streams data from a single channel of an RFSpace SDR "
                     "device to Vita49.\n"
                  << std::endl;
        return ~0;
    }

    bool bw_summary             = vm.count("progress") > 0;
    bool stats                  = vm.count("stats") > 0;
    bool null                   = vm.count("null") > 0;
    bool continue_on_bad_packet = vm.count("continue") > 0;
    // bool vrt                    = vm.count("vrt") > 0;
    // bool zmq                    = vm.count("zmq") > 0;
    bool enable_udp             = vm.count("udp") > 0;

    // SETUP
    #define HEADER_SIZE 2
    #define SEQNUM_SIZE 2

    #define DEFAULT_PORT 50000

    struct sockaddr_in servaddr, cli;
 
    // socket create and verification
    sockfd = socket(AF_INET, SOCK_STREAM, 0);
    if (sockfd == -1) {
        printf("socket creation failed...\n");
        exit(0);
    }

    int sockoptval = 1;
    setsockopt(sockfd, SOL_SOCKET, SO_REUSEADDR, &sockoptval, sizeof(int));
    sockoptval = 1;
    setsockopt(sockfd, IPPROTO_TCP, TCP_NODELAY, &sockoptval, sizeof(int));

    /* don't wait when shutting down */
    linger lngr;
    lngr.l_onoff  = 1;
    lngr.l_linger = 0;
    setsockopt(sockfd, SOL_SOCKET, SO_LINGER, &lngr, sizeof(linger));

    bzero(&servaddr, sizeof(servaddr));

    struct hostent *host;

    if ((host = gethostbyname(sdrhost.c_str())) == NULL)
    {
        printf("SDR hostname not known.\n");
        exit(0);
    }
 
    // assign IP, PORT
    servaddr.sin_family = AF_INET;
    servaddr.sin_addr.s_addr = *(long *)(host->h_addr_list[0]);
    servaddr.sin_port = htons(DEFAULT_PORT);
 
    // connect the client socket to server socket
    if (connect(sockfd, (struct sockaddr*)&servaddr, sizeof(servaddr))
        != 0) {
        printf("connection with the SDR failed.\n");
        exit(0);
    }
    else {
        printf("connected to the SDR.\n");
    }

    // Setup UDP for SDR streaming
    int udp_sockfd; 
 
    // Creating socket file descriptor 
    if ( (udp_sockfd = socket(AF_INET, SOCK_DGRAM, 0)) < 0 ) { 
        perror("UDP socket creation failed"); 
        exit(0); 
    }   

    sockoptval = 1;
    setsockopt(udp_sockfd, SOL_SOCKET, SO_REUSEADDR, &sockoptval, sizeof(int));

    int size = 20*1024*1024; // 20 MB
    setsockopt(udp_sockfd, SOL_SOCKET, SO_RCVBUF, (char *)&size, sizeof(size));  

    /* fill in the hosts's address and data */
    struct sockaddr_in host_sa; /* local address */
    memset(&host_sa, 0, sizeof(host_sa));
    host_sa.sin_family = AF_INET;
    host_sa.sin_addr.s_addr = INADDR_ANY;
    host_sa.sin_port = htons(DEFAULT_PORT);

    if ( bind(udp_sockfd, (const struct sockaddr *)&host_sa, sizeof(host_sa)) < 0 ){
        printf("UDP bind failed.\n");
        exit(0);
    } else {
        printf("UDP port open.\n");
    }

    /* Wait 10 ms before sending queries to device (required for networked radios). */
    std::this_thread::sleep_for(std::chrono::milliseconds(10));

    std::vector<unsigned char> response;

    unsigned char name[] = { 0x04, 0x20, 0x01, 0x00 };
    if ( transaction( name, sizeof(name), response ) )
      std::cerr << "SDR type: " << &response[sizeof(name)] << " " << std::endl;

    unsigned char sern[] = { 0x04, 0x20, 0x02, 0x00 };
    if ( transaction( sern, sizeof(sern), response ) )
      std::cerr << "Serial number: " << &response[sizeof(sern)] << " " << std::endl;;

    // Sample rate
    unsigned char samprate[] = { 0x09, 0x00, 0xB8, 0x00, 0x00, 0x20, 0xA1, 0x07, 0x00 };

    uint32_t u32_rate = rate;
    samprate[sizeof(samprate)-4] = u32_rate >>  0;
    samprate[sizeof(samprate)-3] = u32_rate >>  8;
    samprate[sizeof(samprate)-2] = u32_rate >> 16;
    samprate[sizeof(samprate)-1] = u32_rate >> 24;

    if ( ! transaction( samprate, sizeof(samprate), response ) )
    throw std::runtime_error("set_sample_rate failed");

    u32_rate = 0;
    u32_rate |= response[sizeof(samprate)-4] <<  0;
    u32_rate |= response[sizeof(samprate)-3] <<  8;
    u32_rate |= response[sizeof(samprate)-2] << 16;
    u32_rate |= response[sizeof(samprate)-1] << 24;

    if ( rate != u32_rate )
    std::cerr << "Sample rate set to " << (uint32_t)u32_rate << " Hz"
              << std::endl;

    // Freq
    uint32_t u32_freq = freq;

    unsigned char tune[] = { 0x0A, 0x00, 0x20, 0x00, 0x00, 0xb0, 0x19, 0x6d, 0x00, 0x00 };

    tune[sizeof(tune)-5] = u32_freq >>  0;
    tune[sizeof(tune)-4] = u32_freq >>  8;
    tune[sizeof(tune)-3] = u32_freq >> 16;
    tune[sizeof(tune)-2] = u32_freq >> 24;
    tune[sizeof(tune)-1] = 0;

    transaction( tune, sizeof(tune), response );

    unsigned char getfreq[] = { 0x05, 0x20, 0x20, 0x00, 0x00 };

    if ( ! transaction( getfreq, sizeof(getfreq), response ) )
    throw std::runtime_error("get_center_freq failed");

    uint32_t frequency = 0;
    frequency |= response[response.size()-5] <<  0;
    frequency |= response[response.size()-4] <<  8;
    frequency |= response[response.size()-3] << 16;
    frequency |= response[response.size()-2] << 24;

    printf("Frequency set to: %u\n",frequency);
    
    // Gain

    unsigned char atten[] = { 0x06, 0x00, 0x38, 0x00, 0x00, 0x00 };

    if ( gain <= -30 )
      atten[sizeof(atten)-1] = 0xE2;
    else if ( gain <= -20 )
      atten[sizeof(atten)-1] = 0xEC;
    else if ( gain <= -10 )
      atten[sizeof(atten)-1] = 0xF6;
    else /* 0 dB */
      atten[sizeof(atten)-1] = 0x00;

    transaction( atten, sizeof(atten), response);

    struct timeval time_now{};
    gettimeofday(&time_now, nullptr);

    // seed random generator with seconds and microseconds
    srand(time_now.tv_usec + time_now.tv_sec);

 	// Receive

	size_t samps_per_buff = 10000; // spb

	unsigned long long num_requested_samples = total_num_samps;
    double time_requested = total_time;
    bool int_second             = (bool)vm.count("int-second");

    uint32_t buffer[DIFI_DATA_PACKET_SIZE];
   
    bool first_frame = true;

    struct vrt_packet p;
    vrt_init_packet(&p);

    /* Warn if not standards compliant */
    if (vrt_is_platform_little_endian()) {
        printf("Warning: little endian support is work in progress.\n");
    }

    /* DIFI init */
    difi_init_data_packet(&p);
    
    // Only 1 channel
    p.fields.stream_id = 1;
    
    // ZMQ
    void *zmq_server;
    void *zmq_control;
    void *context = zmq_ctx_new();

    // Stream
    void *responder = zmq_socket(context, ZMQ_PUB);
    int rc = zmq_setsockopt (responder, ZMQ_SNDHWM, &hwm, sizeof hwm);
    assert(rc == 0);

    std::string connect_string = "tcp://*:" + std::to_string(port);
    rc = zmq_bind(responder, connect_string.c_str());
    assert (rc == 0);
    zmq_server = responder;

    // Control
    responder = zmq_socket(context, ZMQ_SUB);
    rc = zmq_bind(responder, "tcp://*:50300");
    assert (rc == 0);
    zmq_control = responder;
    zmq_setsockopt(zmq_control, ZMQ_SUBSCRIBE, "", 0);

    // UDP DI-FI

    int difi_sockfd; 
    struct sockaddr_in difi_servaddr, difi_cliaddr; 
    if (enable_udp) {

        printf("Enable UDP\n");
            
        // Creating socket file descriptor 
        if ( (difi_sockfd = socket(AF_INET, SOCK_DGRAM, 0)) < 0 ) { 
            perror("socket creation failed"); 
            exit(EXIT_FAILURE); 
        } 
            
        memset(&difi_servaddr, 0, sizeof(difi_servaddr)); 
        memset(&difi_cliaddr, 0, sizeof(difi_cliaddr)); 
            
        // Filling server information 
        difi_servaddr.sin_family    = AF_INET; // IPv4 
        difi_servaddr.sin_addr.s_addr = inet_addr(udp_forward.c_str()); /* Server's Address   */
        difi_servaddr.sin_port = htons(50000);  // 4991?
    }

    // Sleep setup time
    std::this_thread::sleep_for(std::chrono::milliseconds(int64_t(1000 * setup_time)));

    if (total_num_samps == 0) {
        std::signal(SIGINT, &sig_int_handler);
        std::cout << "Press Ctrl + C to stop streaming..." << std::endl;
    }

    // time keeping
    auto start_time = std::chrono::steady_clock::now();
    auto stop_time =
        start_time + std::chrono::milliseconds(int64_t(1000 * time_requested));
    
    // Track time and samps between updating the BW summary
    auto last_update                     = start_time;
    auto last_context                    = start_time - std::chrono::milliseconds(2*VRT_CONTEXT_INTERVAL);;
    auto last_keepalive                  = start_time;
    unsigned long long last_update_samps = 0;

    // Hack
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
    uint32_t num_words_read=0;
    uint32_t last_num_rx = 0;

    int16_t bodydata[samps_per_buff*2];

    // Create a circular buffer with a capacity for tbd.
	boost::circular_buffer<int16_t> cb(samps_per_buff*4*2);

    unsigned char data[1024*2];
    uint16_t prev_sequence = 0;

    // Trigger
    if (vm.count("pps")) {
        unsigned char trigger[] = { 0x08, 0x00, 0xb4, 0x00, 0x00, 0x02, 0x00, 0x00 };
        transaction( trigger, sizeof(trigger), response );
    } else {
        unsigned char trigger[] = { 0x08, 0x00, 0xb4, 0x00, 0x00, 0x00, 0x00, 0x00 };
        transaction( trigger, sizeof(trigger), response );
    }

    struct sockaddr_in sa_in;           /* remote address */
    socklen_t addrlen = sizeof(sa_in);  /* length of addresses */

    int64_t sample_counter = 0;
    auto time_first_sample = time_now;

    // START
    unsigned char start[] = { 0x08, 0x00, 0x18, 0x00, 0x80, 0x02, 0x00, 0x00 };
    bool started = transaction( start, sizeof(start), response );

    printf("Started: %s\n", started ? "yes" : "no");

    while (not stop_signal_called) {
 
        const auto now = std::chrono::steady_clock::now();

        // RX
        ssize_t rx_bytes = recvfrom(udp_sockfd, data, sizeof(data), 0, (struct sockaddr *)&sa_in, &addrlen);

        gettimeofday(&time_now, nullptr);

        // No need to check, we configure the SDR for 16 bit
        // if ( (0x04 == data[0] && (0x84 == data[1] || 0x82 == data[1])) )
        // {
        //     //    is_24_bit = false;
        // }
        // else if ( (0xA4 == data[0] && 0x85 == data[1]) ||
        //         (0x84 == data[0] && 0x81 == data[1]) )
        // {
        //     //    is_24_bit = true;
        //     printf("24 bit data not yet supported.\n");
        //     break;
        // }

        uint16_t sequence = *((uint16_t *)(data + HEADER_SIZE));

        uint16_t diff = sequence - prev_sequence;

        if ( diff > 1 )
        {
            std::cerr << "Lost " << (diff-1) << " packets from " << sequence << " "
                  << inet_ntoa(sa_in.sin_addr) << ":" << ntohs(sa_in.sin_port)
                  << std::endl;
            if (not continue_on_bad_packet)
                break;
            // insert zero's in circular buffer to make up for lost packet(s)
            for (uint32_t p = 0; p < (diff-1); p++)
                for (uint32_t i = 0; i < 2*last_num_rx; i++) {
                    cb.push_back(0);
                }
        }

        prev_sequence = (0xffff == sequence) ? 0 : sequence;

        /* get pointer to samples */
        uint16_t *sample = (uint16_t *)(data + HEADER_SIZE + SEQNUM_SIZE);

        size_t rx_samples = (rx_bytes - HEADER_SIZE - SEQNUM_SIZE) / (sizeof(int16_t) * 2);

        last_num_rx = rx_samples;

    	for (uint32_t i = 0; i < 2*rx_samples; i++) {
    		cb.push_back( (int16_t)(sample[i]) );
    	}

        sample_counter += rx_samples;

        const auto time_since_last_context = now - last_context;
        if (time_since_last_context > std::chrono::milliseconds(VRT_CONTEXT_INTERVAL)) {

            last_context = now;

            // VITA 49.2
            /* Initialize to reasonable values */
            struct vrt_packet pc;
            vrt_init_packet(&pc);

            /* DIFI Configure. Note that context packets cannot have a trailer word. */
            difi_init_context_packet(&pc);

            if (context_changed)
                pc.if_context.context_field_change_indicator = true;
            else
                pc.if_context.context_field_change_indicator = false;

            gettimeofday(&time_now, nullptr);
            pc.fields.integer_seconds_timestamp = time_now.tv_sec;
            pc.fields.fractional_seconds_timestamp = 1e3*time_now.tv_usec;

            pc.fields.stream_id = p.fields.stream_id;

            pc.if_context.bandwidth                         = u32_rate;
            pc.if_context.sample_rate                       = u32_rate;
            pc.if_context.rf_reference_frequency            = frequency;
            pc.if_context.rf_reference_frequency_offset     = 0;
            pc.if_context.if_reference_frequency            = 0; // Zero-IF
            pc.if_context.gain.stage1                       = gain;
            pc.if_context.gain.stage2                       = 0;

            pc.if_context.state_and_event_indicators.reference_lock = (bool)(ref=="external") ? true : false;

            pc.if_context.state_and_event_indicators.calibrated_time = (bool)(vm.count("pps")) ? true : false;

            int32_t rv = vrt_write_packet(&pc, buffer, DIFI_DATA_PACKET_SIZE, true);
            if (rv < 0) {
                fprintf(stderr, "Failed to write packet: %s\n", vrt_string_error(rv));
            }

            // ZMQ
            zmq_send (zmq_server, buffer, rv*4, 0);

            if (enable_udp) {
                if (sendto(difi_sockfd, buffer, rv*4, 0,
                     (struct sockaddr *)&servaddr, sizeof(servaddr)) < 0)
                {
                   printf("UDP fail\n");
                }
            }
            context_changed = false;
        }

        if (first_frame) {
                std::cout << boost::format(
                                 "First frame: %u samples, %u full secs, %.09f frac secs (counter %u)")
                                 % (rx_samples) % time_now.tv_sec
                                 % (time_now.tv_usec/1e6)
                                 % sequence
                          << std::endl;
                first_frame = false;
                time_first_sample = time_now;
                // If PPS is enabled, the first packet should be on the PPS
                if ( vm.count("pps") && (sequence==0) )
                    time_first_sample.tv_usec = 0;
        }

    	while (cb.size() > 2*samps_per_buff ) {

        	gettimeofday(&time_now, nullptr);

            num_words_read = samps_per_buff;

            struct timeval interval_time, difi_time;
            int64_t first_sample = sample_counter - cb.size()/2;

            double interval = (double)first_sample/(double)u32_rate;
            interval_time.tv_sec = (time_t)interval;
            interval_time.tv_usec = (interval-(time_t)interval)*1e6;

            timeradd(&time_first_sample, &interval_time, &difi_time);

	        for (uint32_t i = 0; i < 2*samps_per_buff; i++) {
    			bodydata[i] = (int16_t)cb.front();
                cb.pop_front();
    		}	
   
	        num_total_samps += num_words_read;

	        p.body = bodydata;
	        p.header.packet_count = (uint8_t)frame_count%16;
	        p.fields.integer_seconds_timestamp = difi_time.tv_sec;
	        p.fields.fractional_seconds_timestamp = 1e6*difi_time.tv_usec;
	
	        zmq_msg_t msg;
	        int rc = zmq_msg_init_size (&msg, DIFI_DATA_PACKET_SIZE*4);

	        int32_t rv = vrt_write_packet(&p, zmq_msg_data(&msg), DIFI_DATA_PACKET_SIZE, true);

	        // UDP
	        if (enable_udp) {
	            if (sendto(sockfd, zmq_msg_data(&msg), DIFI_DATA_PACKET_SIZE*4, 0,
	                         (struct sockaddr *)&difi_servaddr, sizeof(difi_servaddr)) < 0)
	            {
	               printf("UDP fail\n");
	            }
	        }

	        zmq_msg_send(&msg, zmq_server, 0);

	        zmq_msg_close(&msg);

            frame_count++;

            // Control 
            int len = zmq_recv(zmq_control, buffer, 100000, ZMQ_NOBLOCK);
            if (len > 0) {
                printf("-> Control context received\n");

                struct vrt_header h;
                struct vrt_fields f;

                int32_t offset = 0;
                int32_t size = ZMQ_BUFFER_SIZE;
                int32_t rv = vrt_read_header(buffer + offset, size - offset, &h, true);

                /* Parse header */
                if (rv < 0) {
                    fprintf(stderr, "Failed to parse header: %s\n", vrt_string_error(rv));
                    break;
                }
                offset += rv;

                if (h.packet_type == VRT_PT_IF_CONTEXT) {
                    // Context

                    /* Parse fields */
                    rv = vrt_read_fields(&h, buffer + offset, size - offset, &f, true);
                    if (rv < 0) {
                        fprintf(stderr, "Failed to parse fields section: %s\n", vrt_string_error(rv));
                        break;
                    }
                    offset += rv;

                    struct vrt_if_context c;
                    rv = vrt_read_if_context(buffer + offset, ZMQ_BUFFER_SIZE - offset, &c, true);
                    if (rv < 0) {
                        fprintf(stderr, "Failed to parse IF context section: %s\n", vrt_string_error(rv));
                        break;
                    }

                    // Channel
                    uint32_t ch = 0;
                    while(not (f.stream_id & (1 << ch) ) )
                        ch++;

                    printf("    Channel: %u\n", ch);

                    // Freq
                    if (c.has.rf_reference_frequency) {

                        freq = (double)round(c.rf_reference_frequency);
                        uint32_t u32_freq = freq;

                        unsigned char tune[] = { 0x0A, 0x00, 0x20, 0x00, 0x00, 0xb0, 0x19, 0x6d, 0x00, 0x00 };

                        tune[sizeof(tune)-5] = u32_freq >>  0;
                        tune[sizeof(tune)-4] = u32_freq >>  8;
                        tune[sizeof(tune)-3] = u32_freq >> 16;
                        tune[sizeof(tune)-2] = u32_freq >> 24;
                        tune[sizeof(tune)-1] = 0;

                        transaction( tune, sizeof(tune), response );

                        unsigned char getfreq[] = { 0x05, 0x20, 0x20, 0x00, 0x00 };

                        if ( ! transaction( getfreq, sizeof(getfreq), response ) )
                        throw std::runtime_error("get_center_freq failed");

                        frequency = 0;
                        frequency |= response[response.size()-5] <<  0;
                        frequency |= response[response.size()-4] <<  8;
                        frequency |= response[response.size()-3] << 16;
                        frequency |= response[response.size()-2] << 24;

                        printf("    Frequency set to: %u\n",frequency);
                    }

                    if (c.has.gain) {
                        // Gain
                        gain = c.gain.stage1;
                        int set_gain = 0;
                        unsigned char atten[] = { 0x06, 0x00, 0x38, 0x00, 0x00, 0x00 };

                        if ( gain <= -30 ) {
                          atten[sizeof(atten)-1] = 0xE2;
                          set_gain = -30;
                        }
                        else if ( gain <= -20 ) {
                          atten[sizeof(atten)-1] = 0xEC;
                          set_gain = -20;
                        }
                        else if ( gain <= -10 ) {
                          atten[sizeof(atten)-1] = 0xF6;
                          set_gain = -10;
                        }
                        else /* 0 dB */ {
                          atten[sizeof(atten)-1] = 0x00;
                          set_gain = 0;
                        }
                        printf("    Gain set to: %i\n",set_gain);

                        transaction( atten, sizeof(atten), response);
                    }

                    last_context = start_time - std::chrono::milliseconds(2*VRT_CONTEXT_INTERVAL); // Trigger context update (next)

                    context_changed = true;

                }
            }

            const auto time_since_last_keepalive = now - last_keepalive;
            if (time_since_last_keepalive > std::chrono::seconds(60)) {
                unsigned char status_pkt[] = { 0x04, 0x20, 0x05, 0x00 };
                transaction( status_pkt, sizeof(status_pkt), response );
                last_keepalive = now;
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

	                double datatype_max = 32768.;

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

    const auto actual_stop_time = std::chrono::steady_clock::now();

    unsigned char stop[] = { 0x08, 0x00, 0x18, 0x00, 0x00, 0x01, 0x00, 0x00 };
    transaction( stop, sizeof(stop), response );

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
    // close the socket
    close(sockfd);

    // wait for ZMQ a bit
    std::this_thread::sleep_for(std::chrono::milliseconds(int64_t(1000 * setup_time)));

    // finished
    std::cout << std::endl << "Done!" << std::endl << std::endl;

    return EXIT_SUCCESS;
}
