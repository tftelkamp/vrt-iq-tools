#include <zmq.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/thread/thread.hpp>

#include <chrono>
// #include <complex>
#include <csignal>
#include <fstream>
#include <iostream>
#include <thread>

// VRT
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include <vrt/vrt_read.h>
#include <vrt/vrt_string.h>
#include <vrt/vrt_types.h>
#include <vrt/vrt_util.h>

#include <complex.h>
// #include <fftw3.h>

// TCP
#include <unistd.h>
#include <arpa/inet.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <sys/time.h>
#include <netinet/in.h>
#include <fcntl.h>

// RTL-SDR
#define closesocket close
#define SOCKADDR struct sockaddr
#define SOCKET int
#define SOCKET_ERROR -1

enum RTL_TCP_COMMANDS {
    SET_FREQUENCY             = 0x01,   /* sets frequency - amending hi word of SET_FREQ_HI32 if present */
    SET_FREQ_HI32             = 0x56,   /* in addition to SET_FREQUENCY */
    SET_SAMPLE_RATE           = 0x02,
    SET_GAIN_MODE             = 0x03,
    SET_GAIN                  = 0x04,
    SET_FREQUENCY_CORRECTION  = 0x05,
    SET_IF_STAGE              = 0x06,
    SET_TEST_MODE             = 0x07,
    SET_AGC_MODE              = 0x08,
    SET_DIRECT_SAMPLING       = 0x09,
    SET_OFFSET_TUNING         = 0x0A,
    SET_RTL_CRYSTAL           = 0x0B,
    SET_TUNER_CRYSTAL         = 0x0C,
    SET_TUNER_GAIN_BY_INDEX   = 0x0D,
#if 1
    /* development branch since 2018-10-03 */
    SET_BIAS_TEE              = 0x0E,
    SET_TUNER_BANDWIDTH       = 0x40,
#else
    /* prev code - used in ExtIO - to build compatible rtl_tcp.exe */
    SET_TUNER_BANDWIDTH       = 0x0E,
    SET_BIAS_TEE              = 0x0F
#endif
    UDP_ESTABLISH             = 0x41,
    UDP_TERMINATE             = 0x42,
    SET_I2C_TUNER_REGISTER    = 0x43,   /* for experiments: 32 bit data word:
                                         * 31 .. 20: register (12 bits)
                                         * 19 .. 12: mask (8 bits)
                                         * 11 ..  0: data (12 bits) */
    SET_I2C_TUNER_OVERRIDE    = 0x44,   /* encoding as with SET_I2C_TUNER_REGISTER
                                         * data (bits 11 .. 0) > 255 removes override */
    SET_TUNER_BW_IF_CENTER    = 0x45,   /* freq from SET_FREQUENCY stays in center;
                                         * the bandwidth (from SET_TUNER_BANDWIDTH)
                                         * is set to be centered at given IF frequency */
    SET_TUNER_IF_MODE         = 0x46,   /* set tuner IF mode - or gain */
    SET_SIDEBAND              = 0x47,   /* Mixer Sideband for R820T */
    REPORT_I2C_REGS           = 0x48,   /* perodically report I2C registers
                                         * - if reverse channel is enabled */

    GPIO_SET_OUTPUT_MODE      = 0x49,   /* rtlsdr_set_gpio_output() */
    GPIO_SET_INPUT_MODE       = 0x50,   /* rtlsdr_set_gpio_input() */
    GPIO_GET_IO_STATUS        = 0x51,   /* rtlsdr_set_gpio_status() */
    GPIO_WRITE_PIN            = 0x52,   /* rtlsdr_set_gpio_output() and rtlsdr_set_gpio_bit() */
    GPIO_READ_PIN             = 0x53,   /* rtlsdr_get_gpio_bit() */
    GPIO_GET_BYTE             = 0x54,   /* rtlsdr_get_gpio_byte() */
    
    IS_TUNER_PLL_LOCKED       = 0x55,   /* rtlsdr_is_tuner_PLL_locked() */

    /* SET_FREQ_HI32          = 0x56,    * rtlsdr_set_center_freq64() */
};

#include "difi-tools.h"

namespace po = boost::program_options;

static SOCKET s;

typedef struct { /* structure size must be multiple of 2 bytes */
    char magic[4];
    uint32_t tuner_type;
    uint32_t tuner_gain_count;
} dongle_info_t;


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
    std::string file, type, zmq_address, rtl_address;
    uint16_t port, rtl_port;
    float scale;
    int hwm;
    size_t num_requested_samples;
    double total_time;

    // setup the program options
    po::options_description desc("Allowed options");
    // clang-format off

    desc.add_options()
        ("help", "help message")
        // ("file", po::value<std::string>(&file)->default_value("usrp_samples.dat"), "name of the file to write binary samples to")
        // ("type", po::value<std::string>(&type)->default_value("short"), "sample type: double, float, or short")
        ("nsamps", po::value<size_t>(&num_requested_samples)->default_value(0), "total number of samples to receive")
        ("duration", po::value<double>(&total_time)->default_value(0), "total number of seconds to receive")
        ("progress", "periodically display short-term bandwidth")
        ("int-second", "align start of reception to integer second")
        ("null", "run without writing to file")
        ("continue", "don't abort on a bad packet")
        ("scale", po::value<float>(&scale)->default_value(1.0), "scaling factor for 16 to 8 bit conversion (default 1)")
        ("address", po::value<std::string>(&zmq_address)->default_value("localhost"), "DIFI ZMQ address")
        ("rtl-address", po::value<std::string>(&rtl_address)->default_value("0.0.0.0"), "RTL-TCP address (default 0.0.0.0)")
        ("port", po::value<uint16_t>(&port)->default_value(50100), "DIFI ZMQ port")
        ("rtl-port", po::value<uint16_t>(&rtl_port)->default_value(1234), "RTL-TCP port (default 1234)")
        ("hwm", po::value<int>(&hwm)->default_value(10000), "DIFI ZMQ HWM")
    ;
    // clang-format on
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // print the help message
    if (vm.count("help")) {
        std::cout << boost::format("DIFI samples to nothing. %s") % desc << std::endl;
        std::cout << std::endl
                  << "This application streams data from a DIFI stream "
                     "to nowhwere.\n"
                  << std::endl;
        return ~0;
    }

    bool progress               = vm.count("progress") > 0;
    bool stats                  = vm.count("stats") > 0;
    bool null                   = vm.count("null") > 0;
    bool continue_on_bad_packet = vm.count("continue") > 0;
    bool int_second             = (bool)vm.count("int-second");


    // RTL
    int r;  
    SOCKET listensocket;
    fd_set readfds;
    fd_set writefds;
    socklen_t rlen;
    struct timeval tv = {1,0};
    struct linger ling = {1,0};
    struct sockaddr_in local, remote;

    dongle_info_t dongle_info;

    struct command{
        unsigned char cmd;
        unsigned int param;
    }__attribute__((packed));

    memset(&local,0,sizeof(local));
    local.sin_family = AF_INET;
    local.sin_port = htons(rtl_port);
    local.sin_addr.s_addr = inet_addr(rtl_address.c_str());

    listensocket = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP);

    r = 1;
    setsockopt(listensocket, SOL_SOCKET, SO_REUSEADDR, (char *)&r, sizeof(int));
    setsockopt(listensocket, SOL_SOCKET, SO_LINGER, (char *)&ling, sizeof(ling));
    bind(listensocket,(struct sockaddr *)&local,sizeof(local));

    /* non-blocking socket */
    fcntl(s, F_SETFL, O_NONBLOCK);

    // TODO
    // std::signal(SIGINT, &sig_int_handler);
    // std::cout << "Press Ctrl + C to stop streaming..." << std::endl;

    while(not stop_signal_called) {

        printf("listening...\n");
        printf("Use the device argument 'rtl_tcp=%s:%d' in OsmoSDR "
                   "(gr-osmosdr) source\n"
                   "to receive samples and set "
                   "difi_to_rtl_tcp parameters (frequency, gain, ...).\n",
                   rtl_address.c_str(), rtl_port);
        listen(listensocket,1);

        while(not stop_signal_called) {
            FD_ZERO(&readfds);
            FD_SET(listensocket, &readfds);
            tv.tv_sec = 1;
            tv.tv_usec = 0;
            r = select(listensocket+1, &readfds, NULL, NULL, &tv);
            if(r) {
                rlen = sizeof(remote);
                s = accept(listensocket,(struct sockaddr *)&remote, &rlen);
                break;
            }
        }

        setsockopt(s, SOL_SOCKET, SO_LINGER, (char *)&ling, sizeof(ling));

        printf("client accepted!\n");

        memset(&dongle_info, 0, sizeof(dongle_info));
        memcpy(&dongle_info.magic, "RTL0", 4);

        // 1 = E4000
        // 5 = RTLSDR_TUNER_R820T
        // 6 = RTLSDR_TUNER_R828D
        dongle_info.tuner_type = htonl(5);

        // ????
        dongle_info.tuner_gain_count = htonl(1);

        r = send(s, (const char *)&dongle_info, sizeof(dongle_info), 0);
        if (sizeof(dongle_info) != r)
            printf("failed to send dongle information\n");

        struct command cmd= {0, 0};
        int left, received = 0;

        // Receive commands
        // while(1) {
        //     if(r) {
        //         received = recv(s, (char*)&cmd, sizeof(cmd), 0);
        //          printf("received %d bytes\n", received); 
        //     }
        //     if(received == SOCKET_ERROR)
        //         break;
        // }

        context_type difi_context;
        init_context(&difi_context);

        difi_packet_type difi_packet;

        // std::vector<std::shared_ptr<std::ofstream>> outfiles;
        std::vector<size_t> channel_nums = {0}; // single channel (0)

        // ZMQ
        void *context = zmq_ctx_new();
        void *subscriber = zmq_socket(context, ZMQ_SUB);
        int rc = zmq_setsockopt (subscriber, ZMQ_RCVHWM, &hwm, sizeof hwm);
        std::string connect_string = "tcp://" + zmq_address + ":" + std::to_string(port);
        rc = zmq_connect(subscriber, connect_string.c_str());
        assert(rc == 0);
        zmq_setsockopt(subscriber, ZMQ_SUBSCRIBE, "", 0);

        // time keeping
        auto start_time = std::chrono::steady_clock::now();
        auto stop_time = start_time + std::chrono::milliseconds(int64_t(1000 * total_time));

        uint32_t buffer[ZMQ_BUFFER_SIZE];

        uint8_t rtlbuffer[DIFI_SAMPLES_PER_PACKET*2];
        
        unsigned long long num_total_samps = 0;

        // Track time and samps between updating the BW summary
        auto last_update                     = start_time;
        unsigned long long last_update_samps = 0;

        bool first_frame = true;
        uint64_t last_fractional_seconds_timestamp = 0;

        // set to true to process data before context
        bool start_rx = false;

        uint32_t signal_pointer = 0;

        while (not stop_signal_called
               and (num_requested_samples > num_total_samps or num_requested_samples == 0)
               and (total_time == 0.0 or std::chrono::steady_clock::now() <= stop_time)) {

            int len = zmq_recv(subscriber, buffer, ZMQ_BUFFER_SIZE, 0);

            const auto now = std::chrono::steady_clock::now();

            if (not difi_process(buffer, sizeof(buffer), &difi_context, &difi_packet)) {
                printf("Not a Vita49 packet?\n");
                continue;
            }

            if (not start_rx and difi_packet.context) {
                difi_print_context(&difi_context);
                start_rx = true;
                // Possibly do something with context here
                // difi_context
            }
            
            FD_ZERO(&readfds);
            FD_SET(s, &readfds);
            tv.tv_sec = 0;
            tv.tv_usec = 0;
            r = select(s+1, &readfds, NULL, NULL, &tv);
            if (r) {
                received = recv(s, (char*)&cmd, sizeof(cmd), 0);
                printf("received %d bytes\n", received); 
                if (received == SOCKET_ERROR || received == 0) {
                    printf("exit\n");
                    break;
                }
                int32_t tmp = ntohl(cmd.param);
                // if (cmd.cmd == SET_IF_STAGE) {
                //     printf("set if stage %d gain %.1f dB\n", tmp >> 16, ((short)(tmp & 0xffff))/10.0);
                // }
                if (cmd.cmd == SET_GAIN) {
                    // tmp += 10;
                    printf("set manual scaling gain %.2f dB (%.1f)\n", tmp/10.0, pow(10,tmp/100.0));
                    scale = pow(10,tmp/100.0);
                }

            }

            if (start_rx and difi_packet.data) {

                if (difi_packet.lost_frame)
                   if (not continue_on_bad_packet)
                        break;

                if (int_second) {
                    // check if fractional second has wrapped
                    if (difi_packet.fractional_seconds_timestamp > last_fractional_seconds_timestamp ) {
                            last_fractional_seconds_timestamp = difi_packet.fractional_seconds_timestamp;
                            continue;
                    } else {
                        int_second = false;
                        last_update = now; 
                        start_time = now;
                        stop_time = start_time + std::chrono::milliseconds(int64_t(1000 * total_time));
                    }
                }

                // Process data here
                // Assumes ci16_le

                for (uint32_t i = 0; i < difi_packet.num_rx_samps; i++) {

                    int16_t re;
                    memcpy(&re, (char*)&buffer[difi_packet.offset+i], 2);
                    int16_t img;
                    memcpy(&img, (char*)&buffer[difi_packet.offset+i]+2, 2);

                    re = std::lroundf((float)re/scale);
                    img = std::lroundf((float)img/scale);

                    re += 128;
                    img += 128;

                    // Clip
                    re = re < 0 ? 0 : re;
                    re = re > 255 ? 255 : re;

                    img = img < 0 ? 0 : img;
                    img = img > 255 ? 255 : img;
                    
                    rtlbuffer[i*2] = (uint8_t)re;
                    rtlbuffer[i*2+1] = (uint8_t)img;

                }

                int bytesleft,bytessent;

                bytesleft = difi_packet.num_rx_samps*sizeof(std::complex<uint8_t>);

                FD_ZERO(&writefds);
                FD_SET(s, &writefds);
                r = select(s+1, NULL, &writefds, NULL, &tv);
                if(r) {
                    bytessent = send(s,  (char*)rtlbuffer, bytesleft, 0);
                }
                if(bytessent == SOCKET_ERROR) {
                        printf("worker socket bye\n");
                        break;
                }

                // data: (const char*)&buffer[difi_packet.offset]
                // size (bytes): sizeof(uint32_t)*difi_packet.num_rx_samps
                 
                num_total_samps += difi_packet.num_rx_samps;

                if (start_rx and first_frame) {
                    std::cout << boost::format(
                                     "  First frame: %u samples, %u full secs, %.09f frac secs")
                                     % difi_packet.num_rx_samps
                                     % difi_packet.integer_seconds_timestamp
                                     % ((double)difi_packet.fractional_seconds_timestamp/1e12)
                              << std::endl;
                    first_frame = false;
                }
            }

            if (progress) {
                if (difi_packet.data)
                    last_update_samps += difi_packet.num_rx_samps;
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

                    for (int i=0; i<difi_packet.num_rx_samps; i++ ) {
                        auto sample_i = get_abs_val((std::complex<int16_t>)buffer[difi_packet.offset+i]);
                        sum_i += sample_i;
                        if (sample_i > datatype_max*0.99)
                            clip_i++;
                    }
                    sum_i = sum_i/difi_packet.num_rx_samps;
                    std::cout << boost::format("%.0f") % (100.0*log2(sum_i)/log2(datatype_max)) << "% I (";
                    std::cout << boost::format("%.0f") % ceil(log2(sum_i)+1) << " of ";
                    std::cout << (int)ceil(log2(datatype_max)+1) << " bits), ";
                    std::cout << "" << boost::format("%.0f") % (100.0*clip_i/difi_packet.num_rx_samps) << "% I clip, ";
                    std::cout << std::endl;

                }
            }
        }
        zmq_close(subscriber);
        zmq_ctx_destroy(context);
    }
    printf("Done\n");

    fflush(stdout);

    return 0;

}  
