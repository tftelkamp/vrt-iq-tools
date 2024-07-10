//
// Copyright 2010-2011,2014 Ettus Research LLC
// Copyright 2018 Ettus Research, a National Instruments Company
//
// Copyright 2021/2022 by Thomas Telkamp
//
// SPDX-License-Identifier: GPL-3.0-or-later
//

#include <unistd.h>
#include <uhd/exception.hpp>
#include <uhd/types/tune_request.hpp>
#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/utils/safe_main.hpp>
#include <uhd/utils/thread.hpp>
#include <uhd/types/sensors.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
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

// VRT tools functions
#include "vrt-tools.h"

unsigned long long num_total_samps = 0;

namespace po = boost::program_options;

static bool stop_signal_called = false;
void sig_int_handler(int)
{
    stop_signal_called = true;
}

static UHD_INLINE boost::uint32_t GPIO_BIT(const size_t x)
{
    return (1 << x);
}

typedef std::function<uhd::sensor_value_t(const std::string&)> get_sensor_fn_t;

bool check_locked_sensor(std::vector<std::string> sensor_names,
    const char* sensor_name,
    get_sensor_fn_t get_sensor_fn,
    double setup_time)
{
    if (std::find(sensor_names.begin(), sensor_names.end(), sensor_name)
        == sensor_names.end())
        return false;

    auto setup_timeout = std::chrono::steady_clock::now()
                         + std::chrono::milliseconds(int64_t(setup_time * 1000));
    bool lock_detected = false;

    std::cout << boost::format("Waiting for \"%s\": ") % sensor_name;
    std::cout.flush();

    while (true) {
        if (lock_detected and (std::chrono::steady_clock::now() > setup_timeout)) {
            std::cout << " locked." << std::endl;
            break;
        }
        if (get_sensor_fn(sensor_name).to_bool()) {
            std::cout << "+";
            std::cout.flush();
            lock_detected = true;
        } else {
            if (std::chrono::steady_clock::now() > setup_timeout) {
                std::cout << std::endl;
                throw std::runtime_error(
                    str(boost::format(
                            "timed out waiting for consecutive locks on sensor \"%s\"")
                        % sensor_name));
            }
            std::cout << "_";
            std::cout.flush();
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
    std::cout << std::endl;
    return true;
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


void transmit_worker(uhd::usrp::multi_usrp::sptr usrp,
                     uhd::tx_streamer::sptr tx_streamer,
                     void *zmq_transmit,
                     double tx_freq,
                     double tx_lo_offset,
                     double tx_gain,
                     double sample_rate,
                     bool enable_gpio,
                     double gpio_delay)
   {

    std::vector<std::complex<short>> buff(VRT_SAMPLES_PER_PACKET);
    std::vector<std::complex<short>*> buffs(1, &buff.front());

    size_t samps_per_buff = VRT_SAMPLES_PER_PACKET; // spb
    uint32_t tx_zmq_buffer[VRT_DATA_PACKET_SIZE];

    uhd::tx_metadata_t metadata;
    metadata.start_of_burst = true;
    metadata.end_of_burst   = false;
    metadata.has_time_spec  = false;

    std::string gpio = "FP0";

    //configure GPIO
    uint32_t gpio_bit = 0;
    uint32_t duplex_bit = 1;
    size_t num_bits = 6;
    boost::uint32_t mask = (1 << num_bits) - 1;

    float gpio_start_delay = gpio_delay/1000.0;
    float gpio_stop_delay = gpio_delay/1000.0;

    if (enable_gpio) {
        //set data direction register (DDR)
        usrp->set_gpio_attr(gpio, "DDR", (GPIO_BIT(gpio_bit)|GPIO_BIT(duplex_bit)), mask);

        //set control register
        usrp->set_gpio_attr(gpio, "CTRL", GPIO_BIT(duplex_bit), mask);

        //set output values
        usrp->set_gpio_attr(gpio, "OUT", 0, mask);

        // ATR Duplex for testing on GPIO duplex_bit
        usrp->set_gpio_attr(gpio, "ATR_XX", GPIO_BIT(duplex_bit), mask);
    }

    // send data until the signal handler gets called
    while (not stop_signal_called) {

        // Receive data
        int len = zmq_recv(zmq_transmit, tx_zmq_buffer, 100000, ZMQ_NOBLOCK);

        if (len > 0) {

            // metadata.time_spec = uhd::time_spec_t(usrp->get_time_now() + uhd::time_spec_t(0.5));

            struct vrt_header h;
            struct vrt_fields f;

            int32_t offset = 0;
            int32_t size = ZMQ_BUFFER_SIZE;
            int32_t rv = vrt_read_header(tx_zmq_buffer + offset, size - offset, &h, true);

            /* Parse header */
            if (rv < 0) {
                fprintf(stderr, "Failed to parse header: %s\n", vrt_string_error(rv));
                break;
            }
            offset += rv;

            if (h.packet_type == VRT_PT_IF_DATA_WITH_STREAM_ID) {

                /* Parse fields */
                rv = vrt_read_fields(&h, tx_zmq_buffer + offset, size - offset, &f, true);
                if (rv < 0) {
                    fprintf(stderr, "Failed to parse fields section: %s\n", vrt_string_error(rv));
                    break;
                }
                offset += rv;

                // Add check for missing packets

                uint32_t num_rx_samps = (h.packet_size-offset);

                uint32_t stream_id = f.stream_id;

                if (num_rx_samps <= VRT_SAMPLES_PER_PACKET) {

                    for (uint32_t i = 0; i < num_rx_samps; i++) {
                        int16_t re;
                        memcpy(&re, (char*)&tx_zmq_buffer[offset+i], 2);
                        int16_t img;
                        memcpy(&img, (char*)&tx_zmq_buffer[offset+i]+2, 2);

                        buff[i] = std::complex<short>(re, img);
                    }

                    // send the entire contents of the ZMQ buffer
                    tx_streamer->send(buffs, num_rx_samps, metadata);

                    metadata.start_of_burst = false;
                    metadata.has_time_spec  = false;
                    metadata.end_of_burst   = false;
                }

            } else if (h.packet_type == VRT_PT_IF_CONTEXT) {
                // Context

                /* Parse fields */
                rv = vrt_read_fields(&h, tx_zmq_buffer + offset, size - offset, &f, true);
                if (rv < 0) {
                    fprintf(stderr, "Failed to parse fields section: %s\n", vrt_string_error(rv));
                    break;
                }
                offset += rv;

                struct vrt_if_context c;
                rv = vrt_read_if_context(tx_zmq_buffer + offset, ZMQ_BUFFER_SIZE - offset, &c, true);
                if (rv < 0) {
                    fprintf(stderr, "Failed to parse IF context section: %s\n", vrt_string_error(rv));
                    break;
                }

                if (c.context_field_change_indicator) {

                    double lo_offset;
                    if (c.has.if_band_offset) {
                        lo_offset = c.if_band_offset;
                    } else {
                        lo_offset = tx_lo_offset;
                    }

                    if (c.has.rf_reference_frequency) {
                        if (tx_freq != (double)round(c.rf_reference_frequency)) {
                            tx_freq = (double)round(c.rf_reference_frequency);
                            std::cout << boost::format("    Setting TX Freq: %f MHz...") % (tx_freq / 1e6)
                                      << std::endl;
                            std::cout << boost::format("    Setting TX LO Offset: %f MHz...") % (lo_offset / 1e6)
                                      << std::endl;
                            uhd::tune_request_t tune_request(tx_freq, lo_offset);
                            // if (vm.count("int-n"))
                            //     tune_request.args = uhd::device_addr_t("mode_n=integer");
                            usrp->set_tx_freq(tune_request);
                            std::cout << boost::format("    Actual TX Freq: %f MHz...")
                                             % (usrp->get_tx_freq() / 1e6)
                                      << std::endl;
                        }
                    }
                    if (c.has.gain) {
                        if (tx_gain != c.gain.stage1) {
                            tx_gain = c.gain.stage1;
                            std::cout << boost::format("    Setting TX Gain: %f dB...") % tx_gain << std::endl;
                            usrp->set_tx_gain(tx_gain);
                            std::cout << boost::format("    Actual TX Gain: %f dB...")
                                             % usrp->get_tx_gain()
                                      << std::endl;
                        }
                    }
                }

                if (c.state_and_event_indicators.user_defined == 0x1) {
                    if (c.state_and_event_indicators.has.calibrated_time && c.state_and_event_indicators.calibrated_time) {

                        timeval vrt_time;
                        vrt_time.tv_sec = f.integer_seconds_timestamp;
                        vrt_time.tv_usec = f.fractional_seconds_timestamp/1e6;
                        uhd::time_spec_t start_time(vrt_time.tv_sec, (double)vrt_time.tv_usec / 1e6);

                        metadata.has_time_spec = true;
                        metadata.time_spec = start_time;

                        printf("Timed transmit queued (%ld frac %.09f).\n", vrt_time.tv_sec, (double)vrt_time.tv_usec / 1e6);

                        // GPIO
                        if (enable_gpio) {
                            usrp->set_command_time(start_time - uhd::time_spec_t(gpio_start_delay));
                            usrp->set_gpio_attr(gpio, "OUT", GPIO_BIT(gpio_bit), GPIO_BIT(gpio_bit));
                        }
                    } else {
                        if (enable_gpio) {
                            usrp->set_gpio_attr(gpio, "OUT", GPIO_BIT(gpio_bit), GPIO_BIT(gpio_bit));
                            boost::this_thread::sleep_for(boost::chrono::milliseconds((uint32_t)(gpio_start_delay*1000)));
                        }
                        printf("Start transmit.\n");
                    }


                } else if (c.state_and_event_indicators.user_defined == 0x2) {
                    printf("End transmit.\n");
                    metadata.end_of_burst = true;
                    metadata.start_of_burst = false;
                    metadata.has_time_spec  = false;
                    tx_streamer->send("", 0, metadata);
                    metadata.end_of_burst = false;

                    // GPIO
                    if (enable_gpio) {
                        if (c.state_and_event_indicators.has.calibrated_time && c.state_and_event_indicators.calibrated_time) {
                            timeval vrt_time;
                            vrt_time.tv_sec = f.integer_seconds_timestamp;
                            vrt_time.tv_usec = f.fractional_seconds_timestamp/1e6;
                            uhd::time_spec_t stop_time(vrt_time.tv_sec, (double)vrt_time.tv_usec / 1e6);
                            usrp->set_command_time(stop_time + uhd::time_spec_t(gpio_stop_delay));
                            usrp->set_gpio_attr(gpio, "OUT", 0, GPIO_BIT(gpio_bit));
                        } else {
                            // the '5' is a guess of the usrp buffer depth
                            usrp->set_command_time(usrp->get_time_now() + uhd::time_spec_t(5*(float)VRT_SAMPLES_PER_PACKET/sample_rate) + uhd::time_spec_t(gpio_stop_delay));
                            usrp->set_gpio_attr(gpio, "OUT", 0, GPIO_BIT(gpio_bit));
                        }
                    }

                }

            }
        } else {
            boost::this_thread::sleep_for(boost::chrono::microseconds(10));
        }
    }

    // send a mini EOB packet
    metadata.end_of_burst = true;
    metadata.start_of_burst = false;
    metadata.has_time_spec  = false;
    tx_streamer->send("", 0, metadata);
}


int UHD_SAFE_MAIN(int argc, char* argv[])
{
    // variables to be set by po
    std::string file, type, ant_list, subdev, ref, channel_list, gain_list, freq_list, udp_forward, merge_address;
    size_t total_num_samps, spb;
    uint16_t instance, port, merge_port;
    uint16_t tx_gain;
    int hwm, io_threads;
    uint32_t stream_id;
    double rate, freq, bw, total_time, setup_time, lo_offset, tx_freq, if_freq, pps_offset, gpio_delay;
    uint32_t timestamp_calibration_time = 0;

    bool context_changed = true;
    bool merge;

    // recv_frame_size=1024, num_recv_frames=1024, recv_buff_size
    std::string stdargs = "num_recv_frames=1024";
    std::string args;

    // setup the program options
    po::options_description desc("Allowed options");
    // clang-format off
    desc.add_options()
        ("help,h", "help message")
        ("args", po::value<std::string>(&args)->default_value(""), "multi uhd device address args")
        // ("type", po::value<std::string>(&type)->default_value("short"), "sample type: double, float, or short")
        ("nsamps", po::value<size_t>(&total_num_samps)->default_value(0), "total number of samples to receive")
        ("duration", po::value<double>(&total_time)->default_value(0), "total number of seconds to receive")
        // ("spb", po::value<size_t>(&spb)->default_value(10000), "samples per buffer")
        ("rate", po::value<double>(&rate)->default_value(1e6), "rate of incoming samples")
        ("freq", po::value<std::string>(&freq_list)->required(), "RF center frequency (list) in Hz")
        ("if-freq", po::value<double>(&if_freq)->default_value(0.0), "IF center frequency in Hz")
        ("lo-offset", po::value<double>(&lo_offset)->default_value(0.0),
            "Offset for frontend LO in Hz (optional)")
        ("gain", po::value<std::string>(&gain_list), "gain(s) for the RF chain")
        ("ant", po::value<std::string>(&ant_list), "antenna selection")
        ("subdev", po::value<std::string>(&subdev), "subdevice specification")
        ("usrp-channel", po::value<std::string>(&channel_list)->default_value("0"), "which usrp channel(s) to use (specify \"0\", \"1\", \"0,1\", etc)")
        ("zmq-split", "create a ZeroMQ stream per VRT channel, increasing port number for additional streams")
        ("bw", po::value<double>(&bw), "analog frontend filter bandwidth in Hz")
        ("ref", po::value<std::string>(&ref)->default_value("internal"), "reference source (internal, external, mimo, gpsdo)")
        ("tx", "enable tx")
        ("tx-freq", po::value<double>(&tx_freq)->default_value(0.0), "TX RF center frequency in Hz")
        ("tx-gain", po::value<uint16_t>(&tx_gain)->default_value(0), "gain for the TX RF chain")
        ("gpio", "enable GPIO (TX)")
        ("gpio-delay", po::value<double>(&gpio_delay)->default_value(50), "GPIO advance/delay (ms)")
        ("setup", po::value<double>(&setup_time)->default_value(1.0), "seconds of setup time")
        ("udp", po::value<std::string>(&udp_forward), "VRT UDP forward address")
        ("progress", "periodically display short-term bandwidth")
        ("stats", "show average bandwidth on exit")
        ("pps", "use external pps signal")
        ("pps-offset", po::value<double>(&pps_offset)->default_value(0), "Offset of the PPS pulse in sec.")
        ("temp", "read temperature sensor")
        ("int-second", "align start of reception to integer second")
        ("null", "run without streaming")
        ("continue", "don't abort on a bad packet")
        ("skip-lo", "skip checking LO lock status")
        ("int-n", "tune USRP with integer-N tuning")
        ("port", po::value<uint16_t>(&port), "VRT ZMQ port")
        ("instance", po::value<uint16_t>(&instance)->default_value(0), "VRT ZMQ instance")
        ("merge", po::value<bool>(&merge)->default_value(true), "Merge another VRT ZMQ stream (SUB connect)")
        ("merge-port", po::value<uint16_t>(&merge_port)->default_value(50011), "VRT ZMQ merge port")
        ("merge-address", po::value<std::string>(&merge_address)->default_value("localhost"), "VRT ZMQ merge address")
        ("io-threads", po::value<int>(&io_threads)->default_value(1), "ZMQ IO threads")
        ("hwm", po::value<int>(&hwm)->default_value(10000), "VRT ZMQ HWM")
    ;
    // clang-format on
    po::variables_map vm;
    auto parsed = po::command_line_parser(argc, argv).options(desc).positional({}).run();
    po::store(parsed, vm);

    // print the help message
    if (vm.count("help") || argc < 2) {
        std::cout << boost::format("UHD samples to VRT. %s") % desc << std::endl;
        std::cout << std::endl
                  << "This application streams data from a USRP "
                     "device to VRT.\n"
                  << std::endl;
        return ~0;
    }
    po::notify(vm);

    bool bw_summary             = vm.count("progress") > 0;
    bool stats                  = vm.count("stats") > 0;
    bool null                   = vm.count("null") > 0;
    bool continue_on_bad_packet = vm.count("continue") > 0;
    bool enable_udp             = vm.count("udp") > 0;
    bool enable_temp            = vm.count("temp") > 0;
    bool enable_tx              = vm.count("tx") > 0;
    bool enable_gpio            = vm.count("gpio") > 0;
    bool split                  = vm.count("zmq-split") > 0;

    struct vrt_packet p;
    vrt_init_packet(&p);

    /* Warn if not standards compliant */
    if (vrt_is_platform_little_endian()) {
        printf("Warning: little endian support is work in progress.\n");
    }

    /* VRT init */
    vrt_init_data_packet(&p);

    p.fields.stream_id = 0;

    // detect channels to use
    std::vector<std::string> channel_strings;
    std::vector<size_t> channel_nums;
    boost::split(channel_strings, channel_list, boost::is_any_of("\"',"));

    // ZMQ
    void *zmq_server[MAX_CHANNELS];
    void *zmq_control;
    void *zmq_transmit;

    void *context = zmq_ctx_new();
    void *responder;
    int rc;

    zmq_ctx_set (context, ZMQ_IO_THREADS, io_threads);

    uint16_t main_port;

    if (vm.count("port") > 0) {
        main_port = port;
    } else {
        main_port = DEFAULT_MAIN_PORT + MAX_CHANNELS*instance;
    }

    if (split) {
        for (size_t ch = 0; ch < channel_strings.size(); ch++) {
            responder = zmq_socket(context, ZMQ_PUB);
            rc = zmq_setsockopt (responder, ZMQ_SNDHWM, &hwm, sizeof hwm);
            assert(rc == 0);
            std::string connect_string = "tcp://*:" + std::to_string(main_port+ch);
            rc = zmq_bind(responder, connect_string.c_str());
            assert (rc == 0);
            zmq_server[ch] = responder;
        }
    } else {
        responder = zmq_socket(context, ZMQ_PUB);
        rc = zmq_setsockopt (responder, ZMQ_SNDHWM, &hwm, sizeof hwm);
        assert(rc == 0);
        std::string connect_string = "tcp://*:" + std::to_string(main_port);
        rc = zmq_bind(responder, connect_string.c_str());
        assert (rc == 0);
        zmq_server[0] = responder;
    }

    responder = zmq_socket(context, ZMQ_SUB);
    std::string control_string = "tcp://*:" + std::to_string(main_port+200);
    rc = zmq_bind(responder, control_string.c_str());
    assert (rc == 0);
    zmq_control = responder;
    zmq_setsockopt(zmq_control, ZMQ_SUBSCRIBE, "", 0);

    if (enable_tx) {
        responder = zmq_socket(context, ZMQ_SUB);
        std::string tx_string = "tcp://*:" + std::to_string(main_port+400);
        rc = zmq_bind(responder, tx_string.c_str());
        assert (rc == 0);
        zmq_transmit = responder;
        zmq_setsockopt(zmq_transmit, ZMQ_SUBSCRIBE, "", 0);
    }

    // Merge
    void *merge_zmq = zmq_socket(context, ZMQ_SUB);
    if (merge) {
        std::string connect_string = "tcp://" + merge_address + ":" + std::to_string(merge_port);
        rc = zmq_connect(merge_zmq, connect_string.c_str());
        assert(rc == 0);
        zmq_setsockopt(merge_zmq, ZMQ_SUBSCRIBE, "", 0);
    }

    // UDP VRT

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

    // create a usrp device

    if (not vm["args"].defaulted()) {
        args = stdargs + "," + args;
    } else
        args = stdargs;

    std::cout << std::endl;
    std::cout << boost::format("Creating the usrp device with: %s...") % args
              << std::endl;
    uhd::usrp::multi_usrp::sptr usrp = uhd::usrp::multi_usrp::make(args);

    // Lock mboard clocks
    if (vm.count("ref")) {
        usrp->set_clock_source(ref);
    }

    if (ref == "gpsdo") {
        usrp->set_time_source(ref);
    }

   if (vm.count("pps")) {
        usrp->set_time_source("external");
    }

    std::cout << "Clock source is " << usrp->get_clock_source(0) << std::endl;
    std::cout << "Time source is " << usrp->get_time_source(0) << std::endl;

    // always select the subdevice first, the channel mapping affects the other settings
    if (vm.count("subdev"))
        usrp->set_rx_subdev_spec(subdev);

    std::cout << boost::format("Using Device: %s") % usrp->get_pp_string() << std::endl;

    // detect which channels to use
    for (size_t ch = 0; ch < channel_strings.size(); ch++) {
        size_t chan = std::stoi(channel_strings[ch]);
        if (chan >= usrp->get_rx_num_channels()) {
            throw std::runtime_error("Invalid channel(s) specified.");
        } else
            channel_nums.push_back(std::stoi(channel_strings[ch]));
    }

    // set the sample rate
    if (rate <= 0.0) {
        std::cerr << "Please specify a valid sample rate" << std::endl;
        return ~0;
    }
    std::cout << boost::format("Setting RX Rate: %f Msps...") % (rate / 1e6) << std::endl;
    usrp->set_rx_rate(rate);
    std::cout << boost::format("Actual RX Rate: %f Msps...")
                     % (usrp->get_rx_rate() / 1e6)
              << std::endl
              << std::endl;

    if (enable_tx) {
        std::cout << boost::format("Setting TX Rate: %f Msps...") % (rate / 1e6) << std::endl;
        usrp->set_tx_rate(rate);
        std::cout << boost::format("Actual TX Rate: %f Msps...")
                         % (usrp->get_tx_rate() / 1e6)
                  << std::endl
                  << std::endl;
    }

    std::vector<double> frequencies;
    if (vm.count("freq")) {
        std::vector<std::string> freq_strings;
        boost::split(freq_strings, freq_list, boost::is_any_of("\"',"));
        for (size_t ch = 0; ch < freq_strings.size(); ch++) {
            frequencies.push_back(std::stod(freq_strings[ch]));
        }
    }

    std::vector<size_t> gains;
    if (vm.count("gain")) {
        std::vector<std::string> gain_strings;
        boost::split(gain_strings, gain_list, boost::is_any_of("\"',"));
        for (size_t ch = 0; ch < gain_strings.size(); ch++) {
            gains.push_back(std::stoi(gain_strings[ch]));
        }
    }

    std::vector<std::string> antennas;
    if (vm.count("ant")) {
        std::vector<std::string> ant_strings;
        boost::split(ant_strings, ant_list, boost::is_any_of("\"',"));
        for (size_t ch = 0; ch < ant_strings.size(); ch++) {
            antennas.push_back(ant_strings[ch]);
        }
    }

    for (size_t ch = 0; ch < channel_nums.size(); ch++) {
        size_t channel = channel_nums[ch];
        if (channel_nums.size() > 1) {
            std::cout << "Configuring RX Channel " << channel << std::endl;
        }

        // set the center frequency
        if (vm.count("freq")) {
            freq = (frequencies.size() > ch) ? frequencies[ch] : frequencies[0];
            if (freq < 5e6) {
                throw std::runtime_error("Frequency should be given in Hz.\n" +
                                         std::to_string(freq) + "Hz is probably not what you meant!");
            }
            std::cout << boost::format("Setting RX Freq: %f MHz...") % (freq / 1e6)
                      << std::endl;
            std::cout << boost::format("Setting RX LO Offset: %f MHz...") % (lo_offset / 1e6)
                      << std::endl;
            uhd::tune_request_t tune_request(freq, lo_offset);
            if (vm.count("int-n"))
                tune_request.args = uhd::device_addr_t("mode_n=integer");
            usrp->set_rx_freq(tune_request, channel);
            std::cout << boost::format("Actual RX Freq: %f MHz...")
                             % (usrp->get_rx_freq(channel) / 1e6)
                      << std::endl
                      << std::endl;
        }

        // set the rf gain(s)
        if (gains.size()) {
            size_t gain = (gains.size() > ch) ? gains[ch] : gains[0];
            std::cout << boost::format("Setting RX Gain: %f dB...") % gain << std::endl;
            usrp->set_rx_gain(gain, channel);
            std::cout << boost::format("Actual RX Gain: %f dB...")
                             % usrp->get_rx_gain(channel)
                      << std::endl
                      << std::endl;
        }

        // set the IF filter bandwidth
        if (vm.count("bw")) {
            std::cout << boost::format("Setting RX Bandwidth: %f MHz...") % (bw / 1e6)
                      << std::endl;
            usrp->set_rx_bandwidth(bw, channel);
            std::cout << boost::format("Actual RX Bandwidth: %f MHz...")
                             % (usrp->get_rx_bandwidth(channel) / 1e6)
                      << std::endl
                      << std::endl;
        }

        // set the antenna(s)
        if (antennas.size()) {
            std::string ant = (antennas.size() > ch) ? antennas[ch] : antennas[0];
            std::cout << boost::format("Setting Antenna: %s") % ant << std::endl;
            usrp->set_rx_antenna(ant, channel);
            std::cout << boost::format("Actual Antenna: %s")
                             % usrp->get_rx_antenna(channel)
                      << std::endl
                      << std::endl;
        }
    }

    if (enable_tx) {

        // Freq
        if (freq < 5e6) {
            throw std::runtime_error("TX frequency should be given in Hz.\n" +
                                     std::to_string(tx_freq) + "Hz is probably not what you meant!");
        }
        std::cout << boost::format("Setting TX Freq: %f MHz...") % (tx_freq / 1e6)
                  << std::endl;
        std::cout << boost::format("Setting TX LO Offset: %f MHz...") % (lo_offset / 1e6)
                  << std::endl;
        uhd::tune_request_t tune_request(tx_freq, lo_offset);
        if (vm.count("int-n"))
            tune_request.args = uhd::device_addr_t("mode_n=integer");
        usrp->set_tx_freq(tune_request);
        std::cout << boost::format("Actual TX Freq: %f MHz...")
                         % (usrp->get_tx_freq() / 1e6)
                  << std::endl
                  << std::endl;

        // Gain
        std::cout << boost::format("Setting TX Gain: %f dB...") % tx_gain << std::endl;
        usrp->set_tx_gain(tx_gain);
        std::cout << boost::format("Actual TX Gain: %f dB...")
                         % usrp->get_tx_gain()
                  << std::endl
                  << std::endl;

        // std::cout << boost::format("Setting TX Bandwidth: %f MHz...") % (tx_bw / 1e6)
        //           << std::endl;
        // tx_usrp->set_tx_bandwidth(tx_bw);
        std::cout << boost::format("Actual TX Bandwidth: %f MHz...")
                         % (usrp->get_tx_bandwidth() / 1e6)
                  << std::endl
                  << std::endl;
    }

    std::this_thread::sleep_for(std::chrono::milliseconds(int64_t(1000 * setup_time)));

    size_t channel = channel_nums[0];

    // check Ref and LO Lock detect
    if (not vm.count("skip-lo")) {
        check_locked_sensor(usrp->get_rx_sensor_names(channel),
            "lo_locked",
            [usrp, channel](const std::string& sensor_name) {
                return usrp->get_rx_sensor(sensor_name, channel);
            },
            setup_time);
        if (ref == "mimo") {
            check_locked_sensor(usrp->get_mboard_sensor_names(0),
                "mimo_locked",
                [usrp](const std::string& sensor_name) {
                    return usrp->get_mboard_sensor(sensor_name);
                },
                setup_time);
        }
        if (ref == "external") {
            check_locked_sensor(usrp->get_mboard_sensor_names(0),
                "ref_locked",
                [usrp](const std::string& sensor_name) {
                    return usrp->get_mboard_sensor(sensor_name);
                },
                setup_time);
        }
        if (ref == "gpsdo") {
            check_locked_sensor(usrp->get_mboard_sensor_names(0),
                "ref_locked",
                [usrp](const std::string& sensor_name) {
                    return usrp->get_mboard_sensor(sensor_name);
                },
                setup_time);
        }
    }

    // create a receive streamer
    const std::string& cpu_format = "sc16";
    const std::string& wire_format = "sc16";
    uhd::stream_args_t stream_args(cpu_format, wire_format);
    stream_args.channels             = channel_nums;
    uhd::rx_streamer::sptr rx_stream = usrp->get_rx_stream(stream_args);

    // reset usrp time to prepare for transmit/receive
    std::cout << boost::format("Setting device timestamp to current time...") << std::endl;

    struct timeval time_now{};
    gettimeofday(&time_now, nullptr);

    // seed random generator with seconds and microseconds
    srand(time_now.tv_usec + time_now.tv_sec);

    // Non-PPS
    usrp->set_time_now(uhd::time_spec_t(time_now.tv_sec, (double)time_now.tv_usec / 1e6));

    // PPS
    if (vm.count("pps")) {
        uint32_t usrp_seconds;
        do {
            gettimeofday(&time_now, nullptr);
            time_t integer_time = (time_t)((double)time_now.tv_sec + (double)time_now.tv_usec/1e6 + 2.0 - pps_offset);
            uhd::time_spec_t set_pps_time = uhd::time_spec_t(integer_time + pps_offset);
            std::cout << boost::format("Wait for PPS sync...") << std::endl;
            usrp->set_time_unknown_pps(set_pps_time);
            boost::this_thread::sleep_for(boost::chrono::milliseconds(2100));
            gettimeofday(&time_now, nullptr);
            usrp_seconds = usrp->get_time_now().get_full_secs();
        } while (usrp_seconds != time_now.tv_sec);

        timestamp_calibration_time = (uint32_t)usrp_seconds;
        std::cout << boost::format("Done...") << std::endl;
    }

    if (ref=="gpsdo") {
        // Check PPS and compare UHD device time to GPS time
        uhd::sensor_value_t gps_time   = usrp->get_mboard_sensor("gps_time");
        uhd::time_spec_t last_pps_time = usrp->get_time_last_pps();

        // we only care about the full seconds
        signed gps_seconds;
        long long pps_seconds;

        do {
            std::cout << "\nTrying to align the device time to GPS time..." << std::endl;

            gps_time = usrp->get_mboard_sensor("gps_time");

            // set the device time to the GPS time
            // getting the GPS time returns just after the PPS edge, so just add a
            // second and set the device time at the next PPS edge
            usrp->set_time_next_pps(uhd::time_spec_t(gps_time.to_int() + 1.0));
            // allow some time to make sure the PPS has come…
            std::this_thread::sleep_for(std::chrono::milliseconds(1100));
            //…then ask
            gps_seconds = usrp->get_mboard_sensor("gps_time").to_int();
            pps_seconds = usrp->get_time_last_pps().to_ticks(1.0);
        } while (pps_seconds != gps_seconds);

        timestamp_calibration_time = pps_seconds;

        if (pps_seconds == gps_seconds) {
                std::cout << "GPS and UHD Device time are aligned.\n";
            } else {
                std::cout << "Could not align UHD Device time to GPS time. Giving up.\n";
            }
            std::cout << boost::format("last_pps: %ld vs gps: %ld.") % pps_seconds % gps_seconds
                  << std::endl;

        std::cout << boost::format("GPS Epoch time at last PPS: %.6f seconds\n")
                         % usrp->get_mboard_sensor("gps_time").to_real();
        std::cout << boost::format("UHD Device time last PPS:   %.6f seconds\n")
                         % (usrp->get_time_last_pps().get_real_secs());
        std::cout << boost::format("UHD Device time right now:  %.6f seconds\n")
                         % (usrp->get_time_now().get_real_secs());
        gettimeofday(&time_now, nullptr);
        std::cout << boost::format("PC Clock time:              %.6f seconds\n") % (time_now.tv_sec + (double)time_now.tv_usec / 1e6);  //time(NULL);
    } else {
        std::cout << boost::format("UHD Device time last PPS:   %.6f seconds\n")
                         % (usrp->get_time_last_pps().get_real_secs());
        std::cout << boost::format("UHD Device time right now:  %.6f seconds\n")
                         % (usrp->get_time_now().get_real_secs());
        gettimeofday(&time_now, nullptr);
        std::cout << boost::format("PC Clock time:              %.6f seconds\n") % (time_now.tv_sec + (double)time_now.tv_usec / 1e6);
    }

    // TX
    std::thread transmit_thread;
    uhd::tx_streamer::sptr tx_stream;

    if (enable_tx) {
        std::vector<size_t> tx_channel_nums;
        tx_channel_nums.push_back(0);
        uhd::stream_args_t tx_stream_args("sc16", "sc16");
        tx_stream_args.channels = tx_channel_nums;
        tx_stream = usrp->get_tx_stream(tx_stream_args);

        // start thread
        transmit_thread = std::thread([&]() {
            transmit_worker(usrp, tx_stream, zmq_transmit, tx_freq, lo_offset, tx_gain, rate, enable_gpio, gpio_delay);
        });
    }

    // RX

    if (total_num_samps == 0) {
        std::signal(SIGINT, &sig_int_handler);
        std::cout << "Press Ctrl + C to stop streaming..." << std::endl;
    }

    unsigned long long num_requested_samples = total_num_samps;
    bool int_second             = (bool)vm.count("int-second");

    // fixed buffer size
    size_t samps_per_buff = VRT_SAMPLES_PER_PACKET; // spb

    uint32_t buffer[VRT_DATA_PACKET_SIZE];

    uhd::rx_metadata_t md;
    std::vector<std::vector<std::complex<short>>> buffs(
        channel_nums.size(), std::vector<std::complex<short>>(samps_per_buff));
    // create a vector of pointers to point to each of the channel buffers
    std::vector<std::complex<short>*> buff_ptrs;
    for (size_t i = 0; i < buffs.size(); i++) {
        buff_ptrs.push_back(&buffs[i].front());
    }
    UHD_ASSERT_THROW(buffs.size() == channel_nums.size());

    bool overflow_message = true;
    bool first_frame = true;

    // time keeping
    auto start_time = std::chrono::steady_clock::now();
    // auto stop_time =
    //     start_time + std::chrono::milliseconds(int64_t(1000 * time_requested));

    // setup streaming
    // uhd::stream_cmd_t stream_cmd((num_requested_samples == 0)
    //                                  ? uhd::stream_cmd_t::STREAM_MODE_START_CONTINUOUS
    //                                  : uhd::stream_cmd_t::STREAM_MODE_NUM_SAMPS_AND_DONE);
    // stream_cmd.num_samps  = size_t(num_requested_samples);
    uhd::stream_cmd_t stream_cmd(uhd::stream_cmd_t::STREAM_MODE_START_CONTINUOUS);

    if (int_second || channel_nums.size() > 1) {
        stream_cmd.time_spec  = usrp->get_time_now().get_full_secs() + 1;
        stream_cmd.stream_now = false;
    }
    else {
        stream_cmd.time_spec  = uhd::time_spec_t();
        stream_cmd.stream_now = true;
    }
    rx_stream->issue_stream_cmd(stream_cmd);

    // Track time and samps between updating the BW summary
    auto last_update                     = start_time;
    auto last_context                    = start_time - std::chrono::milliseconds(2*VRT_CONTEXT_INTERVAL);
    unsigned long long last_update_samps = 0;

    // if (int_second) {
    //     stop_time += std::chrono::milliseconds(int64_t(1000.0*(double)(stream_cmd.time_spec.get_real_secs()-usrp->get_time_now().get_real_secs())));
    // }

    if (total_time > 0)
        num_requested_samples = total_time * rate;

    // Run this loop until either time expired (if a duration was given), until
    // the requested number of samples were collected (if such a number was
    // given), or until Ctrl-C was pressed.

    uint32_t frame_count = 0;

    // flush merge queue
    if (merge)
        while ( zmq_recv(merge_zmq, buffer, 100000, ZMQ_NOBLOCK) > 0 ) { }

    while (not stop_signal_called
           and (num_requested_samples > num_total_samps or num_requested_samples == 0)) {
           // and (time_requested == 0.0 or std::chrono::steady_clock::now() <= stop_time)) {
        const auto now = std::chrono::steady_clock::now();

        size_t num_rx_samps =
            rx_stream->recv(buff_ptrs, samps_per_buff, md, 3.0, false);

        if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_TIMEOUT) {
            std::cout << boost::format("Timeout while streaming") << std::endl;
            break;
        }
        if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_OVERFLOW) {
            if (overflow_message) {
                // overflow_message = false;
                std::cerr
                    << boost::format(
                           "Got an overflow indication. Host does not consume data fast enough (%fMB/s).\n")
                           % (usrp->get_rx_rate() * sizeof(std::complex<short>) / 1e6);
                if (!continue_on_bad_packet)
                    break;
            }
            continue;
        }
        if (md.error_code != uhd::rx_metadata_t::ERROR_CODE_NONE) {
            std::string error = str(boost::format("Receiver error: %s") % md.strerror());
            if (continue_on_bad_packet) {
                std::cerr << error << std::endl;
                continue;
            } else
                throw std::runtime_error(error);
        }

        if (first_frame) {
            std::cout << boost::format(
                             "First frame: %u samples, %u full secs, %.09f frac secs")
                             % num_rx_samps % md.time_spec.get_full_secs()
                             % md.time_spec.get_frac_secs()
                      << std::endl;
            first_frame = false;
            stream_cmd.stream_now = false;
            last_update = now;
        }

        const auto time_since_last_context = now - last_context;
        if (time_since_last_context > std::chrono::milliseconds(VRT_CONTEXT_INTERVAL)) {

            last_context = now;

            // VITA 49
            /* Initialize to reasonable values */
            struct vrt_packet pc;
            vrt_init_packet(&pc);


            /* VRT Configure. Note that context packets cannot have a trailer word. */
            vrt_init_context_packet(&pc);

            pc.fields.integer_seconds_timestamp = md.time_spec.get_full_secs();
            pc.fields.fractional_seconds_timestamp = (uint64_t)1e12 * md.time_spec.get_frac_secs();

            for (size_t ch = 0; ch < channel_nums.size(); ch++) {
                size_t channel = channel_nums[ch];

                if (enable_temp) {
                    uhd::sensor_value_t temp = usrp->get_rx_sensor("temp", channel);
                    pc.if_context.has.temperature = true;
                    pc.if_context.temperature = temp.to_real();
                }

                if (timestamp_calibration_time != 0) {
                    pc.if_context.has.timestamp_calibration_time = true;
                    pc.if_context.timestamp_calibration_time = timestamp_calibration_time;
                }

                if (context_changed)
                    pc.if_context.context_field_change_indicator = true;
                else
                    pc.if_context.context_field_change_indicator = false;

                if (split)
                    pc.fields.stream_id = 1;
                else
                    pc.fields.stream_id = 1<<ch;
                pc.if_context.bandwidth                         = usrp->get_rx_bandwidth(channel); // 0.8*usrp->get_rx_rate(); // bandwith is set to 80% of sample rate
                pc.if_context.sample_rate                       = usrp->get_rx_rate(channel);
                pc.if_context.rf_reference_frequency            = usrp->get_rx_freq(channel)+if_freq;
                pc.if_context.rf_reference_frequency_offset     = 0;
                pc.if_context.if_reference_frequency            = if_freq; // 0 for Zero-IF
                pc.if_context.if_band_offset                    = lo_offset;  //todo
                pc.if_context.gain.stage1                       = usrp->get_rx_gain(channel);
                pc.if_context.gain.stage2                       = 0;

                pc.if_context.state_and_event_indicators.has.reference_lock = true;
                pc.if_context.state_and_event_indicators.reference_lock = ((ref == "external") or (ref=="gpsdo"));

                pc.if_context.state_and_event_indicators.has.calibrated_time = true;
                pc.if_context.state_and_event_indicators.calibrated_time = ((vm.count("pps")) or (ref=="gpsdo"));

                int32_t rv = vrt_write_packet(&pc, buffer, VRT_DATA_PACKET_SIZE, true);
                if (rv < 0) {
                    fprintf(stderr, "Failed to write packet: %s\n", vrt_string_error(rv));
                }

                // ZMQ
                if (split)
                    zmq_send (zmq_server[ch], buffer, rv*4, 0);
                else
                    zmq_send (zmq_server[0], buffer, rv*4, 0);

                if (enable_udp) {
                    if (sendto(sockfd, buffer, rv*4, 0,
                         (struct sockaddr *)&servaddr, sizeof(servaddr)) < 0)
                    {
                       printf("UDP fail\n");
                    }
                }
            }
            context_changed = false;
        }

        num_total_samps += num_rx_samps;

        p.fields.integer_seconds_timestamp = md.time_spec.get_full_secs();
        p.fields.fractional_seconds_timestamp = (uint64_t)1e12 * md.time_spec.get_frac_secs();
        p.header.packet_count = (uint8_t)frame_count%16;

        for (size_t i = 0; i < buffs.size(); i++) {
            size_t channel = channel_nums[i];
            p.body = (char*)buff_ptrs[i];
            if (split)
                p.fields.stream_id = 1;
                else
                p.fields.stream_id = 1<<i;
            zmq_msg_t msg;
            int rc = zmq_msg_init_size (&msg, VRT_DATA_PACKET_SIZE*4);
            int32_t rv = vrt_write_packet(&p, zmq_msg_data(&msg), VRT_DATA_PACKET_SIZE, true);

            // VRT
            if (split)
                zmq_msg_send(&msg, zmq_server[i], 0);
            else
                zmq_msg_send(&msg, zmq_server[0], 0);
            zmq_msg_close(&msg);

            // UDP
            if (enable_udp) {
                if (sendto(sockfd, zmq_msg_data(&msg), VRT_DATA_PACKET_SIZE*4, 0,
                             (struct sockaddr *)&servaddr, sizeof(servaddr)) < 0)
                {
                   printf("UDP fail\n");
                }
            }
        }

        // Merge
        if (merge) {
            int mergelen;
            while ( (mergelen = zmq_recv(merge_zmq, buffer, 100000, ZMQ_NOBLOCK)) > 0  ) {

                if (split) {
                    for (size_t ch = 0; ch < channel_nums.size(); ch++) {
                        zmq_msg_t msg;
                        zmq_msg_init_size (&msg, mergelen);
                        memcpy (zmq_msg_data(&msg), buffer, mergelen);
                        zmq_msg_send(&msg, zmq_server[ch], 0);
                        zmq_msg_close(&msg);
                    }
                } else {
                    zmq_msg_t msg;
                    zmq_msg_init_size (&msg, mergelen);
                    memcpy (zmq_msg_data(&msg), buffer, mergelen);
                    zmq_msg_send(&msg, zmq_server[0], 0);
                    zmq_msg_close(&msg);
                }
            }
        }

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
                while(not (f.stream_id & (1 << channel_nums[ch]) ) )
                    ch++;

                uint32_t control_channel = channel_nums[ch];

                printf("    Channel: %u\n", control_channel);

                if (c.has.if_band_offset) {
                    lo_offset = c.if_band_offset;
                }

                if (c.has.rf_reference_frequency || c.has.if_band_offset) {
                    if (c.has.rf_reference_frequency)
                        freq = (double)round(c.rf_reference_frequency);
                    std::cout << boost::format("    Setting RX Freq: %f MHz...") % (freq / 1e6)
                              << std::endl;
                    std::cout << boost::format("    Setting RX LO Offset: %f MHz...") % (lo_offset / 1e6)
                              << std::endl;
                    uhd::tune_request_t tune_request(freq, lo_offset);
                    if (vm.count("int-n"))
                        tune_request.args = uhd::device_addr_t("mode_n=integer");
                    usrp->set_rx_freq(tune_request, control_channel);
                    std::cout << boost::format("    Actual RX Freq: %f MHz...")
                                     % (usrp->get_rx_freq(control_channel) / 1e6)
                              << std::endl;
                }
                if (c.has.gain) {
                    double gain = c.gain.stage1;
                    std::cout << boost::format("    Setting RX Gain: %f dB...") % gain << std::endl;
                    usrp->set_rx_gain(gain, control_channel);
                    std::cout << boost::format("    Actual RX Gain: %f dB...")
                                     % usrp->get_rx_gain(control_channel)
                              << std::endl;
                }

                last_context = start_time - std::chrono::milliseconds(2*VRT_CONTEXT_INTERVAL); // Trigger context update (next)

                context_changed = true;

            }
        }

        if (bw_summary) {
            last_update_samps += num_rx_samps;
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
                if (cpu_format == "sc8" || cpu_format == "s8")
                    datatype_max = 128.;

                for (int i=0; i<num_rx_samps; i++ ) {
                    auto sample_i = get_abs_val(buff_ptrs[0][i]);
                    sum_i += sample_i;
                    if (sample_i > datatype_max*0.99)
                        clip_i++;
                }
                sum_i = sum_i/num_rx_samps;
                std::cout << boost::format("%.0f") % (100.0*log2(sum_i)/log2(datatype_max)) << "% I (";
                std::cout << boost::format("%.0f") % ceil(log2(sum_i)+1) << " of ";
                std::cout << (int)ceil(log2(datatype_max)+1) << " bits), ";
                std::cout << "" << boost::format("%.0f") % (100.0*clip_i/num_rx_samps) << "% I clip.";
                std::cout << std::endl;

            }
        }
    }

    const auto actual_stop_time = std::chrono::steady_clock::now();

    stream_cmd.stream_mode = uhd::stream_cmd_t::STREAM_MODE_STOP_CONTINUOUS;
    rx_stream->issue_stream_cmd(stream_cmd);
    rx_stream.reset();

    // clean up transmit worker
    stop_signal_called = true;
    if (enable_tx) {
        transmit_thread.join();
        tx_stream.reset();
    }

    usrp.reset();

    if (stats) {
        std::cout << std::endl;
        const double actual_duration_seconds =
            std::chrono::duration<float>(actual_stop_time - start_time).count();

        std::cout << boost::format("Received %d samples in %f seconds") % num_total_samps
                         % actual_duration_seconds
                  << std::endl;
        const double rate = (double)num_total_samps / actual_duration_seconds;
        std::cout << (rate / 1e6) << " Msps" << std::endl;
    }

    // wait for ZMQ a bit
    std::this_thread::sleep_for(std::chrono::milliseconds(int64_t(1000 * setup_time)));

    // finished
    std::cout << std::endl << "Done!" << std::endl << std::endl;

    return EXIT_SUCCESS;
}
