//
// Copyright 2010-2011,2014 Ettus Research LLC
// Copyright 2018 Ettus Research, a National Instruments Company
//
// Copyright 2021/2022 by Thomas Telkamp and Tammo Jan Dijkema
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

#include <boost/thread/thread.hpp>

unsigned long long num_total_samps = 0;

namespace po = boost::program_options;

static bool stop_signal_called = false;
void sig_int_handler(int)
{
    stop_signal_called = true;
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

int UHD_SAFE_MAIN(int argc, char* argv[])
{
    // variables to be set by po
    std::string file, type, ant_list, subdev, ref, wirefmt, channel_list, gain_list, udp_forward;
    size_t total_num_samps, spb;
    uint16_t port;
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
        // ("file", po::value<std::string>(&file)->default_value("usrp_samples.dat"), "name of the file to write binary samples to")
        // ("type", po::value<std::string>(&type)->default_value("short"), "sample type: double, float, or short")
        // ("nsamps", po::value<size_t>(&total_num_samps)->default_value(0), "total number of samples to receive")
        // ("duration", po::value<double>(&total_time)->default_value(0), "total number of seconds to receive")
        // ("spb", po::value<size_t>(&spb)->default_value(10000), "samples per buffer")
        ("rate", po::value<double>(&rate)->default_value(1e6), "rate of incoming samples")
        ("freq", po::value<double>(&freq)->default_value(0.0), "RF center frequency in Hz")
        ("lo-offset", po::value<double>(&lo_offset)->default_value(0.0),
            "Offset for frontend LO in Hz (optional)")
        ("gain", po::value<std::string>(&gain_list), "gain(s) for the RF chain")
        ("ant", po::value<std::string>(&ant_list), "antenna selection")
        ("subdev", po::value<std::string>(&subdev), "subdevice specification")
        ("channel", po::value<std::string>(&channel_list)->default_value("0"), "which channel(s) to use (specify \"0\", \"1\", \"0,1\", etc)")
        ("bw", po::value<double>(&bw), "analog frontend filter bandwidth in Hz")
        ("ref", po::value<std::string>(&ref)->default_value("internal"), "reference source (internal, external, mimo, gpsdo)")
        ("wirefmt", po::value<std::string>(&wirefmt)->default_value("sc16"), "wire format (sc8, sc16 or s16)")
        ("setup", po::value<double>(&setup_time)->default_value(1.0), "seconds of setup time")
        ("udp", po::value<std::string>(&udp_forward), "DIFI UDP forward address")
        ("progress", "periodically display short-term bandwidth")
        ("stats", "show average bandwidth on exit")
        ("pps", "use external pps signal")
        ("vrt", "publish IQ using VRT over ZeroMQ (PUB on port 50100")
        // ("zmq", "publish IQ using ZeroMQ (PUB on port 50200 en 50201 (metadata)")
        ("int-second", "align start of reception to integer second")
        ("null", "run without streaming")
        ("continue", "don't abort on a bad packet")
        ("skip-lo", "skip checking LO lock status")
        ("int-n", "tune USRP with integer-N tuning")
        ("port", po::value<uint16_t>(&port)->default_value(50100), "DIFI ZMQ port")
        ("hwm", po::value<int>(&hwm)->default_value(10000), "DIFI ZMQ HWM")
    ;
    // clang-format on
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // print the help message
    if (vm.count("help")) {
        std::cout << boost::format("UHD RX samples to Vita49 %s") % desc << std::endl;
        std::cout << std::endl
                  << "This application streams data from a single channel of a USRP "
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

    // detect which channels to use
    std::vector<std::string> channel_strings;
    std::vector<size_t> channel_nums;
    boost::split(channel_strings, channel_list, boost::is_any_of("\"',"));
    for (size_t ch = 0; ch < channel_strings.size(); ch++) {
        size_t chan = std::stoi(channel_strings[ch]);
        if (chan >= usrp->get_rx_num_channels()) {
            throw std::runtime_error("Invalid channel(s) specified.");
        } else
            channel_nums.push_back(std::stoi(channel_strings[ch]));
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
        if (vm.count("freq")) { // with default of 0.0 this will always be true
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
            size_t gain = (gains.size() > channel) ? gains[ch] : gains[0];
            std::cout << boost::format("Setting RX Gain: %f dB...") % gain << std::endl;
            usrp->set_rx_gain(gain, channel);
            std::cout << boost::format("Actual RX Gain: %f dB...")
                             % usrp->get_rx_gain(channel)
                      << std::endl
                      << std::endl;
        }

        // set the IF filter bandwidth
        if (!vm.count("bw")) {
           bw = rate;
        }
        std::cout << boost::format("Setting RX Bandwidth: %f MHz...") % (bw / 1e6)
                  << std::endl;
        usrp->set_rx_bandwidth(bw, channel);
        std::cout << boost::format("Actual RX Bandwidth: %f MHz...")
                         % (usrp->get_rx_bandwidth(channel) / 1e6)
                  << std::endl
                  << std::endl;

        // set the antenna(s)
        if (antennas.size()) {
            std::string ant = (antennas.size() > channel) ? antennas[ch] : antennas[0];
            std::cout << boost::format("Setting Antenna: %s") % ant << std::endl;
            usrp->set_rx_antenna(ant, channel);
            std::cout << boost::format("Actual Antenna: %s")
                             % usrp->get_rx_antenna(channel)
                      << std::endl
                      << std::endl;
        }
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

    // Create a temporary rx_stream before we set the time.
    uhd::stream_args_t tmp_stream_args("sc16", "sc16");
    tmp_stream_args.channels             = channel_nums;
    uhd::rx_streamer::sptr tmp_rx_stream = usrp->get_rx_stream(tmp_stream_args);

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
        std::cout << boost::format("Wait for PPS sync...") << std::endl;
        gettimeofday(&time_now, nullptr);
        usrp->set_time_unknown_pps(uhd::time_spec_t(time_now.tv_sec + 2.0));
        boost::this_thread::sleep_for(boost::chrono::milliseconds(2100));
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

    if (total_num_samps == 0) {
        std::signal(SIGINT, &sig_int_handler);
        std::cout << "Press Ctrl + C to stop streaming..." << std::endl;
    }

    const std::string& cpu_format = "sc16";
    const std::string& wire_format = wirefmt;
    
    unsigned long long num_requested_samples = total_num_samps;
    double time_requested = total_time;
    bool int_second             = (bool)vm.count("int-second");

    // fixed buffer size
    size_t samps_per_buff = 10000; // spb

    #define SIZE 10007
    uint32_t buffer[SIZE];

    // create a receive streamer
    uhd::stream_args_t stream_args(cpu_format, wire_format);
    stream_args.channels             = channel_nums;
    uhd::rx_streamer::sptr rx_stream = usrp->get_rx_stream(stream_args);

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

    struct vrt_packet p;
    vrt_init_packet(&p);

    /* Warn if not standards compliant */
    if (vrt_is_platform_little_endian()) {
        printf("Warning: little endian support is work in progress.\n");
    }

    /* Configure */
    p.header.packet_type         = VRT_PT_IF_DATA_WITH_STREAM_ID;

    p.header.packet_size         = SIZE;
    p.header.tsm                 = VRT_TSM_FINE;
    p.header.tsi                 = VRT_TSI_OTHER; // unix time
    p.header.tsf                 = VRT_TSF_REAL_TIME;
    p.fields.stream_id           = (uint32_t)rand();
    p.words_body                 = 10000;

    p.header.has.class_id        = true;
    p.fields.class_id.oui        = 0x6A621E; // DIFI OUI
    p.fields.class_id.information_class_code = 0;
    p.fields.class_id.packet_class_code = 0;

    p.header.has.trailer         = false;

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

    // setup streaming
    uhd::stream_cmd_t stream_cmd((num_requested_samples == 0)
                                     ? uhd::stream_cmd_t::STREAM_MODE_START_CONTINUOUS
                                     : uhd::stream_cmd_t::STREAM_MODE_NUM_SAMPS_AND_DONE);
    stream_cmd.num_samps  = size_t(num_requested_samples);
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
    auto last_context                    = start_time;
    unsigned long long last_update_samps = 0;


    if (int_second) {
        stop_time += std::chrono::milliseconds(int64_t(1000.0*(double)(stream_cmd.time_spec.get_real_secs()-usrp->get_time_now().get_real_secs())));
    }

    // Run this loop until either time expired (if a duration was given), until
    // the requested number of samples were collected (if such a number was
    // given), or until Ctrl-C was pressed.

    uint32_t frame_count = 0;

    while (not stop_signal_called) {
           // and (num_requested_samples != num_total_samps or num_requested_samples == 0)
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
                           "Got an overflow indication. Please consider the following:\n"
                           "  Your write medium must sustain a rate of %fMB/s.\n"
                           "  Dropped samples will not be written to the file.\n")
                           % (usrp->get_rx_rate() * sizeof(std::complex<short>) / 1e6);
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
        }

        stream_cmd.stream_now = false;
   
        num_total_samps += num_rx_samps;

        p.body = (char*)buff_ptrs[0];
        p.header.packet_count = (uint8_t)frame_count%16;
        p.fields.integer_seconds_timestamp = md.time_spec.get_full_secs();
        p.fields.fractional_seconds_timestamp = (uint64_t)1e12 * md.time_spec.get_frac_secs();

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

            // VITA
            /* Initialize to reasonable values */
            struct vrt_packet pc;
            vrt_init_packet(&pc);

            /* Configure. Note that context packets cannot have a trailer word. */
            pc.header.packet_type         = VRT_PT_IF_CONTEXT;
            pc.fields.stream_id           = p.fields.stream_id;
            pc.header.has.class_id = true;
            pc.if_context.has.bandwidth   = true;
            pc.if_context.has.sample_rate = true;
            pc.if_context.has.reference_point_identifier = true;
            pc.if_context.has.if_reference_frequency = true;
            pc.if_context.has.rf_reference_frequency = true;
            pc.if_context.has.if_band_offset = true;
            pc.if_context.has.reference_level = true;
            pc.if_context.has.gain = true;
            pc.if_context.has.timestamp_adjustment = true;
            pc.if_context.has.timestamp_calibration_time = true;
            pc.if_context.has.state_and_event_indicators = true;
            pc.if_context.has.data_packet_payload_format = true;

            pc.if_context.data_packet_payload_format.packing_method = VRT_PM_LINK_EFFICIENT;
            pc.if_context.data_packet_payload_format.real_or_complex = VRT_ROC_COMPLEX_CARTESIAN;
            pc.if_context.data_packet_payload_format.data_item_format = VRT_DIF_SIGNED_FIXED_POINT;
            pc.if_context.data_packet_payload_format.sample_component_repeat = false;
            pc.if_context.data_packet_payload_format.item_packing_field_size = 31;
            pc.if_context.data_packet_payload_format.data_item_size = 15;

            pc.header.tsm                 = VRT_TSM_COARSE;
            pc.header.tsi                 = VRT_TSI_OTHER; // unix time (?)
            pc.header.tsf                 = VRT_TSF_REAL_TIME;

            pc.fields.integer_seconds_timestamp = md.time_spec.get_full_secs();
            pc.fields.fractional_seconds_timestamp = (uint64_t)1e12 * md.time_spec.get_frac_secs();

            pc.if_context.bandwidth                         = usrp->get_rx_bandwidth(); // 0.8*usrp->get_rx_rate(); // bandwith is set to 80% of sample rate
            pc.if_context.sample_rate                       = usrp->get_rx_rate();
            pc.if_context.rf_reference_frequency            = usrp->get_rx_freq();
            pc.if_context.if_reference_frequency            = 0; // Zero-IF
            pc.if_context.rf_reference_frequency_offset     = 0;
            pc.if_context.gain.stage1                       = usrp->get_rx_gain();
            pc.if_context.gain.stage2                       = 0;

            pc.if_context.state_and_event_indicators.has.reference_lock = true;
            pc.if_context.state_and_event_indicators.reference_lock = ((ref == "external") or (ref=="gpsdo"));

            pc.if_context.state_and_event_indicators.has.calibrated_time = true;
            pc.if_context.state_and_event_indicators.calibrated_time = ((vm.count("pps")) or (ref=="gpsdo"));

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
  
    // finished
    std::cout << std::endl << "Done!" << std::endl << std::endl;

    return EXIT_SUCCESS;
}
