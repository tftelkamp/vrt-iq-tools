//
// Copyright 2026 by Thomas Telkamp
//
// SPDX-License-Identifier: MIT
//

#include <unistd.h>

#include <chrono>
#include <complex>
#include <csignal>
#include <fstream>
#include <iostream>
#include <thread>

// IIO v0 (libiio-v0)
#include <iio.h>

// AD9361 v0 (libad9361-iio-v0)
#include <ad9361.h>

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
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/circular_buffer.hpp>
#include <boost/algorithm/string.hpp>


// VRT tools functions
#include "vrt-tools.h"

unsigned long long num_total_samps = 0;

namespace po = boost::program_options;

static bool stop_signal_called = false;
void sig_int_handler(int)
{
    stop_signal_called = true;
}

// IIO

#define IIO_ENSURE(expr) { \
    if (!(expr)) { \
        (void) fprintf(stderr, "assertion failed (%s:%d)\n", __FILE__, __LINE__); \
        (void) abort(); \
    } \
}

#define BLOCK_SIZE (4*VRT_SAMPLES_PER_PACKET)

/* RX is input, TX is output */
enum iodev { RX, TX };

/* common RX and TX streaming params */
struct stream_cfg {
    long long bw_hz; // Analog banwidth in Hz
    long long fs_hz; // Baseband sample rate in Hz
    long long lo_hz; // Local oscillator frequency in Hz
    uint8_t gain;// Gain
    const char* rfport; // Port name
};

/* static scratch mem for strings */
static char tmpstr[64];

/* IIO structs required for streaming */
static struct iio_context *ctx   = NULL;
static struct iio_buffer  *rxbuf = NULL;

static struct iio_channel *rx_i[2] = { NULL };
static struct iio_channel *rx_q[2] = { NULL };


/* cleanup and exit */
static void shutdown(void)
{
    printf("* Destroying buffers\n");
    if (rxbuf) { iio_buffer_destroy(rxbuf); }

    printf("* Disabling streaming channels\n");
    if (rx_i[0]) { iio_channel_disable(rx_i[0]); }
    if (rx_q[0]) { iio_channel_disable(rx_q[0]); }
    if (rx_i[1]) { iio_channel_disable(rx_i[1]); }
    if (rx_q[1]) { iio_channel_disable(rx_q[1]); }

    printf("* Destroying context\n");
    if (ctx) { iio_context_destroy(ctx); }
}

/* check return value of attr_write function */
static void errchk(int v, const char* what) {
     if (v < 0) { fprintf(stderr, "Error %d writing to channel \"%s\"\nvalue may not be supported.\n", v, what); shutdown(); }
}

/* write attribute: long long int */
static void wr_ch_lli(struct iio_channel *chn, const char* what, long long val)
{
   errchk(iio_channel_attr_write_longlong(chn, what, val), what);
}

/* write attribute: string */
static void wr_ch_str(struct iio_channel *chn, const char* what, const char* str)
{
    errchk(iio_channel_attr_write(chn, what, str), what);
}

/* helper function generating channel names */
static char* get_ch_name(const char* type, int id)
{
    snprintf(tmpstr, sizeof(tmpstr), "%s%d", type, id);
    return tmpstr;
}

/* returns ad9361 phy device */
static struct iio_device* get_ad9361_phy(void)
{
    struct iio_device *dev =  iio_context_find_device(ctx, "ad9361-phy");
    IIO_ENSURE(dev && "No ad9361-phy found");
    return dev;
}

/* finds AD9361 streaming IIO devices */
static bool get_ad9361_stream_dev(enum iodev d, struct iio_device **dev)
{
    switch (d) {
    case TX: *dev = iio_context_find_device(ctx, "cf-ad9361-dds-core-lpc"); return *dev != NULL;
    case RX: *dev = iio_context_find_device(ctx, "cf-ad9361-lpc");  return *dev != NULL;
    default: IIO_ENSURE(0); return false;
    }
}

/* finds AD9361 streaming IIO channels */
static bool get_ad9361_stream_ch(enum iodev d, struct iio_device *dev, int chid, struct iio_channel **chn)
{
    *chn = iio_device_find_channel(dev, get_ch_name("voltage", chid), d == TX);
    if (!*chn)
        *chn = iio_device_find_channel(dev, get_ch_name("altvoltage", chid), d == TX);
    return *chn != NULL;
}

/* finds AD9361 phy IIO configuration channel with id chid */
static bool get_phy_chan(enum iodev d, int chid, struct iio_channel **chn)
{
    switch (d) {
    case RX: *chn = iio_device_find_channel(get_ad9361_phy(), get_ch_name("voltage", chid), false); return *chn != NULL;
    case TX: *chn = iio_device_find_channel(get_ad9361_phy(), get_ch_name("voltage", chid), true);  return *chn != NULL;
    default: IIO_ENSURE(0); return false;
    }
}

/* finds AD9361 local oscillator IIO configuration channels */
static bool get_lo_chan(enum iodev d, struct iio_channel **chn)
{
    switch (d) {
     // LO chan is always output, i.e. true
    case RX: *chn = iio_device_find_channel(get_ad9361_phy(), get_ch_name("altvoltage", 0), true); return *chn != NULL;
    case TX: *chn = iio_device_find_channel(get_ad9361_phy(), get_ch_name("altvoltage", 1), true); return *chn != NULL;
    default: IIO_ENSURE(0); return false;
    }
}

/* applies streaming configuration through IIO */
bool cfg_ad9361_streaming_ch(struct stream_cfg *cfg, enum iodev type, int chid)
{
    const struct iio_attr *attr;
    struct iio_channel *chn = NULL;

    // Configure phy and lo channels
    printf("* Acquiring AD9361 phy channel %d\n", chid);
    if (!get_phy_chan(type, chid, &chn)) {  return false; }

    wr_ch_str(chn, "rf_port_select",     cfg->rfport);
    wr_ch_lli(chn, "rf_bandwidth",       cfg->bw_hz);
    wr_ch_lli(chn, "sampling_frequency", cfg->fs_hz);
    wr_ch_str(chn, "gain_control_mode", "manual");
    wr_ch_lli(chn, "sampling_frequency", cfg->fs_hz);
    wr_ch_lli(chn, "hardwaregain", cfg->gain);

    // Configure LO channel
    printf("* Acquiring AD9361 %s lo channel\n", type == TX ? "TX" : "RX");
    if (!get_lo_chan(type, &chn)) { return false; }
    wr_ch_lli(chn, "frequency", cfg->lo_hz);

    return true;
}

int main(int argc, char* argv[])
{
    // variables to be set by po
    std::string file, type, ant_list, ref, channel_list, gain_list, merge_address;
    size_t total_num_samps;
    uint16_t instance, port, merge_port;
    int hwm;
    uint32_t stream_id;
    double rate, bw, freq, total_time, setup_time, lo_offset, if_freq, pps_offset;

    bool context_changed = true;
    bool merge;

    std::string uri;

    // setup the program options
    po::options_description desc("Allowed options");
    // clang-format off
    desc.add_options()
        ("help", "help message")
        ("uri", po::value<std::string>(&uri)->default_value("usb:"), "IIO uri")
        ("nsamps", po::value<size_t>(&total_num_samps)->default_value(0), "total number of samples to receive")
        ("duration", po::value<double>(&total_time)->default_value(0), "total number of seconds to receive")
        ("rate", po::value<double>(&rate)->default_value(1e6), "rate of incoming samples")
        ("freq", po::value<double>(&freq)->default_value(0.0), "RF center frequency in Hz")
        ("if-freq", po::value<double>(&if_freq)->default_value(0.0), "IF center frequency in Hz")
        ("gain", po::value<std::string>(&gain_list), "gain(s) for the RF chain")
        // ("ant", po::value<std::string>(&ant_list), "antenna selection")
        ("iio-channel", po::value<std::string>(&channel_list)->default_value("0"), "which usrp channel(s) to use (specify \"0\", \"1\", \"0,1\", etc)")
        ("bw", po::value<double>(&bw), "analog frontend filter bandwidth in Hz")
        ("ref", po::value<std::string>(&ref)->default_value("internal"), "reference source (internal, external) for context only")
        ("setup", po::value<double>(&setup_time)->default_value(1.0), "seconds of setup time")
        ("progress", "periodically display short-term bandwidth")
        ("stats", "show average bandwidth on exit")
        ("temp", "read temperature sensor")
        ("int-second", "align start of reception to integer second")
        ("null", "run without streaming")
        ("continue", "don't abort on a bad packet")
        ("port", po::value<uint16_t>(&port)->default_value(50100), "VRT ZMQ port")
        ("instance", po::value<uint16_t>(&instance)->default_value(0), "VRT ZMQ instance")
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
        std::cout << boost::format("IIO AD936x samples to Vita49 over ZMQ. %s") % desc << std::endl;
        std::cout << std::endl
                  << "This application streams data from an IIO AD936x "
                     "device to Vita49 over ZMQ.\n"
                  << std::endl;
        return ~0;
    }

    bool bw_summary             = vm.count("progress") > 0;
    bool stats                  = vm.count("stats") > 0;
    bool null                   = vm.count("null") > 0;
    bool continue_on_bad_packet = vm.count("continue") > 0; 
    bool enable_temp            = vm.count("temp") > 0;
    bool int_second             = vm.count("int-second");

    struct vrt_packet p;
    vrt_init_packet(&p);

    /* Warn if not standards compliant */
    if (vrt_is_platform_little_endian()) {
        printf("Warning: little endian support is work in progress.\n");
    }

    /* VRT init */
    vrt_init_data_packet(&p);
    
    // Only 1 channel
    p.fields.stream_id = 1;

    // detect channels to use
    std::vector<std::string> channel_strings;
    std::vector<size_t> channel_nums;
    boost::split(channel_strings, channel_list, boost::is_any_of("\"',"));

    std::vector<size_t> gains;
    if (vm.count("gain")) {
        std::vector<std::string> gain_strings;
        boost::split(gain_strings, gain_list, boost::is_any_of("\"',"));
        for (size_t ch = 0; ch < gain_strings.size(); ch++) {
            gains.push_back(std::stoi(gain_strings[ch]));
        }
    }
    
    // ZMQ
    void *zmq_server;
    void *zmq_control;

    uint16_t main_port;

    if (vm.count("port") > 0) {
        main_port = port;
    } else {
        main_port = DEFAULT_MAIN_PORT + MAX_CHANNELS*instance;
    }

    void *context = zmq_ctx_new();
    void *responder = zmq_socket(context, ZMQ_PUB);
    int rc = zmq_setsockopt (responder, ZMQ_SNDHWM, &hwm, sizeof hwm);
    assert(rc == 0);
    std::string connect_string = "tcp://*:" + std::to_string(main_port);
    rc = zmq_bind(responder, connect_string.c_str());
    assert (rc == 0);
    zmq_server = responder;

    responder = zmq_socket(context, ZMQ_SUB);
    rc = zmq_bind(responder, "tcp://*:50300");
    assert (rc == 0);
    zmq_control = responder;
    zmq_setsockopt(zmq_control, ZMQ_SUBSCRIBE, "", 0);

    // Merge
    void *merge_zmq = zmq_socket(context, ZMQ_SUB);
    if (merge) {
        connect_string = "tcp://" + merge_address + ":" + std::to_string(merge_port);
        rc = zmq_connect(merge_zmq, connect_string.c_str());
        assert(rc == 0);
        zmq_setsockopt(merge_zmq, ZMQ_SUBSCRIBE, "", 0);
    }

    // detect which channels to use
    for (size_t ch = 0; ch < channel_strings.size(); ch++) {
        size_t chan = std::stoi(channel_strings[ch]);
        if (chan >= 2) { // max 2 channels
            throw std::runtime_error("Invalid channel(s) specified.");
        } else
            channel_nums.push_back(std::stoi(channel_strings[ch]));
    }
    
    // INIT IIO

    uint32_t sample_rate = rate;

    // Streaming devices
    struct iio_device *rx;

    // RX sample counters
    size_t nrx = 0;

    // RX sample size
    size_t rx_sample_sz;

    // Stream configurations
    struct stream_cfg rxcfg;

    int err;

    uint32_t bandwidth;
    if (vm.count("bw") > 0)
        bandwidth = bw;
    else
        bandwidth = 0.85 * rate;
  
    printf("* Acquiring IIO context\n");
    
    IIO_ENSURE((ctx = iio_create_context_from_uri(uri.c_str())) && "No context");
    IIO_ENSURE(iio_context_get_devices_count(ctx) > 0 && "No devices");


    printf("* Acquiring AD9361 streaming devices\n");
    IIO_ENSURE(get_ad9361_stream_dev(RX, &rx) && "No rx dev found");

    printf("* Configuring AD9361 for streaming\n");

    ad9361_set_bb_rate(get_ad9361_phy(),(unsigned long)sample_rate);

    for (size_t ch = 0; ch < channel_nums.size(); ch++) {

        size_t channel = channel_nums[ch];

        std::cout << "Configuring RX Channel " << channel << std::endl;

        size_t gain = (gains.size() > ch) ? gains[ch] : gains[0];

        // RX stream config
        rxcfg.bw_hz  = (uint64_t)bandwidth;
        rxcfg.fs_hz  = (uint64_t)sample_rate;
        rxcfg.lo_hz  = (uint64_t)freq;
        rxcfg.gain   = (uint8_t)gain;
        rxcfg.rfport = "A_BALANCED"; // port A (locked on Pluto)
  
        IIO_ENSURE(cfg_ad9361_streaming_ch(&rxcfg, RX, channel) && "RX port not found");

        printf("* Initializing AD9361 IIO streaming channels\n");
        IIO_ENSURE(get_ad9361_stream_ch(RX, rx, 2*channel, &rx_i[ch]) && "RX chan i not found");
        IIO_ENSURE(get_ad9361_stream_ch(RX, rx, 2*channel+1, &rx_q[ch]) && "RX chan q not found");
    }
   
    printf("* Enabling IIO streaming channels\n");
    for (size_t ch = 0; ch < channel_nums.size(); ch++) {
        iio_channel_enable(rx_i[ch]);
        iio_channel_enable(rx_q[ch]);
    }   

    printf("* Creating non-cyclic IIO buffers\n");
    rxbuf = iio_device_create_buffer(rx, BLOCK_SIZE, false);
    if (!rxbuf) {
        perror("Could not create RX buffer");
        shutdown();
    }

    rx_sample_sz = iio_device_get_sample_size(rx);

    // END IIO

    struct timeval time_now{};
    gettimeofday(&time_now, nullptr);

    // if (vm.count("pps")) {
    // }

    gettimeofday(&time_now, nullptr);
    std::cout << boost::format("PC Clock time: %.6f seconds\n") % (time_now.tv_sec + (double)time_now.tv_usec / 1e6);

    // RX

    if (total_num_samps == 0) {
        std::signal(SIGINT, &sig_int_handler);
        std::cout << "Press Ctrl + C to stop streaming..." << std::endl;
    }
    
    unsigned long long num_requested_samples = total_num_samps;

    uint32_t buffer[VRT_DATA_PACKET_SIZE];

    if (total_time > 0)  
        num_requested_samples = total_time * rate;

    // Run this loop until either time expired (if a duration was given), until
    // the requested number of samples were collected (if such a number was
    // given), or until Ctrl-C was pressed.

    uint64_t frame_count = 0;
    uint64_t sys_timestamp = 0;

    // Create a circular buffer with a capacity for 4 * IIO buffer (*2 for IQ and *2 for 2 channels)

    const ssize_t cb_size = 4*BLOCK_SIZE*2*2;

    typedef boost::circular_buffer<int16_t> iq_cb;
    std::vector<iq_cb> per_channel_cb;

    for (size_t ch = 0; ch < channel_nums.size(); ch++) {
        per_channel_cb.emplace_back(cb_size);
    }

    bool first_frame = true;
    size_t num_rx_samps;

    int16_t bodydata[2][VRT_SAMPLES_PER_PACKET*2];

    const struct iio_block *rxblock;

    int64_t time_per_block = (1e12*VRT_SAMPLES_PER_PACKET)/sample_rate;
    
    if (int_second) {
        gettimeofday(&time_now, nullptr);
        struct timeval new_time{};
        gettimeofday(&new_time, nullptr);
        while (time_now.tv_sec==new_time.tv_sec)
            gettimeofday(&new_time, nullptr);
    }
    
    // flush merge queue
    if (merge)
        while ( zmq_recv(merge_zmq, buffer, 100000, ZMQ_NOBLOCK) > 0 ) { }

    // time keeping
    auto start_time = std::chrono::steady_clock::now();

    // Track time and samps between updating the BW summary
    auto last_update                     = start_time;
    auto last_context                    = start_time - std::chrono::milliseconds(2*VRT_CONTEXT_INTERVAL);
    unsigned long long last_update_samps = 0;

    while (not stop_signal_called 
           and (num_requested_samples > num_total_samps or num_requested_samples == 0)) {

        ssize_t nbytes_rx;
        int16_t *p_dat, *p_end;
        ptrdiff_t p_inc;

        const auto now = std::chrono::steady_clock::now();

        gettimeofday(&time_now, nullptr);

        // Refill RX buffer
        nbytes_rx = iio_buffer_refill(rxbuf);
        if (nbytes_rx < 0) { printf("Error refilling buf %d\n",(int) nbytes_rx); shutdown(); }

        num_rx_samps = 0;

        /* READ: Get pointers to RX buf and read IQ from RX buf port 0 */
        p_inc = iio_buffer_step(rxbuf);
        p_end = (int16_t*)iio_buffer_end(rxbuf);

        for (size_t ch = 0; ch < channel_nums.size(); ch++) {
            size_t channel = channel_nums[ch];

            for (p_dat = (int16_t*)iio_buffer_first(rxbuf, rx_i[ch]); p_dat < p_end;
                 p_dat += p_inc / sizeof(*p_dat)) {
                const int16_t i = ((int16_t*)p_dat)[0]; // Real (I)
                const int16_t q = ((int16_t*)p_dat)[1]; // Imag (Q)

                per_channel_cb[ch].push_back( (int16_t)(i) );
                per_channel_cb[ch].push_back( (int16_t)(q) );
            }
        }

        num_rx_samps = nbytes_rx / p_inc;

        if (first_frame) {
            std::cout << boost::format(
                             "First frame: %u samples, %u full secs, %.09f frac secs")
                             % (num_rx_samps) % (time_now.tv_sec)
                             % ((double)(time_now.tv_usec/1e6))
                      << std::endl;
            first_frame = false;
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

            // Use host time, not very accurate
            pc.fields.integer_seconds_timestamp = time_now.tv_sec;
            pc.fields.fractional_seconds_timestamp = 1e6*time_now.tv_usec;

            for (size_t ch = 0; ch < channel_nums.size(); ch++) {
                size_t channel = channel_nums[ch];

                size_t gain = (gains.size() > ch) ? gains[ch] : gains[0];

                if (enable_temp) {
                    double temp;
                    iio_channel_attr_read_double(
                        iio_device_find_channel( get_ad9361_phy(), get_ch_name("temp", 0), RX ),
                        "input",
                        &temp );
                    pc.if_context.has.temperature = true;
                    pc.if_context.temperature = (float)temp/1000;
                }

                if (context_changed)
                    pc.if_context.context_field_change_indicator = true;
                else
                    pc.if_context.context_field_change_indicator = false;
       
                pc.fields.stream_id = 1<<ch;
                pc.if_context.bandwidth                         = rxcfg.bw_hz;
                pc.if_context.sample_rate                       = rxcfg.fs_hz;
                pc.if_context.rf_reference_frequency            = rxcfg.lo_hz;
                pc.if_context.rf_reference_frequency_offset     = 0;
                pc.if_context.if_reference_frequency            = if_freq; // 0 for Zero-IF
                pc.if_context.if_band_offset                    = 0;
                pc.if_context.gain.stage1                       = gain;
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
                zmq_send (zmq_server, buffer, rv*4, 0);

            }
            
            context_changed = false;
        }

        ssize_t block_count = 0;

        while (per_channel_cb[0].size() > 2*VRT_SAMPLES_PER_PACKET ) {

            for (size_t ch = 0; ch < channel_nums.size(); ch++) {
                size_t channel = channel_nums[ch];
        
                for (uint32_t i = 0; i < 2*VRT_SAMPLES_PER_PACKET; i++) {
                    bodydata[ch][i] = (int16_t)per_channel_cb[ch].front();
                    per_channel_cb[ch].pop_front();
                }   
       
                num_total_samps += VRT_SAMPLES_PER_PACKET;

                p.fields.stream_id = 1<<ch;
                p.body = bodydata[ch];
                p.header.packet_count = (uint8_t)frame_count%16;

                // Use host time, not very accurate
                p.fields.integer_seconds_timestamp = time_now.tv_sec;
                p.fields.fractional_seconds_timestamp = 1e6*time_now.tv_usec + block_count*time_per_block;

                if (p.fields.fractional_seconds_timestamp > 1e12) {
                    p.fields.fractional_seconds_timestamp -= 1e12;
                    p.fields.integer_seconds_timestamp += 1;      
                }
        
                zmq_msg_t msg;
                int rc = zmq_msg_init_size (&msg, VRT_DATA_PACKET_SIZE*4);

                int32_t rv = vrt_write_packet(&p, zmq_msg_data(&msg), VRT_DATA_PACKET_SIZE, true);

                zmq_msg_send(&msg, zmq_server, 0);

                zmq_msg_close(&msg);

            }

            block_count++;
            frame_count++;
        }

        // Merge
        if (merge) {
            int mergelen;
            while ( (mergelen = zmq_recv(merge_zmq, buffer, 100000, ZMQ_NOBLOCK)) > 0  ) {
                zmq_msg_t msg;
                zmq_msg_init_size (&msg, mergelen);
                memcpy (zmq_msg_data(&msg), buffer, mergelen);
                zmq_msg_send(&msg, zmq_server, 0);
                zmq_msg_close(&msg);
            }
        }

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

                if (c.has.rf_reference_frequency) {
                    freq = (double)round(c.rf_reference_frequency);

                    // TODO set freq
                    // std::cout << boost::format("    Setting RX Freq: %f MHz...") % (freq / 1e6)
                    //           << std::endl;
                }
                if (c.has.gain) {
                    uint8_t rx_gain = c.gain.stage1;

                    // TODO set gain
                    // printf("Info: set gain index to %" PRIu8 "\n", rx_gain);
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

                for (size_t ch = 0; ch < channel_nums.size(); ch++) {

                    double max_iq = 0;
                    uint32_t clip_iq = 0;

                    double datatype_max = 32767.;

                    for (int i=0; i < VRT_SAMPLES_PER_PACKET; i++ ) {
                        std::complex<int16_t> sample (bodydata[ch][2*i], bodydata[ch][2*i+1]);
                        max_iq = fmax(max_iq, fmax(fabs(sample.real()), fabs(sample.imag())));
                        if (fabs(sample.real()) > datatype_max*0.99 || fabs(sample.imag()) > datatype_max*0.99)
                            clip_iq++;
                    }

                    std::cout << "CH" << boost::format("%u") % ch << ": ";
                    std::cout << boost::format("%3.0f") % (20*log10(max_iq/datatype_max)) << " dBFS (";
                    std::cout << boost::format("%2.0f") % ceil(log2(max_iq)+1) << "/";
                    std::cout << (int)ceil(log2(datatype_max)+1) << " bits), ";
                    std::cout << "" << boost::format("%2.0f") % (100.0*clip_iq/num_rx_samps) << "% clip. ";
                }

                std::cout << std::endl;

            }
        }


    }

    const auto actual_stop_time = std::chrono::steady_clock::now();

    // clean up transmit worker
    stop_signal_called = true;

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

    shutdown();

    // wait for ZMQ a bit
    std::this_thread::sleep_for(std::chrono::milliseconds(int64_t(1000 * setup_time)));

    zmq_close(zmq_server);
    zmq_close(zmq_control);
    zmq_close(merge_zmq);
    zmq_close(responder);
  
    // finished
    std::cout << std::endl << "Done!" << std::endl << std::endl;

    return EXIT_SUCCESS;
}
