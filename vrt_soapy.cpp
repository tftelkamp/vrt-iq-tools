#include <SoapySDR/Device.hpp>
#include <SoapySDR/Registry.hpp>
#include <SoapySDR/Formats.hpp>
#include <iostream>

#include <zmq.h>
#include "vrt-tools.h"
#include <memory>
#include <utility>
#include <vector>

#define VRT_SOAPY_DEBUG

#ifdef VRT_SOAPY_DEBUG

#define vrt_soapy_debug(message)        \
do {                                    \
    std::cout << message << std::endl;  \
} while(0);                             \

#else

#define vrt_soapy_debug(message)

#endif

struct VrtStream {
    VrtStream(void *subscriber, std::string format) :
    zmq_subscriber(subscriber), format(std::move(format)), cur_packet_idx(0)
    {
        init_context(&vrt_context);
    }
    void *zmq_subscriber;
    std::string format;
    context_type vrt_context;
    uint32_t zmq_buffer[ZMQ_BUFFER_SIZE];
    size_t cur_packet_idx;
};

void* create_zmq_subscriber(void* zmq_ctx, int port, const int hwm = 10000)
{
    vrt_soapy_debug(std::string("Create subscriber"));
    void* subscriber = zmq_socket(zmq_ctx, ZMQ_SUB);
    int rc = zmq_setsockopt (subscriber, ZMQ_RCVHWM, &hwm, sizeof hwm);
    if (rc != 0) {
        throw std::runtime_error("VrtDevice failed to set RCVHWM on ZMQ socket");
    }
    zmq_setsockopt(subscriber, ZMQ_SUBSCRIBE, "", 0);

    std::string addr = "tcp://localhost:" + std::to_string(port);
    zmq_connect(subscriber, addr.c_str());
    return subscriber;
};

template<typename T>
void iterate_packet_buffer(
    VrtStream *stream,
    void * const *buffer,
    const packet_type &packet,
    bool &start_rx,
    unsigned long &cur_buff_idx,
    const size_t numElems)
{
    int16_t re;
    int16_t img;
    auto *buff0 = static_cast<std::complex<T>*>(buffer[0]);
    for (; stream->cur_packet_idx < packet.num_rx_samps; stream->cur_packet_idx++) {
        memcpy(&re, (char*)&stream->zmq_buffer[packet.offset+stream->cur_packet_idx], sizeof(int16_t));
        memcpy(&img, (char*)&stream->zmq_buffer[packet.offset+stream->cur_packet_idx]+sizeof(int16_t), sizeof(int16_t));
        buff0[cur_buff_idx] = std::complex<T>(re, img);
        cur_buff_idx++;

        // reached max elements to return in this call
        if (cur_buff_idx >= numElems-1) {
            start_rx = false;
            return;
        }
    }
}

class VrtDevice : public SoapySDR::Device {
    public:

    VrtDevice(int vrt_instance, float freq, float rate): SoapySDR::Device(),
        vrt_instance_(vrt_instance), freq_(freq), rate_(rate)
    {
        vrt_soapy_debug(std::string("VrtDevice create"));
        zmq_ctx_ = zmq_ctx_new();
        port_ = DEFAULT_MAIN_PORT + MAX_CHANNELS * vrt_instance_;
    };

    ~VrtDevice() override
    {
        vrt_soapy_debug(std::string("VrtDevice destroy"));
        zmq_ctx_destroy(zmq_ctx_);
    }

    size_t getNumChannels(const int dir) const override
    {
        vrt_soapy_debug(std::string("getNumChannels"));
        return (dir == SOAPY_SDR_RX) ? 1 : 0;
    }

    std::vector<std::string> getStreamFormats(int, size_t) const override
    {
        vrt_soapy_debug(std::string("getStreamFormats"));
        return {SOAPY_SDR_CF32, SOAPY_SDR_CS16};
    }

    std::string getNativeStreamFormat(const int direction, const size_t channel, double &fullScale) const override
    {
        vrt_soapy_debug(std::string("getNativeStreamFormat"));
        //check that direction is SOAPY_SDR_RX
        if (direction != SOAPY_SDR_RX) {
            throw std::runtime_error("VrtDevice is RX only, use SOAPY_SDR_RX");
        }

        // fullScale = INT16_MAX;
        return SOAPY_SDR_CS16;
    }

    double getSampleRate(int, size_t) const override
    {
        vrt_soapy_debug(std::string("getSampleRate"));
        return rate_;
    }

    std::string getDriverKey() const override {
        vrt_soapy_debug(std::string("getDriverKey"));
        return "vrt_device";
    }

    std::string getHardwareKey() const override {
        vrt_soapy_debug(std::string("getHardwareKey"));
        return "VRTZMQ";
    }


    std::vector<std::string> listFrequencies(
        const int direction, const size_t channel) const override
    {
        vrt_soapy_debug(std::string("listFrequencies"));
        return {"RF"};
    }

    SoapySDR::RangeList getFrequencyRange(int, size_t) const override
    {
        vrt_soapy_debug(std::string("getFrequencyRange"));
        return { SoapySDR::Range(freq_, freq_) };
    }

    void setSampleRate(int, size_t, double) override
    {
        vrt_soapy_debug(std::string("setSampleRate"));
    }

    size_t getStreamMTU(SoapySDR::Stream *s) const
    {
        return VRT_SAMPLES_PER_PACKET;
    }

    SoapySDR::Stream *setupStream(
        const int direction,
        const std::string &format,
        const std::vector<size_t> &,
        const SoapySDR::Kwargs &args
    ) override
    {
        vrt_soapy_debug(std::string("setupStream"));
        auto formats = getStreamFormats(0, 0);
        if (std::find(formats.begin(), formats.end(), format) == formats.end()) {
            throw std::runtime_error(std::string("VrtDevice does not support format: ") + format);
        }
        void* subscriber = create_zmq_subscriber(zmq_ctx_, port_);
        return reinterpret_cast<SoapySDR::Stream *>(
            new VrtStream(subscriber, format)
        );
    }

    void closeStream(SoapySDR::Stream *stream) override
    {
        vrt_soapy_debug(std::string("closeStream"));
        delete reinterpret_cast<VrtStream *>(stream);
    }

    int activateStream(
        SoapySDR::Stream *, int, long long, size_t) override
    {
        vrt_soapy_debug(std::string("activateStream"));
        return 0;
    }

    int deactivateStream(
        SoapySDR::Stream *, int, long long) override
    {
        vrt_soapy_debug(std::string("deactivateStream"));
        return 0;
    }

    int readStream(
        SoapySDR::Stream *s, void * const *buffs, size_t numElems,
        int &flags, long long &timeNs, long timeoutUs
    ) override
    {
        flags = 0;
        timeNs = 0;
        bool start_rx = true;
        size_t cur_buff_idx = 0;
        long timeoutRemaining = timeoutUs;
        auto deadline = zmq_stopwatch_start();
        VrtStream *stream = reinterpret_cast<VrtStream*>(s);

        while (start_rx) {
            // TODO(Corne): Support multiple / selectable channels
            packet_type cur_packet;
            cur_packet.channel_filt = 1;

            // Detect timeout and return with current number of elements
            timeoutRemaining -= zmq_stopwatch_intermediate(deadline);
            if (timeoutRemaining <= 0) {
                // vrt_soapy_debug(std::string("readStream: timeout"));
                break;
            }

            // Detect if end of current packet reached and reset
            if (stream->cur_packet_idx+1 >= cur_packet.num_rx_samps) {
                stream->cur_packet_idx = 0;
            }

            // Receive next packet, noblock for timeout
            if (stream->cur_packet_idx == 0) {
                if (zmq_recv(
                    stream->zmq_subscriber, stream->zmq_buffer,
                    ZMQ_BUFFER_SIZE, ZMQ_NOBLOCK
                ) < 0) {
                    continue;
                }
            }

            if (!vrt_process(
                stream->zmq_buffer,
                sizeof(stream->zmq_buffer),
                &stream->vrt_context, &cur_packet)
            ) {
                std::cerr << "VrtDevice received and invalid Vita49 packet" << std::endl;
                continue;
            }

            if (not cur_packet.data) {
                continue;
            }

            // Go through the packets from the current sample (cur_packet_idx) until numElems reached or end of packet
            if (stream->format == SOAPY_SDR_CF32) {
                iterate_packet_buffer<float>(stream, buffs, cur_packet, start_rx, cur_buff_idx, numElems);
            }
            else if (stream->format == SOAPY_SDR_CS16) {
                iterate_packet_buffer<int16_t>(stream, buffs, cur_packet, start_rx, cur_buff_idx, numElems);
            }
        }

        zmq_stopwatch_stop(deadline);
        // vrt_soapy_debug(std::string("readStream: ") << ((cur_buff_idx == 0) ? 0 : cur_buff_idx+1));
        // if (stream->cur_packet_idx-1 < cur_packet.num_rx_samps) flags |= SOAPY_SDR_MORE_FRAGMENTS;
        return static_cast<int>(cur_buff_idx);
    }

    bool hasHardwareTime(const std::string &) const override
    {
        vrt_soapy_debug(std::string("hasHardwareTime"));
        return false;
    }

    std::vector<std::string> listAntennas(int, size_t) const override
    {
        vrt_soapy_debug(std::string("listAntennas"));
        return {};
    }

    void setAntenna(int, size_t, const std::string &) override
    {
        vrt_soapy_debug(std::string("setAntenna"));
    }

    std::string getAntenna(int, size_t) const override
    {
        vrt_soapy_debug(std::string("getAntenna"));
        return "";
    }

    SoapySDR::ArgInfoList getSettingInfo() const override
    {
        vrt_soapy_debug(std::string("getSettingInfo"));
        return {};
    }

    double getFrequency(int, size_t) const override
    {
        vrt_soapy_debug(std::string("getFrequency"));
        return freq_;
    }

    void setFrequency(
        int dir,
        size_t channel,
        double frequency,
        const SoapySDR::Kwargs &args
    ) override
    {
        vrt_soapy_debug(std::string("setFrequency"));
    }

    std::vector<double> listSampleRates(int dir, size_t channel) const override
    {
        vrt_soapy_debug(std::string("listSampleRates"));
        return {rate_};
    }
protected:
    int vrt_instance_;
    float freq_;
    float rate_;
    void *zmq_ctx_;
    int port_;
};



SoapySDR::KwargsList findVrtDevice(const SoapySDR::Kwargs &args) {
    std::vector<SoapySDR::Kwargs> instances;
    SoapySDR::Kwargs dev;

    std::vector<int> ports;
    for (size_t instance=0; instance<10; instance++) {
        ports.push_back(DEFAULT_MAIN_PORT + MAX_CHANNELS * instance);
    }

    const int TIMEOUT_MS = 200;  // Expect at least 10 packets per second
    const int hwm = 10000;

    void* ctx = zmq_ctx_new();

    std::vector<void*> subscribers;
    for (int port : ports) {
        subscribers.push_back(create_zmq_subscriber(ctx, port, hwm));
    }

    std::vector<zmq_pollitem_t> items(subscribers.size());
    for (size_t i = 0; i < subscribers.size(); i++) {
        items[i] = {subscribers[i], 0, ZMQ_POLLIN, 0};
    }

    std::vector<bool> found(ports.size(), false);
    int remaining = ports.size();

    auto deadline = zmq_stopwatch_start();

    uint32_t buffer[ZMQ_BUFFER_SIZE];
    context_type vrt_context;
    packet_type vrt_packet;
    // TODO(Corne): Support multiple / selectable channels
    vrt_packet.channel_filt = 1;
    init_context(&vrt_context);

    while (remaining > 0) {
        long elapsed_us = zmq_stopwatch_intermediate(deadline);
        long remaining_ms = TIMEOUT_MS - (elapsed_us / 1000);
        if (remaining_ms <= 0) break;
        int rc = zmq_poll(items.data(), items.size(), remaining_ms);
        if (rc <= 0) break;
        for (size_t i = 0; i < items.size(); i++) {
            if (items[i].events && (items[i].revents & ZMQ_POLLIN)) {
                remaining--;
                items[i].events = 0;

                if (args.count("vrt_instance") and args.at("vrt_instance") != std::to_string(i)) {
                    continue;
                }

                while (true) {
                    int len = zmq_recv(subscribers[i], buffer, ZMQ_BUFFER_SIZE, 0);
                    if (len <= 0) {
                        continue;
                    }
                    if (not vrt_process(buffer, sizeof(buffer), &vrt_context, &vrt_packet)) {
                        std::cerr << "Not a Vita49 packet?" << std::endl;
                        break;
                    }
                    if (vrt_packet.context) {
                        float rate = vrt_context.sample_rate;
                        float freq = vrt_context.rf_freq;
                        dev["freq"] = std::to_string(freq);
                        dev["driver"] = "vrt_device";

                        std::string rate_str = (boost::format("%.1f Msps") % (rate / 1e6)).str();
                        if (rate < 1e6) {
                            std::string rate_str = (boost::format("%.0f ksps") % (rate / 1e3)).str();
                        }
                        dev["label"]  = "VRT Instance " + std::to_string(i) + ": " +
                                        (boost::format("%.1f MHz") % (freq / 1e6)).str() +
                                        " (" + (boost::format("%.1f Msps") % (rate / 1e6)).str() + ")";
                        dev["rate"] = std::to_string(rate);
                        dev["vrt_instance"] = std::to_string(i);
                        instances.push_back(dev);
                        break;
                    }
                }
            }
        }
    }
    zmq_stopwatch_stop(deadline);

    for (void* subscriber : subscribers) zmq_close(subscriber);
    zmq_ctx_destroy(ctx);

    return instances;
}


SoapySDR::Device *makeVrtDevice(const SoapySDR::Kwargs &args)
{
    int instance = std::stoi(args.at("vrt_instance"));
    float freq = std::strtod(args.at("freq").c_str(), NULL);
    float rate = std::strtod(args.at("rate").c_str(), NULL);
    return new VrtDevice(instance, freq, rate);
}

static SoapySDR::Registry registerVrtDevice("vrt_device", &findVrtDevice, &makeVrtDevice, SOAPY_SDR_ABI_VERSION);
