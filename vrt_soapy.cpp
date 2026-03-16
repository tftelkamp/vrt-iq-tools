#include <SoapySDR/Device.hpp>
#include <SoapySDR/Registry.hpp>
#include <SoapySDR/Formats.hpp>
#include <iostream>

#include <zmq.h>
#include "vrt-tools.h"
#include <memory>
#include <vector>

struct DummyStream {};

class VrtDevice : public SoapySDR::Device
{
    public:

    VrtDevice(int vrt_instance, float freq, float rate): SoapySDR::Device(),
        vrt_instance_(vrt_instance), freq_(freq), rate_(rate) {
        int port = DEFAULT_MAIN_PORT + MAX_CHANNELS * vrt_instance_;
        std::cout<<"Tammo says creating VrtDevice at port " <<port<<std::endl;
    };

    size_t getNumChannels(const int dir) const override
    {
        return (dir == SOAPY_SDR_RX) ? 1 : 0;
    }

    std::vector<std::string> getStreamFormats(int, size_t) const override
    {
        std::cout<<"Tammo says getStreamFormats"<<std::endl;
        return {SOAPY_SDR_CS16, SOAPY_SDR_CF32};
    }

    SoapySDR::Stream *setupStream(int, const std::string &,
                                  const std::vector<size_t> &,
                                  const SoapySDR::Kwargs &
                                 ) override
    {
        std::cout<<"Tammo says setupStream"<<std::endl;
        return reinterpret_cast<SoapySDR::Stream *>(new DummyStream());
    }

    void closeStream(SoapySDR::Stream *stream) override {
        std::cout<<"Tammo says closeStream"<<std::endl;
        delete reinterpret_cast<DummyStream *>(stream);
    }

/*
    SoapySDR::RangeList getSampleRateRange(int, size_t) const override
    {
        std::cout<<"Tammo says getSampleRateRange"<<std::endl;
        return { SoapySDR::Range(de1e6, 1e6) };
    }
*/

    double getSampleRate(int, size_t) const override
    {
        return rate_;
    }

    std::string getDriverKey() const override {
        return "vrt_device";
    }

    std::string getHardwareKey() const override {
        std::cout<<"Tammo says getHardwareKey"<<std::endl;
        return "VRTZMQ";
    }


    std::vector<std::string> listFrequencies(
    const int direction,
    const size_t channel) const override
    {
        std::cout<<"Tammo says listFrequencies"<<std::endl;
        return {"RF"};
    }

    SoapySDR::RangeList getFrequencyRange(int, size_t) const override
    {
        std::cout<<"Tammo says getFrequencyRange"<<std::endl;
        return { SoapySDR::Range(freq_, freq_) };
    }

    void setSampleRate(int, size_t, double) override {}

    int activateStream(
        SoapySDR::Stream *, int, long long, size_t) override
    {
        std::cout << "Tammo says activateStream" << std::endl;
        return 0;
    }

    int deactivateStream(
        SoapySDR::Stream *, int, long long) override
    {
        std::cout << "Tammo says deactivateStream" << std::endl;
        return 0;
    }

    int readStream(
        SoapySDR::Stream *, void * const *buffs, size_t numElems,
        int &flags, long long &timeNs, long timeoutUs
    ) override
    {
        //std::cout << "Tammo says readStream" << std::endl;

        flags = 0;
        timeNs = 0;
        std::memset(buffs[0], 0, numElems * sizeof(int16_t) * 2);
        return numElems;
    }

    bool hasHardwareTime(const std::string &) const override
    {
        return false;
    }

    std::vector<std::string> listAntennas(int, size_t) const override
    {
        return {};
    }

    void setAntenna(int, size_t, const std::string &) override {}

    std::string getAntenna(int, size_t) const override
    {
        return "";
    }

    SoapySDR::ArgInfoList getSettingInfo() const override
    {
        return {};
    }

    double getFrequency(int, size_t) const override
    {
        return freq_;
    }

    void setFrequency(
        int dir,
        size_t channel,
        double frequency,
        const SoapySDR::Kwargs &args
    ) override
    {
    }

    std::vector<double> listSampleRates(int dir, size_t channel) const override
    {
        std::cout << "Tammo says listSampleRates" << std::endl;
        return {rate_};
    }
private:
    float freq_;
    float rate_;
    int vrt_instance_;
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
        void* subscriber = zmq_socket(ctx, ZMQ_SUB);
        int rc = zmq_setsockopt (subscriber, ZMQ_RCVHWM, &hwm, sizeof hwm);
        assert(rc == 0);
        zmq_setsockopt(subscriber, ZMQ_SUBSCRIBE, "", 0);

        std::string addr = "tcp://localhost:" + std::to_string(port);
        zmq_connect(subscriber, addr.c_str());  // Asynchronous
        subscribers.push_back(subscriber);
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
                    std::cout<<"Not interested in this instance!"<<std::endl;
                    continue;
                }

                while (true) {
                    int len = zmq_recv(subscribers[i], buffer, ZMQ_BUFFER_SIZE, 0);
                    if (len <= 0) {std::cerr<<"Small len"<<std::endl; continue;}
                    if (not vrt_process(buffer, sizeof(buffer), &vrt_context, &vrt_packet)) {
                        std::cerr<<"Not a Vita49 packet?"<<std::endl;
                        break;
                    }
                    if (vrt_packet.context) {
                        float rate = vrt_context.sample_rate;
                        float freq = vrt_context.rf_freq;
                        dev["freq"] = std::to_string(freq);
                        dev["driver"] = "vrt_device";

                        std::string rate_str = (boost::format("%.1f Msps") % (rate / 1e6)).str();
                        if (rate < 1e6)
                            std::string rate_str = (boost::format("%.0f ksps") % (rate / 1e3)).str();

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
    for (auto &kv : args) {
        std::cerr << "makeVrtDevice arg: " << kv.first << " = " << kv.second << std::endl;
    }

    int instance = std::stoi(args.at("vrt_instance"));
    float freq = std::strtod(args.at("freq").c_str(), NULL);
    float rate = std::strtod(args.at("rate").c_str(), NULL);
    std::cout<<"Tammo says makeVrtDevice "<< args.at("label") << std::endl;
    std::cout<<"MakeVrtDevice freq: "<< freq << std::endl;
    std::cout<<"MakeVrtDevice rate: "<< rate << std::endl;
    return new VrtDevice(instance, freq, rate);
}


static SoapySDR::Registry registerVrtDevice("vrt_device", &findVrtDevice, &makeVrtDevice, SOAPY_SDR_ABI_VERSION);
