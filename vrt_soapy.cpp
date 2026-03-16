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

    VrtDevice(int instance): SoapySDR::Device() {
        int port = DEFAULT_MAIN_PORT + MAX_CHANNELS * instance;
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


    SoapySDR::RangeList getSampleRateRange(int, size_t) const override
    {
        std::cout<<"Tammo says getSampleRateRange"<<std::endl;
        return { SoapySDR::Range(1e6, 1e6) };
    }

    double getSampleRate(int, size_t) const override
    {
        std::cout<<"Tammo says getSampleRate"<<std::endl;
        return 1e6;
    }

    std::string getDriverKey() const override {
        std::cout<<"Tammo says getDriverKey"<<std::endl;
        return "vrt_device";
    }

    std::string getHardwareKey() const override {
        std::cout<<"Tammo says getHardwareKey"<<std::endl;
        return "VRTZMQ";
    }


    SoapySDR::RangeList getFrequencyRange(int, size_t) const override
    {
        std::cout<<"Tammo says getFrequencyRange"<<std::endl;
        return { SoapySDR::Range(0, 6e9) };
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
        std::cout << "Tammo says readStream" << std::endl;

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
        return 100e6;
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
        return {1e6};
    }
};



std::vector<int> discover_instances() {
    std::vector<int> found_instances;
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
                found_instances.push_back(i);
                remaining--;
                items[i].events = 0;

                while (true) {
                    int len = zmq_recv(subscribers[i], buffer, ZMQ_BUFFER_SIZE, 0);
                    if (len <= 0) {std::cerr<<"Small len"<<std::endl; continue;}
                    if (not vrt_process(buffer, sizeof(buffer), &vrt_context, &vrt_packet)) {
                        std::cerr<<"Not a Vita49 packet?"<<std::endl;
                        break;
                    }
                    if (vrt_packet.context) {
                        std::cout<<"Tammo says sample rate: "<<(float)vrt_context.sample_rate<<std::endl;
                        break;
                    }
                }
            }
        }
    }
    zmq_stopwatch_stop(deadline);

    for (void* subscriber : subscribers) zmq_close(subscriber);
    zmq_ctx_destroy(ctx);

    std::sort(found_instances.begin(), found_instances.end());
    return found_instances;
}

SoapySDR::KwargsList findVrtDevice(const SoapySDR::Kwargs &args)
{
    (void)args;
    std::vector<SoapySDR::Kwargs> instances;

    std::vector<int> vrt_instances = discover_instances();

    for (auto vrt_instance : vrt_instances) {
        SoapySDR::Kwargs dev;
        dev["driver"] = "vrt_device";
        dev["label"]  = "VRT Instance " + std::to_string(vrt_instance);
        dev["vrt_instance"] = std::to_string(vrt_instance);
        if (!args.count("vrt_instance") || args.at("vrt_instance") == std::to_string(vrt_instance)) {
            instances.push_back(dev);
        }
    }
    return instances;
}


SoapySDR::Device *makeVrtDevice(const SoapySDR::Kwargs &args)
{
    for (auto &kv : args) {
        std::cerr << "makeVrtDevice arg: " << kv.first << " = " << kv.second << std::endl;
    }
    std::cout<<"Tammo says makeVrtDevice "<< args.at("label") << std::endl;
    return new VrtDevice(std::stoi(args.at("vrt_instance")));
}


static SoapySDR::Registry registerVrtDevice("vrt_device", &findVrtDevice, &makeVrtDevice, SOAPY_SDR_ABI_VERSION);
