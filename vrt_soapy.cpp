#include <SoapySDR/Device.hpp>
#include <SoapySDR/Registry.hpp>
#include <SoapySDR/Formats.hpp>
#include <iostream>
#include <memory>

struct DummyStream {};

class VrtDevice : public SoapySDR::Device
{
    public:

    VrtDevice(): SoapySDR::Device() {
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


SoapySDR::KwargsList findVrtDevice(const SoapySDR::Kwargs &args)
{
    (void)args;
    std::vector<SoapySDR::Kwargs> instances;
    SoapySDR::Kwargs dev;
    dev["driver"] = "vrt_device";
    dev["label"]  = "VRT Instance 1";
    instances.push_back(dev);
    return instances;
}


SoapySDR::Device *makeVrtDevice(const SoapySDR::Kwargs &args)
{
    (void)args;  // Translate args into constructor arguments here
    std::cout<<"Tammo says makeVrtDevice"<<std::endl;
    return new VrtDevice();
}


static SoapySDR::Registry registerVrtDevice("vrt_device", &findVrtDevice, &makeVrtDevice, SOAPY_SDR_ABI_VERSION);
