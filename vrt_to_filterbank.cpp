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
#include <fftw3.h>

#include "vrt-tools.h"
#include "dt-extended-context.h"

namespace po = boost::program_options;

#define REAL 0
#define IMAG 1

#define SCALE_MAX 32768

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

    // FFTW
    fftw_complex *signal, *result;
    fftw_plan plan;

    float *magnitudes;

    FILE *write_ptr;

    // variables to be set by po
    std::string file, type, zmq_address, source_name, coords;
    uint16_t instance, main_port, port;
    uint32_t channel;
    uint32_t integrations;
    uint32_t num_bins = 0;
    uint32_t threads = 1;
    int32_t machine_id, telescope_id, data_type;
    int hwm;
    bool power2;
    size_t num_requested_samples;
    double total_time;
    float bin_size, integration_time;

    bool stream_has_dt_extended_context = false;

    // setup the program options
    po::options_description desc("Allowed options");
    // clang-format off

    desc.add_options()
        ("help", "help message")
        ("nsamps", po::value<size_t>(&num_requested_samples)->default_value(0), "total number of samples to receive")
        ("duration", po::value<double>(&total_time)->default_value(0), "total number of seconds to receive")
        ("file", po::value<std::string>(&file)->default_value("vrt.fil"), "name of the file to write binary filterbank data to")
        ("channel", po::value<uint32_t>(&channel)->default_value(0), "VRT channel")
        ("source-name", po::value<std::string>(&source_name)->default_value("not defined"), "Source name")
        ("coordinates", po::value<std::string>(&coords), "Coordinates (ra,dec,az,el)")
        // ("fft-duration", po::value<uint32_t>(&fft_len), "number of seconds to integrate")
        ("num-bins", po::value<uint32_t>(&num_bins)->default_value(10000), "number of bins")
        ("bin-size", po::value<float>(&bin_size), "size of bin in Hz")
        ("power2", po::value<bool>(&power2)->default_value(true), "Round number of bins to nearest power of two")
        ("integrations", po::value<uint32_t>(&integrations)->default_value(1), "number of integrations")
        ("integration-time", po::value<float>(&integration_time), "integration time (seconds)")
        ("threads", po::value<uint32_t>(&threads)->default_value(1), "enable multi-threading")
        ("machine-id", po::value<int32_t>(&machine_id)->default_value(0), "set filterbank machine_id (0=FAKE)")
        ("telescope-id", po::value<int32_t>(&telescope_id)->default_value(0), "set filterbank telescope_id (0=FAKE)")
        ("data-type", po::value<int32_t>(&data_type)->default_value(1), "set filterbank data_type (1=filterbank)")
        ("negative-foff", "negative foff frequency and set fch1 as highest channel")
        ("progress", "periodically display short-term bandwidth")
        // ("stats", "show average bandwidth on exit")
        ("int-second", "align start of reception to integer second")
        ("dt-trace", "use coordinates from DT trace data in VRT stream")
        ("null", "run without writing to file")
        ("continue", "don't abort on a bad packet")
        // ("ignore-dc", "Ignore  DC bin")
        ("address", po::value<std::string>(&zmq_address)->default_value("localhost"), "VRT ZMQ address")
        ("zmq-split", "create a ZeroMQ stream per VRT channel, increasing port number for additional streams")
        ("instance", po::value<uint16_t>(&instance)->default_value(0), "VRT ZMQ instance")
        ("port", po::value<uint16_t>(&port), "VRT ZMQ port")
        ("hwm", po::value<int>(&hwm)->default_value(10000), "VRT ZMQ HWM")
    ;
    // clang-format on
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // print the help message
    if (vm.count("help")) {
        std::cout << boost::format("VRT samples to filterbank. %s") % desc << std::endl;
        std::cout << std::endl
                  << "This application processes data from a VRT stream "
                     "to to filterbank format.\n"
                  << std::endl;
        return ~0;
    }

    bool progress               = vm.count("progress") > 0;
    bool stats                  = vm.count("stats") > 0;
    bool null                   = vm.count("null") > 0;
    bool continue_on_bad_packet = vm.count("continue") > 0;
    bool neg_foff               = vm.count("negative-foff") > 0;
    bool int_second             = (bool)vm.count("int-second");
    bool dt_trace               = vm.count("dt-trace") > 0;
    bool zmq_split              = vm.count("zmq-split") > 0;
    // bool ignore_dc              = (bool)vm.count("ignore-dc");

    std::vector<std::string> coord_strings;
    char *ptr;
    if (vm.count("coordinates")) {
        boost::split(coord_strings, coords, boost::is_any_of(", "),boost::token_compress_on);
        if (coord_strings.size()!=4) {
            printf("Incorrect number of coordinates. Exiting.\n");
            exit(1);
        }
        // for (size_t ch = 0; ch < coord_strings.size(); ch++) {
        //     printf("val: %lf\n", strtod(coord_strings[ch].c_str(), &ptr) );
        // }
    }

    write_ptr = fopen(file.c_str(),"wb");  // w for write, b for binary

    context_type vrt_context;
    dt_ext_context_type dt_ext_context;
    init_context(&vrt_context);

    packet_type vrt_packet;

    if (vm.count("port") > 0) {
        main_port = port;
    } else {
        main_port = DEFAULT_MAIN_PORT + MAX_CHANNELS*instance;
    }

    if (zmq_split) {
        main_port += channel;
        vrt_packet.channel_filt = 1;
    } else {
        vrt_packet.channel_filt = 1<<channel;
    }

    // ZMQ
    void *context = zmq_ctx_new();
    void *subscriber = zmq_socket(context, ZMQ_SUB);
    int rc = zmq_setsockopt (subscriber, ZMQ_RCVHWM, &hwm, sizeof hwm);
    std::string connect_string = "tcp://" + zmq_address + ":" + std::to_string(main_port);
    rc = zmq_connect(subscriber, connect_string.c_str());
    assert(rc == 0);
    zmq_setsockopt(subscriber, ZMQ_SUBSCRIBE, "", 0);

    // time keeping
    auto start_time = std::chrono::steady_clock::now();

    uint32_t buffer[ZMQ_BUFFER_SIZE];

    unsigned long long num_total_samps = 0;

    // Track time and samps between updating the BW summary
    auto last_update                     = start_time;
    unsigned long long last_update_samps = 0;

    bool first_frame = true;
    bool start_rx = false;
    uint64_t last_fractional_seconds_timestamp = 0;

    uint32_t signal_pointer = 0;
    uint32_t integration_counter = 0;

    int exit_code = EXIT_SUCCESS;
    while (not stop_signal_called
           and (num_requested_samples > num_total_samps or num_requested_samples == 0)) {

        int len = zmq_recv(subscriber, buffer, ZMQ_BUFFER_SIZE, 0);

        const auto now = std::chrono::steady_clock::now();

        if (not vrt_process(buffer, sizeof(buffer), &vrt_context, &vrt_packet)) {
            printf("Not a Vita49 packet?\n");
            continue;
        }

        if (not start_rx and vrt_packet.context) {
            vrt_print_context(&vrt_context);
            start_rx = true;

            if (vm.count("bin-size")) {
                if (power2) {
                    num_bins = (uint32_t)((float)vrt_context.sample_rate/(float)bin_size);
                    uint32_t pow2 = (uint32_t)(log2(num_bins)+0.8);
                    num_bins = pow(2,pow2);
                } else {
                    num_bins = (uint32_t)((float)vrt_context.sample_rate/(float)bin_size);
                }
            }

            if (vm.count("integration-time")) {
                integrations = (uint32_t)((double)integration_time/((double)num_bins/(double)vrt_context.sample_rate));
            }

            if (total_time > 0)
                num_requested_samples = total_time * vrt_context.sample_rate;

            fftw_plan_with_nthreads(threads);

            signal = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_bins);
            result = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_bins);
            plan = fftw_plan_dft_1d(num_bins, signal, result, FFTW_FORWARD, FFTW_ESTIMATE);
            magnitudes = (float*)malloc(num_bins * sizeof(float));
            memset(magnitudes, 0, num_bins*sizeof(float));

            printf("# Filterbank parameters:\n");
            printf("#    Bins: %u\n", num_bins);
            printf("#    Bin size [Hz]: %.0f\n", ((double)vrt_context.sample_rate)/((double)num_bins));
            printf("#    Integrations: %u\n", integrations);
            printf("#    Integration Time [sec]: %.4f\n", (double)integrations*(double)num_bins/(double)vrt_context.sample_rate);
        }

        if (start_rx and vrt_packet.data and (dt_ext_context.dt_ext_context_received or not dt_trace)) {

            if (vrt_packet.lost_frame)
               if (not continue_on_bad_packet) {
                    exit_code = 1;
                    break;
               }

            if (int_second) {
                // check if fractional second has wrapped
                if (vrt_packet.fractional_seconds_timestamp > last_fractional_seconds_timestamp ) {
                        last_fractional_seconds_timestamp = vrt_packet.fractional_seconds_timestamp;
                        continue;
                } else {
                    int_second = false;
                    last_update = now;
                    start_time = now;
                }
            }

            if (first_frame) {
                std::cout << boost::format(
                                 "# First frame: %u samples, %u full secs, %.09f frac secs")
                                 % vrt_packet.num_rx_samps
                                 % vrt_packet.integer_seconds_timestamp
                                 % ((double)vrt_packet.fractional_seconds_timestamp/1e12)
                          << std::endl;
                first_frame = false;

                // Filterbank Header
                const char* keyword;
                const char* string;
                int32_t len;
                int32_t int_value;
                double double_value;

                keyword = "HEADER_START";
                len = strlen(keyword);
                fwrite( &len, sizeof(len), 1, write_ptr);
                fwrite( (char*)keyword, len, 1, write_ptr);

                keyword = "machine_id";
                int_value = machine_id;
                len = strlen(keyword);
                fwrite( &len, sizeof(len), 1, write_ptr);
                fwrite( (char*)keyword, len, 1, write_ptr);
                fwrite( &int_value, sizeof(int_value), 1, write_ptr);

                keyword = "telescope_id";
                int_value = telescope_id;
                len = strlen(keyword);
                fwrite( &len, sizeof(len), 1, write_ptr);
                fwrite( (char*)keyword, len, 1, write_ptr);
                fwrite( &int_value, sizeof(int_value), 1, write_ptr);

                keyword = "data_type";
                int_value = data_type;
                len = strlen(keyword);
                fwrite( &len, sizeof(len), 1, write_ptr);
                fwrite( (char*)keyword, len, 1, write_ptr);
                fwrite( &int_value, sizeof(int_value), 1, write_ptr);

                keyword = "ibeam";
                int_value = 1;
                len = strlen(keyword);
                fwrite( &len, sizeof(len), 1, write_ptr);
                fwrite( (char*)keyword, len, 1, write_ptr);
                fwrite( &int_value, sizeof(int_value), 1, write_ptr);

                keyword = "source_name";
                len = strlen(keyword);
                fwrite( &len, sizeof(len), 1, write_ptr);
                fwrite( (char*)keyword, len, 1, write_ptr);
                len = strlen(source_name.c_str());
                fwrite( &len, sizeof(len), 1, write_ptr);
                fwrite( (char*)source_name.c_str(), len, 1, write_ptr);

                keyword = "nchans";
                int_value = num_bins;
                len = strlen(keyword);
                fwrite( &len, sizeof(len), 1, write_ptr);
                fwrite( (char*)keyword, len, 1, write_ptr);
                fwrite( &int_value, sizeof(int_value), 1, write_ptr);

                keyword = "nbeams";
                int_value = 1;
                len = strlen(keyword);
                fwrite( &len, sizeof(len), 1, write_ptr);
                fwrite( (char*)keyword, len, 1, write_ptr);
                fwrite( &int_value, sizeof(int_value), 1, write_ptr);

                keyword = "nbits";
                int_value = 32;
                len = strlen(keyword);
                fwrite( &len, sizeof(len), 1, write_ptr);
                fwrite( (char*)keyword, len, 1, write_ptr);
                fwrite( &int_value, sizeof(int_value), 1, write_ptr);

                keyword = "nifs";
                int_value = 1;
                len = strlen(keyword);
                fwrite( &len, sizeof(len), 1, write_ptr);
                fwrite( (char*)keyword, len, 1, write_ptr);
                fwrite( &int_value, sizeof(int_value), 1, write_ptr);

                keyword = "fch1";
                if (neg_foff)
                    double_value = (double)vrt_context.rf_freq/1e6+(double)vrt_context.sample_rate/2e6;
                else
                    double_value = (double)vrt_context.rf_freq/1e6-(double)vrt_context.sample_rate/2e6;
                len = strlen(keyword);
                fwrite( &len, sizeof(len), 1, write_ptr);
                fwrite( (char*)keyword, len, 1, write_ptr);
                fwrite( &double_value, sizeof(double_value), 1, write_ptr);

                keyword = "foff";
                if (neg_foff)
                    double_value = -((double)vrt_context.sample_rate/1e6)/((double)num_bins);
                else
                    double_value = ((double)vrt_context.sample_rate/1e6)/((double)num_bins);
                len = strlen(keyword);
                fwrite( &len, sizeof(len), 1, write_ptr);
                fwrite( (char*)keyword, len, 1, write_ptr);
                fwrite( &double_value, sizeof(double_value), 1, write_ptr);

                keyword = "tstart";
                double_value = ((double)vrt_packet.integer_seconds_timestamp+(double)vrt_packet.fractional_seconds_timestamp/1e12)/86400.0 + 40587.0;
                len = strlen(keyword);
                fwrite( &len, sizeof(len), 1, write_ptr);
                fwrite( (char*)keyword, len, 1, write_ptr);
                fwrite( &double_value, sizeof(double_value), 1, write_ptr);

                keyword = "tsamp";
                double_value = (double)integrations*(double)num_bins/(double)vrt_context.sample_rate;
                len = strlen(keyword);
                fwrite( &len, sizeof(len), 1, write_ptr);
                fwrite( (char*)keyword, len, 1, write_ptr);
                fwrite( &double_value, sizeof(double_value), 1, write_ptr);

                if (dt_trace) {
                    keyword = "src_raj";
                    double ra_h = ((12.0/M_PI)*dt_ext_context.ra_current);
                    int ra_hours = (int)ra_h;
                    int ra_minutes = (int)(ra_h*60)%60;
                    double ra_seconds = fmod(ra_h*3600.0, 60.0);
                    double_value = ra_hours*1e4 + ra_minutes*1e2 + ra_seconds;
                    len = strlen(keyword);
                    fwrite( &len, sizeof(len), 1, write_ptr);
                    fwrite( (char*)keyword, len, 1, write_ptr);
                    fwrite( &double_value, sizeof(double_value), 1, write_ptr);

                    keyword = "src_dej";
                    double dec_deg = ((180.0/M_PI)*dt_ext_context.dec_current);
                    int dec_degrees = (int)dec_deg;
                    int dec_minutes = (int)(dec_deg*60.0)%60;
                    double dec_seconds = fmod(dec_deg*3600, 60.0);
                    double_value = dec_degrees*1e4 + dec_minutes*1e2 + dec_seconds;
                    len = strlen(keyword);
                    fwrite( &len, sizeof(len), 1, write_ptr);
                    fwrite( (char*)keyword, len, 1, write_ptr);
                    fwrite( &double_value, sizeof(double_value), 1, write_ptr);

                    keyword = "az_start";
                    double_value = ((180.0/M_PI)*dt_ext_context.azimuth);
                    len = strlen(keyword);
                    fwrite( &len, sizeof(len), 1, write_ptr);
                    fwrite( (char*)keyword, len, 1, write_ptr);
                    fwrite( &double_value, sizeof(double_value), 1, write_ptr);

                    keyword = "za_start";
                    double_value = 90.0 - ((180.0/M_PI)*dt_ext_context.elevation);
                    len = strlen(keyword);
                    fwrite( &len, sizeof(len), 1, write_ptr);
                    fwrite( (char*)keyword, len, 1, write_ptr);
                    fwrite( &double_value, sizeof(double_value), 1, write_ptr);

                } else if (vm.count("coordinates")) {
                    keyword = "src_raj";
                    double_value = strtod(coord_strings[0].c_str(), &ptr);
                    len = strlen(keyword);
                    fwrite( &len, sizeof(len), 1, write_ptr);
                    fwrite( (char*)keyword, len, 1, write_ptr);
                    fwrite( &double_value, sizeof(double_value), 1, write_ptr);

                    keyword = "src_dej";
                    double_value = strtod(coord_strings[1].c_str(), &ptr);
                    len = strlen(keyword);
                    fwrite( &len, sizeof(len), 1, write_ptr);
                    fwrite( (char*)keyword, len, 1, write_ptr);
                    fwrite( &double_value, sizeof(double_value), 1, write_ptr);

                    keyword = "az_start";
                    double_value = strtod(coord_strings[2].c_str(), &ptr);
                    len = strlen(keyword);
                    fwrite( &len, sizeof(len), 1, write_ptr);
                    fwrite( (char*)keyword, len, 1, write_ptr);
                    fwrite( &double_value, sizeof(double_value), 1, write_ptr);

                    keyword = "za_start";
                    double_value = strtod(coord_strings[3].c_str(), &ptr);
                    len = strlen(keyword);
                    fwrite( &len, sizeof(len), 1, write_ptr);
                    fwrite( (char*)keyword, len, 1, write_ptr);
                    fwrite( &double_value, sizeof(double_value), 1, write_ptr);
                }

                keyword = "HEADER_END";
                len = strlen(keyword);
                fwrite( &len, sizeof(len), 1, write_ptr);
                fwrite( (char*)keyword, len, 1, write_ptr);
                // end header
            }

            int mult = 1;
            for (uint32_t i = 0; i < vrt_packet.num_rx_samps; i++) {

                int16_t re;
                memcpy(&re, (char*)&buffer[vrt_packet.offset+i], 2);
                int16_t img;
                memcpy(&img, (char*)&buffer[vrt_packet.offset+i]+2, 2);
                signal[signal_pointer][REAL] = mult*re;
                signal[signal_pointer][IMAG] = mult*img;
                mult *= -1;

                signal_pointer++;

                if (signal_pointer >= num_bins) {

                    signal_pointer = 0;

                    fftw_execute(plan);

                    for (uint32_t i = 0; i < num_bins; ++i) {
                        size_t index;
                        if (neg_foff)
                            index = num_bins-1-i;
                        else
                            index = i;
                        magnitudes[index] += (result[i][REAL] * result[i][REAL] +
                                  result[i][IMAG] * result[i][IMAG]);
                    }
                    integration_counter++;
                    if (integration_counter == integrations) {
                        for (uint32_t i = 0; i < num_bins; ++i)
                            magnitudes[i] /= (float)integrations;
                        fwrite(magnitudes, num_bins*sizeof(float), 1, write_ptr);
                        integration_counter = 0;
                        memset(magnitudes, 0, num_bins*sizeof(float));
                    }
                }
            }

            fflush(write_ptr);

            num_total_samps += vrt_packet.num_rx_samps;

        }

        if (vrt_packet.extended_context) {
            if (stream_has_dt_extended_context and not dt_trace) {
                std::cerr << "WARNING: DT metadata is present in the stream, but it is ignored. Did you forget --dt-trace?" << std::endl;
            }
            stream_has_dt_extended_context = dt_process(buffer, sizeof(buffer), &vrt_packet, &dt_ext_context);
        }

        if (progress) {
            if (vrt_packet.data)
                last_update_samps += vrt_packet.num_rx_samps;
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

                for (int i=0; i<vrt_packet.num_rx_samps; i++ ) {
                    auto sample_i = get_abs_val((std::complex<int16_t>)buffer[vrt_packet.offset+i]);
                    sum_i += sample_i;
                    if (sample_i > datatype_max*0.99)
                        clip_i++;
                }
                sum_i = sum_i/vrt_packet.num_rx_samps;
                std::cout << boost::format("%.0f") % (100.0*log2(sum_i)/log2(datatype_max)) << "% I (";
                std::cout << boost::format("%.0f") % ceil(log2(sum_i)+1) << " of ";
                std::cout << (int)ceil(log2(datatype_max)+1) << " bits), ";
                std::cout << "" << boost::format("%.0f") % (100.0*clip_i/vrt_packet.num_rx_samps) << "% I clip, ";
                std::cout << std::endl;

            }
        }
    }

    zmq_close(subscriber);
    zmq_ctx_destroy(context);

    fclose(write_ptr);

    return exit_code;
}
