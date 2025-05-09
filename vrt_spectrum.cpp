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

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/posix_time_io.hpp>

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
#include "tracker-extended-context.h"

#ifdef __APPLE__
#define DEFAULT_GNUPLOT_TERMINAL "qt"
#else
#define DEFAULT_GNUPLOT_TERMINAL "x11"
#endif

namespace po = boost::program_options;

#define NUM_POINTS 10000
#define REAL 0
#define IMAG 1

const double pi = std::acos(-1.0);

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

float haversine(float dec1, float dec2, float ra1, float ra2) {

    float dec_delta = dec2 - dec1;
    float ra_delta = ra2 - ra1;

    float a =
      pow(sin(dec_delta / 2), 2) + cos(dec1) * cos(dec2) * pow(sin(ra_delta / 2), 2);
    float c = 2 * atan2(sqrt(a), sqrt(1 - a));
    return(c);
}

float bearing(float dec1, float dec2, float ra1, float ra2) {

    float b = atan2(cos(dec1)*sin(dec2)-sin(dec1)*cos(dec2)*cos(ra2-ra1), sin(ra2-ra1)*cos(dec2));
    return b;
}

int main(int argc, char* argv[])
{

    // FFTW
    fftw_complex *signal, *result;
    fftw_plan plan;
    double *magnitudes, *filter_out;

    std::complex<float> *wola_buffer;
    float *wola_taps;

    uint32_t num_points = 0;
    uint32_t num_bins = 0;
    uint32_t wola_partitions;

    bool power2;
    float bin_size, integration_time = 0.0;
    float binsize;
    double alpha, tau;
    double min_offset, max_offset;
    uint32_t output_counter = 0;
    int32_t min_bin, max_bin;

    // variables to be set by po
    std::string file, type, zmq_address, gnuplot_terminal, gnuplot_commands, source;
    size_t num_requested_samples;
    uint32_t bins, updates_per_second;
    double total_time;
    uint32_t integrations, num_integrations;
    uint16_t instance, main_port, port;
    uint32_t channel;
    int hwm;

    bool dt_trace_warning_given = false;

    float min_y = 1e10;
    float max_y = -1e10;

    FILE *outfile;

    std::vector<double> poly;

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
        ("channel", po::value<uint32_t>(&channel)->default_value(0), "VRT channel")
        // ("stats", "show average bandwidth on exit")
        ("int-second", "align start of reception to integer second")
        ("int-interval", "align start of reception to integer number of integration intervals (implies --int-second)")
        ("num-bins", po::value<uint32_t>(&num_bins)->default_value(1000), "number of bins")
        ("bin-size", po::value<float>(&bin_size), "size of bin in Hz")
        ("power2", po::value<bool>(&power2)->default_value(false), "Round number of bins to nearest power of two")
        ("integrations", po::value<uint32_t>(&integrations), "number of integrations")
        ("integration-time", po::value<float>(&integration_time)->default_value(1.0), "integration time (seconds)")
        ("num-integrations", po::value<uint32_t>(&num_integrations)->default_value(0), "exit after number of integrations")
        ("tau", po::value<double>(&tau), "Exponential weighted moving average time constant (sec)")
        ("poly", po::value<std::vector<double> >(&poly)->multitoken(), "Polynomal coefficients to compensate bandpass")
        ("source", po::value<std::string>(&source), "Source description (ECSV and gnuplot)")
        ("gnuplot", "Gnuplot mode")
        ("fftmax", "fftmax mode")
        ("two", "square signal before processing (to detect BPSK signals)")
        ("four", "square-square signal before processing (to detect QPSK signals")
        ("wola", "apply Weighted OverLap Add method")
        ("wola-partitions", po::value<uint32_t>(&wola_partitions)->default_value(4), "number of WOLA partitions")
        ("min-offset", po::value<double>(&min_offset), "min. freq. offset to track (Hz)")
        ("max-offset", po::value<double>(&max_offset), "max. freq. offset to track (Hz)")
        ("gnuplot-commands", po::value<std::string>(&gnuplot_commands)->default_value(""), "Extra gnuplot commands like \"set yr [ymin:ymax];\"")
        ("term", po::value<std::string>(&gnuplot_terminal)->default_value(DEFAULT_GNUPLOT_TERMINAL), "Gnuplot terminal (x11 or qt)")
        ("minmax", "min/max hold for y-axis scale (gnuplot)")
        ("db", "output power in dB")
        ("dc", "suppress DC peak")
        ("ecsv", "output in ECSV format (Astropy)")
        ("bin-file", po::value<std::string>(&file), "output binary data to file")
        ("center-freq", "output center frequency")
        ("temperature", "output temperature")
        ("null", "run without writing to file")
        ("continue", "don't abort on a bad packet")
        ("dt-trace", "use DT trace data in VRT stream")
        ("address", po::value<std::string>(&zmq_address)->default_value("localhost"), "VRT ZMQ address")
        ("zmq-split", "create a ZeroMQ stream per VRT channel, increasing port number for additional streams")
        ("instance", po::value<uint16_t>(&instance)->default_value(0), "VRT ZMQ instance")
        ("port", po::value<uint16_t>(&port), "VRT ZMQ port")
        ("hwm", po::value<int>(&hwm)->default_value(10000), "VRT ZMQ HWM")
    ;
    // clang-format on
    po::variables_map vm;
    // po::store(po::parse_command_line(argc, argv, desc), vm);
    auto parsed = po::command_line_parser(argc, argv).options(desc).positional({}).style(po::command_line_style::unix_style ^ po::command_line_style::allow_short).run();
    po::store(parsed, vm);

    // print the help message
    if (vm.count("help")) {
        std::cout << boost::format("VRT samples to spectra. %s") % desc << std::endl;
        std::cout << std::endl
                  << "This application generates spectra from a VRT stream.\n"
                  << std::endl;
        return ~0;
    }
    po::notify(vm);

    bool progress               = vm.count("progress") > 0;
    bool stats                  = vm.count("stats") > 0;
    bool null                   = vm.count("null") > 0;
    bool continue_on_bad_packet = vm.count("continue") > 0;
    bool int_interval           = (bool)vm.count("int-interval");
    bool int_second             = int_interval || (bool)vm.count("int-second");
    bool dt_trace               = vm.count("dt-trace") > 0;
    bool db                     = vm.count("db") > 0;
    bool log_freq               = vm.count("center-freq") > 0;
    bool log_temp               = vm.count("temperature") > 0;
    bool gnuplot                = vm.count("gnuplot") > 0;
    bool fftmax                 = vm.count("fftmax") > 0;
    bool poly_calib             = vm.count("poly") > 0;
    bool iir                    = vm.count("tau") > 0;
    bool minmax                 = vm.count("minmax") > 0;
    bool ecsv                   = vm.count("ecsv") > 0;
    bool binary                 = vm.count("bin-file") > 0;
    bool dc                     = vm.count("dc") > 0;
    bool has_source             = vm.count("source") > 0;
    bool zmq_split              = vm.count("zmq-split") > 0;
    bool wola                   = vm.count("wola") > 0;
    bool flag_x2                = vm.count("two") > 0;
    bool flag_x4                = vm.count("four") > 0;  

    if (iir) {
        alpha = (1.0 - exp(-1/(tau/integration_time)));
    }

    if (int_interval && int(integration_time) == 0) {
        throw(std::runtime_error("--int-interval requires --integration_time > 1"));
    }

    double freq_div;

    if (flag_x2 || flag_x4) {
        freq_div = (flag_x4) ? 4 : 2;
    } else {
        freq_div = 1;
    }

    context_type vrt_context;
    dt_ext_context_type dt_ext_context;
    tracker_ext_context_type tracker_ext_context;
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

    bool first_frame = true;

    // time keeping
    auto start_time = std::chrono::steady_clock::now();
    auto stop_time =
        start_time + std::chrono::milliseconds(int64_t(1000 * total_time));

    uint32_t buffer[ZMQ_BUFFER_SIZE];

    unsigned long long num_total_samps = 0;

    // Track time and samps between updating the BW summary
    auto last_update                     = start_time;
    unsigned long long last_update_samps = 0;

    bool start_rx = false;
    uint64_t last_fractional_seconds_timestamp = 0;

    uint32_t signal_pointer = 0;
    uint32_t integration_counter = 0;
    uint32_t num_integrations_counter = 0;

    if (binary) {
        outfile=fopen(file.c_str(),"w");
    }

    while (not stop_signal_called
           and (num_requested_samples > num_total_samps or num_requested_samples == 0)
           and (total_time == 0.0 or std::chrono::steady_clock::now() <= stop_time)) {

        int len = zmq_recv(subscriber, buffer, ZMQ_BUFFER_SIZE, 0);

        const auto now = std::chrono::steady_clock::now();

        if (not vrt_process(buffer, sizeof(buffer), &vrt_context, &vrt_packet)) {
            printf("Not a Vita49 packet?\n");
            continue;
        }

        if (not start_rx and vrt_packet.context) {
            if (!ecsv)
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

            if (not vm.count("integrations")) {
                integrations = (uint32_t)round((double)integration_time/((double)num_bins/(double)vrt_context.sample_rate));
            }

            if (total_time > 0)
                num_requested_samples = total_time * vrt_context.sample_rate;

            min_bin = 0;
            max_bin = num_bins;

            binsize = ((double)vrt_context.sample_rate)/((double)num_bins);

            if (vm.count("min-offset")) {
                min_bin = (min_offset/binsize)+num_bins/2;
                min_bin = min_bin < 0 ? 0 : min_bin;
                min_bin = min_bin > num_bins ? num_bins : min_bin;
            }

            if (vm.count("max-offset")) {
                max_bin = (max_offset/binsize)+num_bins/2;
                max_bin = max_bin < 0 ? 0 : max_bin;
                max_bin = max_bin > num_bins ? num_bins : max_bin;
            }

            signal = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_bins);
            result = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_bins);
            plan = fftw_plan_dft_1d(num_bins, signal, result, FFTW_FORWARD, FFTW_ESTIMATE);
            magnitudes = (double*)malloc(num_bins * sizeof(double));
            memset(magnitudes, 0, num_bins*sizeof(double));
            filter_out = (double*)malloc(num_bins * sizeof(double));
            memset(filter_out, 0, num_bins*sizeof(double));

            if (wola) {
                int wola_len = wola_partitions*num_bins;
                wola_buffer = (std::complex<float>*)malloc(sizeof(std::complex<float>) * wola_len);
                wola_taps = (float*)malloc(sizeof(float) * wola_len);

                for (uint32_t i = 0; i < wola_len; i++)
                    wola_buffer[i] = std::complex<float>(0,0);

                // Blackman-Harris window
                double a0 = 0.35875;
                double a1 = 0.48829;
                double a2 = 0.14128;
                double a3 = 0.01168;

                for (int i=0;i<wola_len;i++) {
                    int j = i - (wola_len)/2;

                    double blackman_window = a0 - a1*cos(2*pi*(double)i/((double)wola_len-1)) + 
                                                a2*cos(4*pi*(double)i/((double)wola_len-1)) + 
                                                a3*cos(6*pi*(double)i/((double)wola_len-1));
                 
                    double x = ((double)wola_partitions*pi)*((double)j/(double)wola_len);
                    if (x!=0)
                        wola_taps[i] = blackman_window*sin(x)/(x);
                    else
                        wola_taps[i] = blackman_window*1.0;
                }
            }

            if (!ecsv) {
                printf("# Spectrum parameters:\n");
                printf("#    Bins: %u\n", num_bins);
                printf("#    Bin size [Hz]: %.2f\n", binsize);
                printf("#    Integrations: %u\n", integrations);
                printf("#    Integration Time [sec]: %.2f\n", (double)integrations*(double)num_bins/(double)vrt_context.sample_rate);
            } else {
                uint32_t first_col = 1;
                if (log_freq) first_col++;
                if (log_temp) first_col++;
                if (dt_trace) first_col += 15;

                printf("# %%ECSV 1.0\n");
                printf("# ---\n");

                uint32_t ch=0;
                while(not (vrt_context.stream_id & (1 << ch) ) )
                    ch++;
                printf("# delimiter: \',\'\n");
                printf("# meta: !!omap\n");
                printf("# - vrt: !!omap\n");
                printf("#   - {stream_id: %u}\n", vrt_context.stream_id);
                printf("#   - {channel: %u}\n", ch);
                printf("#   - {sample_rate: %.1f}\n", (float)vrt_context.sample_rate);
                printf("#   - {frequency: %.1f}\n", (double)vrt_context.rf_freq);
                printf("#   - {bandwidth: %.1f}\n", (float)vrt_context.bandwidth);
                printf("#   - {rx_gain: %.1f}\n", (float)vrt_context.gain);
                printf("#   - {reference: %s}\n", vrt_context.reflock == 1 ? "external" : "internal");
                printf("#   - {time_source: %s}\n", vrt_context.time_cal == 1? "pps" : "internal");
                printf("# - spectrum: !!omap\n");
                printf("#   - {db: %s}\n", db ? "True" : "False");
                printf("#   - {bins: %u}\n", num_bins);
                printf("#   - {col_first_bin: %u}\n", first_col);
                printf("#   - {bin_size: %.2f}\n", ((double)vrt_context.sample_rate)/((double)num_bins));
                printf("#   - {integrations: %u}\n", integrations);
                printf("#   - {integration_time: %.2f}\n", (double)integrations*(double)num_bins/(double)vrt_context.sample_rate);
                if (has_source) {
                    printf("# - description: !!omap\n");
                    printf("#   - {source: %s}\n", source.c_str());
                }

                printf("# datatype:\n");
                printf("# - {name: timestamp, datatype: float64}\n");
                if (log_freq) {
                    printf("# - {name: center_freq_hz, unit: Hz, datatype: float64}\n");
                }
                if (log_temp) {
                    printf("# - {name: temperature_deg_c, datatype: float64}\n");
                }
                if (dt_trace) {
                    printf("# - {name: current_az_deg, unit: deg, datatype: float64}\n");
                    printf("# - {name: current_el_deg, unit: deg, datatype: float64}\n");
                    printf("# - {name: current_az_error_deg, unit: deg, datatype: float64}\n");
                    printf("# - {name: current_el_error_deg, unit: deg, datatype: float64}\n");
                    printf("# - {name: current_az_speed_deg, unit: deg, datatype: float64}\n");
                    printf("# - {name: current_el_speed_deg, unit: deg, datatype: float64}\n");
                    printf("# - {name: current_az_offset_deg, unit: deg, datatype: float64}\n");
                    printf("# - {name: current_el_offset_deg, unit: deg, datatype: float64}\n");
                    printf("# - {name: current_ra_h, unit: h, datatype: float64}\n");
                    printf("# - {name: current_dec_deg, unit: deg, datatype: float64}\n");
                    printf("# - {name: setpoint_ra_h, unit: deg, datatype: float64}\n");
                    printf("# - {name: setpoint_dec_deg, unit: deg, datatype: float64}\n");
                    printf("# - {name: radec_error_angle_deg, unit: deg, datatype: float64}\n");
                    printf("# - {name: radec_error_bearing_deg, unit: deg, datatype: float64}\n");
                    printf("# - {name: focusbox_mm, unit: mm, datatype: float64}\n");
                }
                if (fftmax) {
                    printf("# - {name: max_frequency, unit: Hz, datatype: float64}\n");
                    printf("# - {name: max_power, datatype: float64}\n");
                } else {
                    for (uint32_t i = 0; i < num_bins; ++i) {
                            printf("# - {name: \'%.0f\', datatype: float64}\n", (double)((double)vrt_context.rf_freq + (i*binsize - vrt_context.sample_rate/2)/freq_div));
                    }
                }
                printf("# schema: astropy-2.0\n");
            }

            // Header
            if (!gnuplot) {
                printf("timestamp");
                if (log_freq)
                    printf(", center_freq_hz");
                if (log_temp)
                    printf(", temperature_deg_c");
                if (dt_trace)
                    printf(", current_az_deg, current_el_deg, current_az_error_deg, current_el_error_deg, current_az_speed_deg, current_el_speed_deg, current_az_offset_deg, current_el_offset_deg, current_ra_h, current_dec_deg, setpoint_ra_h, setpoint_dec_deg, radec_error_angle_deg, radec_error_bearing_deg, focusbox_mm");
                if (fftmax) {
                    printf(", max_frequency, max_power");
                } else {
                    for (uint32_t i = 0; i < num_bins; ++i) {
                            printf(", %.0f", (double)((double)vrt_context.rf_freq + (i*binsize - vrt_context.sample_rate/2)/freq_div));
                    }
                }
                printf("\n");
                fflush(stdout);
            }
        }

        if (start_rx and vrt_packet.data) {

            if (vrt_packet.lost_frame)
               if (not continue_on_bad_packet)
                    break;

            if (int_second) {
                // check if fractional second has wrapped
                if (vrt_packet.fractional_seconds_timestamp > last_fractional_seconds_timestamp ) {
                        last_fractional_seconds_timestamp = vrt_packet.fractional_seconds_timestamp;
                        continue;
                } else {
                    if (int_interval && int(vrt_packet.integer_seconds_timestamp) % int(integration_time) != 0) {
                        continue;
                    }
                    int_second = false;
                    last_update = now;
                    start_time = now;
                }
            }

            int mult = 1;
            for (uint32_t i = 0; i < vrt_packet.num_rx_samps; i++) {

                int16_t re;
                memcpy(&re, (char*)&buffer[vrt_packet.offset+i], 2);
                int16_t img;
                memcpy(&img, (char*)&buffer[vrt_packet.offset+i]+2, 2);

                if (wola) {
                    wola_buffer[signal_pointer+((wola_partitions-1)*num_bins)] = std::complex<float>(mult*re,mult*img);
                    mult *= -1;
                } else {
                    signal[signal_pointer][REAL] = mult*re;
                    signal[signal_pointer][IMAG] = mult*img;
                    mult *= -1;
                }

                signal_pointer++;

                if (signal_pointer >= num_bins) {

                    signal_pointer = 0;

                    uint64_t seconds = vrt_packet.integer_seconds_timestamp;
                    uint64_t frac_seconds = vrt_packet.fractional_seconds_timestamp;
                    frac_seconds += (i+1)*1e12/vrt_context.sample_rate;
                    if (frac_seconds > 1e12) {
                        frac_seconds -= 1e12;
                        seconds++;
                    }

                    // (double) square signal
                    if (flag_x2 || flag_x4) {
                        int mult = 1;
                        for (uint32_t i = 0; i < num_bins; i++) {

                            double real = signal[i][REAL];
                            double imag = signal[i][IMAG];

                            double real2 = real * real - imag * imag;
                            double imag2 = 2 * real * imag;

                            if (flag_x4) {
                                signal[i][REAL] = mult*(real2 * real2 - imag2 * imag2);
                                signal[i][IMAG] = mult*(2 * real2 * imag2);
                            } else {
                                signal[i][REAL] = mult*real2;
                                signal[i][IMAG] = mult*imag2;
                            }
                            mult *= -1;
                        }
                    }

                    if (num_bins > 1) {
                        if (wola) {
                            for (int j=0; j < num_bins; j++) {
                                signal[j][REAL] = 0;
                                signal[j][IMAG] = 0;
                            }
                            for (int p=0; p<wola_partitions; p++) {
                                for (int j=0; j < num_bins; j++) {
                                    signal[j][REAL] += wola_taps[p*num_bins+j] * wola_buffer[p*num_bins+j].real();
                                    signal[j][IMAG] += wola_taps[p*num_bins+j] * wola_buffer[p*num_bins+j].imag();
                                }
                            }
                            fftw_execute(plan);

                            // shift wola buffer
                            memcpy(&wola_buffer[0], &wola_buffer[num_bins], (wola_partitions-1)*num_bins*sizeof(std::complex<float>));

                        } else {
                            fftw_execute(plan);
                        }

                        for (uint32_t i = 0; i < num_bins; ++i) {
                            magnitudes[i] += (result[i][REAL] * result[i][REAL] +
                                      result[i][IMAG] * result[i][IMAG]);
                        }
                    } else {
                        magnitudes[0] += (signal[0][REAL] * signal[0][REAL] +
                                      signal[0][IMAG] * signal[0][IMAG]);
                    }

                    if (dc) {
                        size_t dcbin = num_bins/2;
                        magnitudes[dcbin] = (magnitudes[dcbin-1]+magnitudes[dcbin+1])/2;
                    }

                    integration_counter++;
                    if (integration_counter == integrations) {
                        num_integrations_counter++;
                        if (!gnuplot) {
                            if (binary) {
                                double timestamp = (double)seconds + (double)(frac_seconds/1e12);
                                fwrite(&timestamp,sizeof(double),1,outfile);
                            } else {
                                printf("%lu.%09li", static_cast<unsigned long>(seconds), static_cast<long>(frac_seconds/1e3));
                            }
                            if (log_freq) {
                                if (not binary) {
                                    printf(", %li", static_cast<long>(vrt_context.rf_freq));
                                }
                                else {
                                    double freq = vrt_context.rf_freq;
                                    fwrite(&freq,sizeof(double),1,outfile);
                                }
                            }
                            if (log_temp) {
                                if (not binary) {
                                    printf(", %.2f", vrt_context.temperature);
                                } else {
                                    double temp = vrt_context.temperature;
                                    fwrite(&temp,sizeof(double),1,outfile);
                                }
                            }
                            if (dt_trace) {
                                if (not binary) {
                                    printf(", %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f",
                                        ((180.0/M_PI)*dt_ext_context.azimuth),
                                        ((180.0/M_PI)*dt_ext_context.elevation),
                                        ((180.0/M_PI)*dt_ext_context.azimuth_error),
                                        ((180.0/M_PI)*dt_ext_context.elevation_error),
                                        ((180.0/M_PI)*dt_ext_context.azimuth_speed),
                                        ((180.0/M_PI)*dt_ext_context.elevation_speed),
                                        ((180.0/M_PI)*dt_ext_context.azimuth_offset),
                                        ((180.0/M_PI)*dt_ext_context.elevation_offset),
                                        ((12.0/M_PI)*dt_ext_context.ra_current),
                                        ((180.0/M_PI)*dt_ext_context.dec_current),
                                        ((12.0/M_PI)*dt_ext_context.ra_setpoint),
                                        ((180.0/M_PI)*dt_ext_context.dec_setpoint),
                                        ((180.0/M_PI)*haversine(dt_ext_context.dec_setpoint, dt_ext_context.dec_current, dt_ext_context.ra_setpoint, dt_ext_context.ra_current)),
                                        ((180.0/M_PI)*bearing(dt_ext_context.dec_setpoint, dt_ext_context.dec_current, dt_ext_context.ra_setpoint, dt_ext_context.ra_current)),
                                        dt_ext_context.focusbox);
                                } else {
                                    double trace_values[15];
                                    trace_values[0] = ((180.0/M_PI)*dt_ext_context.azimuth);
                                    trace_values[1] = ((180.0/M_PI)*dt_ext_context.elevation);
                                    trace_values[2] = ((180.0/M_PI)*dt_ext_context.azimuth_error);
                                    trace_values[3] = ((180.0/M_PI)*dt_ext_context.elevation_error);
                                    trace_values[4] = ((180.0/M_PI)*dt_ext_context.azimuth_speed);
                                    trace_values[5] = ((180.0/M_PI)*dt_ext_context.elevation_speed);
                                    trace_values[6] = ((180.0/M_PI)*dt_ext_context.azimuth_offset);
                                    trace_values[7] = ((180.0/M_PI)*dt_ext_context.elevation_offset);
                                    trace_values[8] = ((12.0/M_PI)*dt_ext_context.ra_current);
                                    trace_values[9] = ((180.0/M_PI)*dt_ext_context.dec_current);
                                    trace_values[10] = ((12.0/M_PI)*dt_ext_context.ra_setpoint);
                                    trace_values[11] = ((180.0/M_PI)*dt_ext_context.dec_setpoint);
                                    trace_values[12] = ((180.0/M_PI)*haversine(dt_ext_context.dec_setpoint, dt_ext_context.dec_current, dt_ext_context.ra_setpoint, dt_ext_context.ra_current));
                                    trace_values[13] = ((180.0/M_PI)*bearing(dt_ext_context.dec_setpoint, dt_ext_context.dec_current, dt_ext_context.ra_setpoint, dt_ext_context.ra_current));
                                    trace_values[14] = dt_ext_context.focusbox;
                                    fwrite(&trace_values,15*sizeof(double),1,outfile);
                                }
                            }

                            int N = poly.size();
                            output_counter++;

                            double max_power = 0; // change this to minimal double
                            int32_t max_i = -1;

                            double value;
                            // uint32_t dc = num_points/2;

                            for (uint32_t i = 0; i < num_bins; ++i) {
                                magnitudes[i] /= (double)integrations;

                                if (iir) {
                                    double current_alpha = (1.0/(float)output_counter > alpha) ? 1.0/(float)output_counter : alpha;
                                    filter_out[i] += (double)current_alpha*(magnitudes[i]-filter_out[i]);
                                } else {
                                    filter_out[i] = magnitudes[i];
                                }

                                double offset = i*binsize - vrt_context.sample_rate/2;

                                double correction = 1;

                                if  (poly_calib) {
                                    correction = 0;
                                    for (int32_t p = 0; p < N; p++) {
                                        correction += poly[p] * pow(offset, int(N-p-1));
                                    }
                                }
                                if (!fftmax) {
                                    if (db) {
                                        correction = 10*log10(correction);
                                        value = 10*log10(filter_out[i])-correction;
                                        if (not binary) {
                                            printf(", %.3f", value);
                                        } else {
                                            fwrite(&value,sizeof(double),1,outfile);
                                        }
                                    } else {
                                        value = filter_out[i]/correction;
                                        if (not binary) {
                                            printf(", %.3f", value);
                                        } else {
                                            fwrite(&value,sizeof(double),1,outfile);
                                        }
                                    }
                                } else {
                                    if (db) {
                                        correction = 10*log10(correction);
                                        value = 10*log10(filter_out[i])-correction;
                                    } else {
                                        value = filter_out[i]/correction;
                                    }
                                    if ( (value > max_power) and (i >= min_bin) and (i <= max_bin) and not (dc && i==num_bins/2)) {
                                        max_power = value;
                                        max_i = i;
                                    }
                        
                                }
                            }
                            if (fftmax) {
                                printf(", %.2f", (double)vrt_context.rf_freq + (max_i*binsize - vrt_context.sample_rate/2)/freq_div);
                                printf(", %.3f", max_power);
                            }
                            if (not binary)
                                printf("\n");
                        } else {
                            // gnuplot
                            double max_power = -1e10; // change this to minimal double
                            double max_freq = -1;

                            float scale = 1e6; // MHz

                            float ticks = (vrt_context.sample_rate/freq_div)/(4*scale);
                            printf("set term %s 1 noraise; set xtics %f; set xlabel \"Frequency (MHz)\"; set ylabel \"Power (dB)\"; ", gnuplot_terminal.c_str(), ticks);
                            printf("%s; ", gnuplot_commands.c_str());
                            if (has_source) {
                                printf("set title \"%s (%sZ)\" font \",14\"\n;", source.c_str(), 
                                    boost::posix_time::to_iso_extended_string(boost::posix_time::from_time_t(seconds)).c_str());
                            } else if (tracker_ext_context.tracker_ext_context_received) {
                                printf("set title \"%s (%sZ)\" font \",14\"\n;", tracker_ext_context.object_name, 
                                    boost::posix_time::to_iso_extended_string(boost::posix_time::from_time_t(seconds)).c_str());
                            }
                                
                            if (vrt_context.sample_rate <= 100e3)
                                printf("set format x \"%%.4f\";\n");
                            else
                                printf("set format x \"%%.3f\";\n");
                            if (minmax && (min_y < max_y))
                                printf("set yr [%f:%f];", min_y, max_y);
                            else
                                printf("set offsets 0, 0, 0.2, 0.2;");
                            if (!fftmax)
                                printf("plot \"-\" u 1:2 with lines title \"signal\";\n");
                            else
                                printf("plot \"-\" u 1:2 with lines title \"signal\", \"-\" u 1:2 with points title \"max\";\n");

                            int N = poly.size();
                            output_counter++;

                            for (uint32_t i = 0; i < num_bins; ++i) {
                                magnitudes[i] /= (double)integrations;

                                if (iir) {
                                    double current_alpha = (1.0/(float)output_counter > alpha) ? 1.0/(float)output_counter : alpha;
                                    filter_out[i] += (double)current_alpha*(magnitudes[i]-filter_out[i]);
                                } else {
                                    filter_out[i] = magnitudes[i];
                                }
                                double offset = i*binsize - vrt_context.sample_rate/2;
                                double freq = ((double)vrt_context.rf_freq + offset/freq_div)/scale;

                                double correction = 0;

                                if  (poly_calib) {
                                    for (int32_t p = 0; p < N; p++) {
                                        correction += poly[p] * pow(offset, int(N-p-1));
                                    }
                                    correction = 10*log10(correction);
                                }
                                float value = 10*log10(filter_out[i])-correction;
                                if (minmax) {
                                    min_y = (value < min_y) ? value : min_y;
                                    max_y = (value > max_y) ? value : max_y;
                                }
                                printf("%.6f, %.6f\n", freq, value);
                                if ( (value > max_power) and (i >= min_bin) and (i <= max_bin) and not (dc && i==num_bins/2)) {
                                        max_power = value;
                                        max_freq = freq;
                                }
                            }
                            printf("e\n");
                            if (fftmax) {
                                int numdigits = 6-int(log10(binsize));
                                numdigits = numdigits < 0 ? 0 : numdigits;
                                printf("%.6f, %.6f\n", max_freq, max_power);
                                printf("e\n");
                                printf("snr_str = sprintf(\"%.1f\");\n", max_power);
                                printf("set label 1 snr_str at %.6f, %.6f offset 1,-1;\n", max_freq, max_power);
                                printf("freq_str = sprintf(\"%.*f\");\n", numdigits, max_freq);
                                printf("set label 2 freq_str at %.6f, %.6f offset 1,-2;\n", max_freq, max_power);
                            }
                        }

                        integration_counter = 0;
                        memset(magnitudes, 0, num_bins*sizeof(double));
                        if (binary)
                            fflush(outfile);
                        else
                            fflush(stdout);
                    }
                    if ( (num_integrations > 0) && (num_integrations_counter == num_integrations))
                        stop_signal_called = true;
                }
            }

            num_total_samps += vrt_packet.num_rx_samps;

        }

        if (vrt_packet.extended_context) {
            if (not dt_trace_warning_given and not gnuplot and dt_ext_context.dt_ext_context_received and not dt_trace) {
                std::cerr << "WARNING: DT metadata is present in the stream, but it is ignored. Did you forget --dt-trace?" << std::endl;
                dt_trace_warning_given = true;
            }
            dt_process(buffer, sizeof(buffer), &vrt_packet, &dt_ext_context);
            tracker_process(buffer, sizeof(buffer), &vrt_packet, &tracker_ext_context);   
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

    if (binary)
        fclose(outfile);

    zmq_close(subscriber);
    zmq_ctx_destroy(context);

    return 0;

}
