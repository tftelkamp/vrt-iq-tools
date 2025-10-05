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

#ifdef __APPLE__
#define DEFAULT_GNUPLOT_TERMINAL "qt"
#else
#define DEFAULT_GNUPLOT_TERMINAL "x11"
#endif

namespace po = boost::program_options;

#define REAL 0
#define IMAG 1

#define SQUELCH_THRESHOLD (0.02)

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

inline void swap(float *p,float *q) {
   float t;

   t=*p;
   *p=*q;
   *q=t;
}

inline void sort(float a[],int n) {
   int i,j;

   for(i = 0;i < n-1;i++) {
      for(j = 0;j < n-i-1;j++) {
         if(a[j] > a[j+1])
            swap(&a[j],&a[j+1]);
      }
   }
}

float dm_time(float dm, float freq_mhz) {
    // for a DM of 1 we expect 4148.8 usec delay at 1 GHZ.
    return 4148.8 * dm / (freq_mhz*freq_mhz);
}

int main(int argc, char* argv[])
{

    // FFTW
    fftw_complex **signal, *result;
    fftw_plan plan[2];

    float **mean_freq;
    float **mean_time;
    float ***data_block;

    float *median_freq;
    float *median_time;

    float **dedisp;
    int *dispersion;
    float **plotbuffer;

    FILE *audio_pipe;

    uint32_t num_bins = 0;

    float block_time;
    uint32_t block_size;

    float f_threshold;
    float t_threshold;

    // variables to be set by po
    std::string file, type, zmq_address, channel_list, gnuplot_terminal;
    size_t num_requested_samples;
    uint32_t bins;
    int gain;
    double total_time;
    uint16_t instance, main_port, port, pub_port;
    uint32_t channel;
    int hwm;
    float dm, period, agg_time;
    uint64_t seqno[] = {0, 0};
    float mean_block[] = {0, 0};
    int time_integrations;
    int buffer_size;
    float period_samples_float, amplitude;
    int period_samples_int;

    // setup the program options
    po::options_description desc("Allowed options");
    // clang-format off

    desc.add_options()
        ("help", "help message")
        ("nsamps", po::value<size_t>(&num_requested_samples)->default_value(0), "total number of samples to receive")
        ("duration", po::value<double>(&total_time)->default_value(0), "total number of seconds to receive")
        ("progress", "periodically display short-term bandwidth")
        // ("channel", po::value<uint32_t>(&channel)->default_value(0), "VRT channel")
        ("channel", po::value<std::string>(&channel_list)->default_value("0"), "which VRT channel(s) to use (specify \"0\", \"1\", \"0,1\", etc)")
        ("int-second", "align start of reception to integer second")
        ("num-bins", po::value<uint32_t>(&num_bins)->default_value(2000), "number of bins")
        ("block-time", po::value<float>(&block_time)->default_value(0.1), "block time (seconds)")
        ("f-threshold", po::value<float>(&f_threshold)->default_value(1.15), "frequency cut threshold")
        ("t-threshold", po::value<float>(&t_threshold)->default_value(1.2), "time cut threshold")
        ("dm", po::value<float>(&dm)->default_value(26.8), "PSR Dispersion Measure")
        ("period", po::value<float>(&period)->default_value(0.7145197), "PSR Period")
        ("agg-time", po::value<float>(&agg_time)->default_value(1), "Aggregation time in milliseconds")
        ("amplitude", po::value<float>(&amplitude)->default_value(1), "amplitude correction of second channel")
        ("term", po::value<std::string>(&gnuplot_terminal)->default_value(DEFAULT_GNUPLOT_TERMINAL), "Gnuplot terminal (x11 or qt)")
        ("zmq-pub", "enable zmq pub")
        ("no-stdout", "disable stdout")
        ("pub-port", po::value<uint16_t>(&pub_port)->default_value(60001), "ZMQ publisher port")
        ("quiet", "no data output")
        ("sum", "sum polarizations")
        ("audio", "enable audio")
        ("squelch", "audio squelch")
        ("gain", po::value<int>(&gain)->default_value(8), "audio gain")
        ("gnuplot", "enable gnuplot mode")
        ("null", "run without writing to file")
        ("continue", "don't abort on a bad packet")
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
        std::cout << boost::format("VRT pulsar processing. %s") % desc << std::endl;
        std::cout << std::endl
                  << "This application processes pulsar data from "
                     "a VRT stream.\n"
                  << std::endl;
        return ~0;
    }

    bool progress               = vm.count("progress") > 0;
    bool stats                  = vm.count("stats") > 0;
    bool null                   = vm.count("null") > 0;
    bool audio                  = vm.count("audio") > 0;
    bool gnuplot                = vm.count("gnuplot") > 0;
    bool sum                    = vm.count("sum") > 0;
    bool continue_on_bad_packet = vm.count("continue") > 0;
    bool quiet                  = vm.count("quiet") > 0;
    bool squelch                = vm.count("squelch") > 0;
    bool int_second             = (bool)vm.count("int-second");
    bool zmq_split              = vm.count("zmq-split") > 0;
    bool no_stdout              = vm.count("no-stdout") > 0;
    bool zmq_pub                = vm.count("zmq-pub") > 0;

    context_type vrt_context;
    init_context(&vrt_context);

    packet_type vrt_packet;

    if (vm.count("port") > 0) {
        main_port = port;
    } else {
        main_port = DEFAULT_MAIN_PORT + MAX_CHANNELS*instance;
    }

    vrt_packet.channel_filt = 0;

     // detect which channels to use
    std::vector<std::string> channel_strings;
    std::vector<size_t> channel_nums;
    boost::split(channel_strings, channel_list, boost::is_any_of("\"',"));
    for (size_t ch = 0; ch < channel_strings.size(); ch++) {
        size_t chan = std::stoi(channel_strings[ch]);
        channel_nums.push_back(std::stoi(channel_strings[ch]));
        vrt_packet.channel_filt |= 1<<std::stoi(channel_strings[ch]);
    }

    if (channel_nums.size() > 2) {
        printf("More than 2 channels not supported.\n");
        exit(1);
    }

    if (zmq_split) {
        if (channel_nums.size()>1) {
            printf("Multiple channels with --zmq-split is not supported.\n");
            exit(EXIT_FAILURE);
        }
        main_port += channel_nums[0];
        channel_nums[0] = 0;
        vrt_packet.channel_filt = 1;
    }

    // FILE *write_ptr;
    // write_ptr = fopen("dedisp.fc32","wb");  // w for write, b for binary

    // ZMQ

    void *context = zmq_ctx_new();
    void *zmq_server;
    void *subscriber = zmq_socket(context, ZMQ_SUB);
    int rc = zmq_setsockopt (subscriber, ZMQ_RCVHWM, &hwm, sizeof hwm);
    std::string connect_string = "tcp://" + zmq_address + ":" + std::to_string(main_port);
    rc = zmq_connect(subscriber, connect_string.c_str());
    assert(rc == 0);
    zmq_setsockopt(subscriber, ZMQ_SUBSCRIBE, "", 0);

    if (zmq_pub) {
        zmq_server = zmq_socket(context, ZMQ_PUB);
        rc = zmq_setsockopt(zmq_server, ZMQ_SNDHWM, &hwm, sizeof hwm);
        assert(rc == 0);
        connect_string = "tcp://*:" + std::to_string(pub_port);
        rc = zmq_bind(zmq_server, connect_string.c_str()) ;
        assert(rc == 0);
    }

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

    uint32_t signal_pointer[] = {0, 0};
    uint32_t block_counter[] = {0, 0};
    uint32_t integration_counter[] = {0, 0};

    bool first_block = true;

    while (not stop_signal_called
           and (num_requested_samples > num_total_samps or num_requested_samples == 0) ) {

        int len = zmq_recv(subscriber, buffer, ZMQ_BUFFER_SIZE, 0);

        const auto now = std::chrono::steady_clock::now();

        if (not vrt_process(buffer, sizeof(buffer), &vrt_context, &vrt_packet)) {
            printf("Not a Vita49 packet?\n");
            continue;
        }

        if (not vrt_packet.context and not vrt_packet.data)
            continue;

        uint32_t ch = 0;
        for(ch = 0; ch<channel_nums.size(); ch++)
            if (vrt_packet.stream_id & (1 << channel_nums[ch]) )
                break;

        uint32_t channel = channel_nums[ch];

        if (not start_rx and vrt_packet.context) {
            vrt_print_context(&vrt_context);
            start_rx = true;

            if (total_time > 0)
                num_requested_samples = total_time * vrt_context.sample_rate;

            block_size = block_time*vrt_context.sample_rate/num_bins;

            time_integrations = agg_time*(vrt_context.sample_rate/num_bins)/1000;

            signal = (fftw_complex **)malloc(sizeof(fftw_complex*)*channel_nums.size());

            for (size_t ch=0; ch < channel_nums.size(); ch++)
                signal[ch] = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_bins);

            result = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_bins);

            for (size_t ch=0; ch < channel_nums.size(); ch++)
                plan[ch] = fftw_plan_dft_1d(num_bins, signal[ch], result, FFTW_FORWARD, FFTW_ESTIMATE);

            data_block = (float ***)malloc(sizeof(float *)*channel_nums.size());

            for (size_t ch=0; ch < channel_nums.size(); ch++) {
                data_block[ch] = (float **)malloc(sizeof(float *)*num_bins);

                for(size_t i=0; i < num_bins; i++) {
                    data_block[ch][i] = (float *)malloc( 2 * sizeof(float)*block_size);
                }
            }

            mean_freq = (float **)malloc(sizeof(float *)*channel_nums.size());
            mean_time = (float **)malloc(sizeof(float *)*channel_nums.size());

            for (size_t ch=0; ch < channel_nums.size(); ch++) {
                mean_freq[ch] = (float*)malloc(num_bins * sizeof(float));
                mean_time[ch] = (float*)malloc(block_size * sizeof(float));
            }

            median_freq = (float*)malloc(num_bins * sizeof(float));
            median_time = (float*)malloc(block_size * sizeof(float));

            dedisp = (float **)malloc(sizeof(float *)*channel_nums.size());

            for (size_t ch=0; ch < channel_nums.size(); ch++)
                dedisp[ch] = (float*)malloc( (block_size/time_integrations) * sizeof(float));

            dispersion = (int*)malloc(num_bins*sizeof(int));

            // gnuplot
            buffer_size = 4 * (1000.0/agg_time); // 4 seconds

            plotbuffer = (float **)malloc(sizeof(float *)*channel_nums.size());
            for (size_t ch=0; ch < channel_nums.size(); ch++)
                plotbuffer[ch] = (float*)malloc(buffer_size * sizeof(float));

            // data
            period_samples_float = (1000.0/agg_time)*period;
            period_samples_int = floor(period_samples_float);

            // create dispersion table
            float freq_bin0 = (double)(vrt_context.rf_freq - vrt_context.sample_rate/2);
            float disp_bin0 = dm_time(dm,freq_bin0/1e6) * (float)vrt_context.sample_rate/(float(num_bins));

            for(size_t chan=0; chan < num_bins; chan++) {
                float freq = (double)(vrt_context.rf_freq + (chan*(double)vrt_context.sample_rate/(double)num_bins) - vrt_context.sample_rate/2);
                float disp = dm_time(dm,freq/1e6) * (float)vrt_context.sample_rate/(float(num_bins)) - disp_bin0;
                dispersion[chan] = (int)disp;
            }

            printf("# Spectrum parameters:\n");
            printf("#    Bins: %u\n", num_bins);
            printf("#    Bin size [Hz]: %.2f\n", ((double)vrt_context.sample_rate)/((double)num_bins));
            printf("#    Block size: %u\n", block_size);
            printf("#    Aggregations: %u\n", time_integrations);

            // Gnuplot
            if (gnuplot)
                printf("set terminal %s noraise; unset mouse; set grid;\n", gnuplot_terminal.c_str());
            // set terminal x11;
            // set yrange [0:200000000] set xtics 1; set ytics 1;

            if (audio) {

                uint32_t audio_rate = round(1000.0/agg_time);

                uint32_t num_chans = sum ? 1 : channel_nums.size();

                std::string sox_command = "play -q -r " + std::to_string(audio_rate)
                        + " --input-buffer 4000 --buffer 500 -c "
                        + std::to_string(num_chans)
                        + " -b 16 -e signed-integer -t raw - lowpass 200 rate 40k gain -l "
                        + std::to_string(gain);

                audio_pipe = popen(sox_command.c_str(), "w");

                if (!audio_pipe)
                {
                  printf("Error starting Sox play.\n");
                  return EXIT_FAILURE;
                }
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
                    int_second = false;
                    last_update = now;
                    start_time = now;
                }
            }

            if (first_block and (ch!=0) ) {
                continue;
            } else {
                first_block = false;
            }

            int mult = 1;
            for (uint32_t i = 0; i < vrt_packet.num_rx_samps; i++) {

                int16_t re;
                memcpy(&re, (char*)&buffer[vrt_packet.offset+i], 2);
                int16_t img;
                memcpy(&img, (char*)&buffer[vrt_packet.offset+i]+2, 2);
                if (ch==1) {
                    signal[ch][signal_pointer[ch]][REAL] = amplitude*mult*re;
                    signal[ch][signal_pointer[ch]][IMAG] = amplitude*mult*img;
                } else {
                    signal[ch][signal_pointer[ch]][REAL] = mult*re;
                    signal[ch][signal_pointer[ch]][IMAG] = mult*img;
                }
                mult *= -1; // fftshift

                signal_pointer[ch]++;

                if (signal_pointer[ch] >= num_bins) {

                    signal_pointer[ch] = 0;

                    fftw_execute(plan[ch]);

                    uint64_t seconds = vrt_packet.integer_seconds_timestamp;
                    uint64_t frac_seconds = vrt_packet.fractional_seconds_timestamp;
                    frac_seconds += (i+1)*1e12/vrt_context.sample_rate;
                    if (frac_seconds > 1e12) {
                        frac_seconds -= 1e12;
                        seconds++;
                    }

                    float sum_channels = 0;
                    for (uint32_t i = 0; i < num_bins; ++i) {
                        float mag = sqrt(result[i][REAL] * result[i][REAL] +
                                            result[i][IMAG] * result[i][IMAG]);
                        data_block[ch][i][block_size+block_counter[ch]] = mag;
                        mean_freq[ch][i] += mag/(float)block_size;
                        sum_channels += mag;
                    }

                    mean_time[ch][block_counter[ch]] = sum_channels/(float)num_bins;

                    block_counter[ch]++;

                    if (block_counter[ch] == block_size) {
                        // mow the lawn!
                        memcpy(median_freq, mean_freq[ch], num_bins * sizeof(float));
                        memcpy(median_time, mean_time[ch], block_size * sizeof(float));
                        sort(median_freq,num_bins);
                        sort(median_time,block_size);

                        float thresh_freq = f_threshold*median_freq[(num_bins+1)/2-1];
                        float thresh_time = t_threshold*median_time[(block_size+1)/2-1];

                        float freq_med = median_freq[(num_bins+1)/2-1];
                        float time_med = median_time[(block_size+1)/2-1];

                        int clean = 0;

                        for (size_t chan = 0; chan < num_bins; chan++)
                            for (size_t block = 0; block < block_size; block++) {
                                if ( mean_freq[ch][chan] > thresh_freq ) {
                                    data_block[ch][chan][block_size+block] = freq_med;
                                    clean++;
                                    continue;
                                }
                                if ( mean_time[ch][block] > thresh_time) {
                                    data_block[ch][chan][block_size+block] = time_med;
                                    clean++;
                                }
                            }

                        // now what?
                        // dedisperse and aggregate

                        for (size_t index = 0; index < block_size/time_integrations; index++) {
                            dedisp[ch][index] = 0;
                        }

                        for(size_t chan=0; chan < num_bins; chan++) {
                            for (size_t index = 0; index < block_size/time_integrations; index++) {
                                for (size_t j=0; j<time_integrations; j++) {
                                     dedisp[ch][index] += data_block[ch][chan][block_size+index*time_integrations+j+dispersion[chan]];
                                }
                            }
                        }

                        // for data analysis:
                        // fwrite(dedisp,sizeof(float)*block_size/time_integrations,1,write_ptr);

                        float max_block = 0;
                        mean_block[ch] = 0;

                        for (size_t index = 0; index < block_size/time_integrations; index++) {
                            // sum to avg
                            dedisp[ch][index] /= num_bins*(block_size/time_integrations);
                            // mean of block
                            mean_block[ch] += dedisp[ch][index];
                            if (dedisp[ch][index]>max_block)
                                max_block = dedisp[ch][index];
                        }
                        mean_block[ch] = mean_block[ch]/(block_size/time_integrations);

                        // if (!first_block) {
                            for (size_t index = 0; index < block_size/time_integrations; index++) {
                                plotbuffer[ch][seqno[ch] % buffer_size] = dedisp[ch][index];
                                if (!gnuplot and !quiet) {
                                    char message[512] = "";
                                    zmq_msg_t msg;
                                    if (channel_nums.size()==2) {
                                        if (ch==1) {
                                            if (sum) {
                                                snprintf(message, 512, "%i %i %f\n",period_samples_int,(int)floor(fmod(seqno[ch],period_samples_float)), dedisp[0][index] + dedisp[1][index]); 
                                            } else {
                                                snprintf(message, 512, "%i %i %f %f\n",period_samples_int,(int)floor(fmod(seqno[ch],period_samples_float)), dedisp[0][index], dedisp[1][index]);
                                            }
                                        }
                                    } else {
                                        snprintf(message, 512, "%i %i %f\n",period_samples_int,(int)floor(fmod(seqno[ch],period_samples_float)), dedisp[ch][index]);
                                    }
                                    if (strlen(message)>0) {
                                        // stdout
                                        if (!no_stdout)
                                            printf("%s",message);
                                        // ZMQ
                                        if (zmq_pub) {
                                            zmq_msg_init_size(&msg, strlen(message));
                                            memcpy(zmq_msg_data(&msg), message, strlen(message));
                                            zmq_msg_send(&msg, zmq_server, 0);
                                            zmq_msg_close(&msg);
                                        }
                                    }
                                }
                                if (audio){
                                    if (channel_nums.size()==2) {
                                        if (ch==1) {
                                            if (sum) {
                                                // write sum on single channel
                                                int16_t sample = 32768.0*(dedisp[0][index]-mean_block[0])/mean_block[0] +
                                                                 32768.0*(dedisp[1][index]-mean_block[1])/mean_block[1];
                                                if (squelch and sample < 2*SQUELCH_THRESHOLD*32768.0)
                                                    sample = 0;
                                                fwrite(&sample, sizeof(sample), 1, audio_pipe);
                                            } else {
                                                // write 2 channels
                                                int16_t sample = 32768.0*(dedisp[0][index]-mean_block[0])/mean_block[0];
                                                if (squelch and sample < SQUELCH_THRESHOLD*32768.0)
                                                    sample = 0;
                                                fwrite(&sample, sizeof(sample), 1, audio_pipe);
                                                sample = 32768.0*(dedisp[1][index]-mean_block[1])/mean_block[1];
                                                if (squelch and sample < SQUELCH_THRESHOLD*32768.0)
                                                    sample = 0;
                                                fwrite(&sample, sizeof(sample), 1, audio_pipe);
                                            }
                                        }
                                    } else {
                                        // single channel
                                        int16_t sample = 32768.0* (dedisp[ch][index]-mean_block[ch])/mean_block[ch];
                                        if (squelch and sample < SQUELCH_THRESHOLD*32768.0)
                                                    sample = 0;
                                        fwrite(&sample, sizeof(sample), 1, audio_pipe);
                                    }
                                }
                                seqno[ch]++;
                            }
                        // }
                        // first_block = false;

                        if (audio)
                            fflush(audio_pipe);
                        fflush(stdout);

                        // copy current data block to first position (for dedispersion of the next block)
                        for(size_t i=0; i < num_bins; i++) {
                            memcpy(&data_block[ch][i][0], &data_block[ch][i][block_size], sizeof(float)*block_size);
                        }

                        // gnuplot
                        if (gnuplot) {

                            float mean_plot_buffer = 0;
                            for (i = 0; i< buffer_size; i++) {
                                mean_plot_buffer += plotbuffer[ch][i];
                            }
                            mean_plot_buffer /= buffer_size;

                            float time_per_sample = vrt_context.sample_rate/(num_bins*time_integrations);

                            if (channel_nums.size()==2) {
                                if (seqno[0] == seqno[1]) {
                                    printf("set xrange [%.2lf:%.2lf];\n", seqno[0]/time_per_sample, (seqno[0] + buffer_size)/time_per_sample);
                                    printf("set yrange [%.2lf:%.2lf];\n", mean_plot_buffer*0.97, mean_plot_buffer*1.5);
                                    printf("plot '-' u 1:2 notitle w l, '-' u 1:3 notitle w l\n");
                                    for (i = 0; i< buffer_size; i++) {
                                        printf("%lf\t%lf\t%lf\n",(seqno[0]+i)/time_per_sample, plotbuffer[0][(seqno[0]+i)%buffer_size], plotbuffer[1][(seqno[1]+i)%buffer_size] );
                                    }
                                    printf("e\n");
                                }
                            } else {
                                printf("set xrange [%.2lf:%.2lf];\n",seqno[ch]/time_per_sample, (seqno[ch] + buffer_size)/time_per_sample);
                                printf("set yrange [%.2lf:%.2lf];\n", mean_plot_buffer*0.97, mean_plot_buffer*1.5);
                                printf("plot '-' u 1:2 notitle w l\n");
                                for (i = 0; i< buffer_size; i++)
                                    printf("%lf\t%lf\n",(seqno[ch]+i)/time_per_sample, plotbuffer[ch][(seqno[ch]+i)%buffer_size]);
                                printf("e\n");
                            }
                        }

                        // clean-up
                        memset(mean_freq[ch], 0 , num_bins * sizeof(float));
                        block_counter[ch] = 0;
                    }

                }
            }

            num_total_samps += vrt_packet.num_rx_samps;

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

    return 0;

}
