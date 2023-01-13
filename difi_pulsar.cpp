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

#include "difi-tools.h"

namespace po = boost::program_options;

#define REAL 0
#define IMAG 1

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
    fftw_complex *signal, *result;
    fftw_plan plan;
    float *magnitudes;
    float *mean_freq;
    float *mean_time;
    float **data_block;

    float *median_freq;
    float *median_time;

    float *dedisp;
    int *dispersion;
    float *plotbuffer;

    FILE *audio_pipe;

    uint32_t num_bins = 0;

    float block_time;
    uint32_t block_size;

    float f_threshold;
    float t_threshold;

    // variables to be set by po
    std::string file, type, zmq_address;
    size_t num_requested_samples;
    uint32_t bins;
    double total_time;
    uint16_t port;
    uint32_t channel;
    int hwm;
    float dm, period, agg_time;
    uint64_t seqno;
    int time_integrations;
    int buffer_size;
    float period_samples_float;
    int period_samples_int;

    // setup the program options
    po::options_description desc("Allowed options");
    // clang-format off

    desc.add_options()
        ("help", "help message")
        ("nsamps", po::value<size_t>(&num_requested_samples)->default_value(0), "total number of samples to receive")
        ("duration", po::value<double>(&total_time)->default_value(0), "total number of seconds to receive")
        ("progress", "periodically display short-term bandwidth")
        ("channel", po::value<uint32_t>(&channel)->default_value(0), "DIFI channel")
        ("int-second", "align start of reception to integer second")
        ("num-bins", po::value<uint32_t>(&num_bins)->default_value(2000), "number of bins")
        ("block-time", po::value<float>(&block_time)->default_value(0.2), "block time (seconds)")
        ("f-threshold", po::value<float>(&f_threshold)->default_value(1.15), "frequency cut threshold")
        ("t-threshold", po::value<float>(&t_threshold)->default_value(1.2), "time cut threshold")
        ("dm", po::value<float>(&dm)->default_value(26.8), "PSR Dispersion Measure")
        ("period", po::value<float>(&period)->default_value(0.7145197), "PSR Period")
        ("agg-time", po::value<float>(&agg_time)->default_value(1), "Aggregation time in milliseconds")
        ("audio", "enable audio")
        ("gnuplot", "enable gnuplot mode")
        ("null", "run without writing to file")
        ("continue", "don't abort on a bad packet")
        ("address", po::value<std::string>(&zmq_address)->default_value("localhost"), "DIFI ZMQ address")
        ("port", po::value<uint16_t>(&port)->default_value(50100), "DIFI ZMQ port")
        ("hwm", po::value<int>(&hwm)->default_value(10000), "DIFI ZMQ HWM")

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
    bool continue_on_bad_packet = vm.count("continue") > 0;
    bool int_second             = (bool)vm.count("int-second");

    context_type difi_context;
    init_context(&difi_context);

    difi_packet_type difi_packet;

    difi_packet.channel_filt = 1<<channel;

    // FILE *write_ptr;
    // write_ptr = fopen("dedisp.fc32","wb");  // w for write, b for binary

    // ZMQ

    void *context = zmq_ctx_new();
    void *subscriber = zmq_socket(context, ZMQ_SUB);
    int rc = zmq_setsockopt (subscriber, ZMQ_RCVHWM, &hwm, sizeof hwm);
    std::string connect_string = "tcp://" + zmq_address + ":" + std::to_string(port);
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
    uint32_t block_counter = 0;
    uint32_t integration_counter = 0;

    bool first_block = true;

    while (not stop_signal_called
           and (num_requested_samples > num_total_samps or num_requested_samples == 0) ) {

        int len = zmq_recv(subscriber, buffer, ZMQ_BUFFER_SIZE, 0);

        const auto now = std::chrono::steady_clock::now();

        if (not difi_process(buffer, sizeof(buffer), &difi_context, &difi_packet)) {
            printf("Not a Vita49 packet?\n");
            continue;
        }

        if (not start_rx and difi_packet.context) {
            difi_print_context(&difi_context);
            start_rx = true;

            if (total_time > 0)  
                num_requested_samples = total_time * difi_context.sample_rate;

            block_size = block_time*difi_context.sample_rate/num_bins;

            time_integrations = agg_time*(difi_context.sample_rate/num_bins)/1000; 

            signal = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_bins);
            result = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_bins);
            plan = fftw_plan_dft_1d(num_bins, signal, result, FFTW_FORWARD, FFTW_ESTIMATE);

            data_block = (float **)malloc(sizeof(float *)*num_bins);

            for(size_t i=0; i < num_bins; i++) {
                data_block[i] = (float *)malloc( 2 * sizeof(float)*block_size); 
            }

            magnitudes = (float*)malloc(num_bins * sizeof(float));
            mean_freq = (float*)malloc(num_bins * sizeof(float));
            mean_time = (float*)malloc(block_size * sizeof(float));
            median_freq = (float*)malloc(num_bins * sizeof(float));
            median_time = (float*)malloc(block_size * sizeof(float));

            dedisp = (float*)malloc( (block_size/time_integrations) * sizeof(float));

            dispersion = (int*)malloc(num_bins*sizeof(int));

            // gnuplot
            buffer_size = 5 * (1000.0/agg_time)*period; // 5 intervals
            seqno = 0; //buffer_size;
            plotbuffer = (float*)malloc(buffer_size * sizeof(float));

            // data
            period_samples_float = (1000.0/agg_time)*period;
            period_samples_int = floor(period_samples_float);

            // create dispersion table
            float freq_bin0 = (double)(difi_context.rf_freq - difi_context.sample_rate/2);
            float disp_bin0 = dm_time(dm,freq_bin0/1e6) * (float)difi_context.sample_rate/(float(num_bins));

            for(size_t chan=0; chan < num_bins; chan++) {
                float freq = (double)(difi_context.rf_freq + (chan*(double)difi_context.sample_rate/(double)num_bins) - difi_context.sample_rate/2);
                float disp = dm_time(dm,freq/1e6) * (float)difi_context.sample_rate/(float(num_bins)) - disp_bin0;
                dispersion[chan] = (int)disp;
            }

            printf("# Spectrum parameters:\n");
            printf("#    Bins: %u\n", num_bins);
            printf("#    Bin size [Hz]: %.2f\n", ((double)difi_context.sample_rate)/((double)num_bins));
            printf("#    Block size: %u\n", block_size);
            printf("#    Aggregations: %u\n", time_integrations);

            // Gnuplot
            if (gnuplot)
                printf("set terminal x11; unset mouse; set grid;\n");
            // set terminal x11; 
            // set yrange [0:200000000] set xtics 1; set ytics 1;

            if (audio) {

                uint32_t audio_rate = round(1000.0/agg_time);

                std::string sox_command = "play -q -r " + std::to_string(audio_rate) + " --buffer 200 -c 1 -b 16 -e signed-integer -v 30 -t raw - lowpass 40 gain 10";
        
                audio_pipe = popen(sox_command.c_str(), "w");
            
                if (!audio_pipe)
                {
                  printf("Error starting Sox play.\n");
                  return EXIT_FAILURE;
                }
            }
        }

        if (start_rx and difi_packet.data) {

            if (difi_packet.lost_frame)
               if (not continue_on_bad_packet)
                    break;

            if (int_second) {
                // check if fractional second has wrapped
                if (difi_packet.fractional_seconds_timestamp > last_fractional_seconds_timestamp ) {
                        last_fractional_seconds_timestamp = difi_packet.fractional_seconds_timestamp;
                        continue;
                } else {
                    int_second = false;
                    last_update = now; 
                    start_time = now;
                }
            }

            int mult = 1;
            for (uint32_t i = 0; i < difi_packet.num_rx_samps; i++) {

                int16_t re;
                memcpy(&re, (char*)&buffer[difi_packet.offset+i], 2);
                int16_t img;
                memcpy(&img, (char*)&buffer[difi_packet.offset+i]+2, 2);
                signal[signal_pointer][REAL] = mult*re;
                signal[signal_pointer][IMAG] = mult*img;
                mult *= -1; // fftshift

                signal_pointer++;

                if (signal_pointer >= num_bins) { 

                    signal_pointer = 0;

                    fftw_execute(plan);

                    uint64_t seconds = difi_packet.integer_seconds_timestamp;
                    uint64_t frac_seconds = difi_packet.fractional_seconds_timestamp;
                    frac_seconds += (i+1)*1e12/difi_context.sample_rate;
                    if (frac_seconds > 1e12) {
                        frac_seconds -= 1e12;
                        seconds++;
                    }

                    float sum_channels = 0;
                    for (uint32_t i = 0; i < num_bins; ++i) {
                        float mag = sqrt(result[i][REAL] * result[i][REAL] +
                                            result[i][IMAG] * result[i][IMAG]);
                        data_block[i][block_size+block_counter] = mag;
                        mean_freq[i] += mag/(float)block_size;
                        sum_channels += mag;
                    }

                    mean_time[block_counter] = sum_channels/(float)num_bins;

                    block_counter++;

                    if (block_counter == block_size) {
                        // mow the lawn!
                        memcpy(median_freq, mean_freq, num_bins * sizeof(float));
                        memcpy(median_time, mean_time, block_size * sizeof(float));
                        sort(median_freq,num_bins);
                        sort(median_time,block_size);

                        float thresh_freq = f_threshold*median_freq[(num_bins+1)/2-1];
                        float thresh_time = t_threshold*median_time[(block_size+1)/2-1];

                        float freq_med = median_freq[(num_bins+1)/2-1];
                        float time_med = median_time[(block_size+1)/2-1];

                        int clean = 0;

                        for (size_t chan = 0; chan < num_bins; chan++)
                            for (size_t block = 0; block < block_size; block++) {
                                if ( mean_freq[chan] > thresh_freq ) {
                                    data_block[chan][block_size+block] = freq_med;
                                    clean++;
                                    continue;
                                }
                                if ( mean_time[block] > thresh_time) {
                                    data_block[chan][block_size+block] = time_med;
                                    clean++;
                                }
                            }

                        // now what?
                        // dedisperse and aggregate

                        for (size_t index = 0; index < block_size/time_integrations; index++) {
                            dedisp[index] = 0;
                        }

                        for(size_t chan=0; chan < num_bins; chan++) {
                            for (size_t index = 0; index < block_size/time_integrations; index++) {
                                for (size_t j=0; j<time_integrations; j++) {
                                     dedisp[index] += data_block[chan][block_size+index*time_integrations+j+dispersion[chan]];
                                }   
                            }
                        }

                        // for data analysis:
                        // fwrite(dedisp,sizeof(float)*block_size/time_integrations,1,write_ptr); 

                        float mean_block = 0;
                        float max_block = 0;

                        for (size_t index = 0; index < block_size/time_integrations; index++) {
                            // sum to avg 
                            dedisp[index] /= num_bins*(block_size/time_integrations);
                            // mean of block
                            mean_block += dedisp[index];
                            if (dedisp[index]>max_block)
                                max_block = dedisp[index];
                        }
                        mean_block = mean_block/(block_size/time_integrations);

                        if (!first_block) {
                            for (size_t index = 0; index < block_size/time_integrations; index++) {
                                plotbuffer[seqno % buffer_size] = dedisp[index];
                                if (!gnuplot)
                                    printf("%i %i %f\n",period_samples_int,(int)floor(fmod(seqno,period_samples_float)),dedisp[index]);
                                if (audio){
                                    int16_t sample = 2*32768.0* (dedisp[index]-mean_block)/mean_block;
                                    fwrite(&sample, sizeof(sample), 1, audio_pipe);
                                }
                                seqno++;
                            }
                        }
                        first_block = false;

                        if (audio)
                            fflush(audio_pipe);
                        fflush(stdout);

                        // copy current data block to first position (for dedispersion of the next block)
                        for(size_t i=0; i < num_bins; i++) {
                            memcpy(&data_block[i][0], &data_block[i][block_size], sizeof(float)*block_size);
                        }

                        // gnuplot
                        if (gnuplot) {
                            printf("set xrange [%.0lf:%.0lf];\n",seqno/period, (seqno + buffer_size)/period);
                            printf("plot '-' notitle w l\n");

                            for (i = 0; i< buffer_size; i++) {
                                // integration_counter = (seqno-BUFFER_SIZE+i)/BUFFER_SIZE;
                                printf("%lf\t%lf\n",(seqno+i)/period, plotbuffer[(seqno +i)% buffer_size] );
                            }
                            printf("e\n");
                        }

                        // clean-up
                        memset(mean_freq, 0 , num_bins * sizeof(float));
                        block_counter = 0;
                    }

                }
            }

            num_total_samps += difi_packet.num_rx_samps;

        }

        if (progress) {
            if (difi_packet.data)
                last_update_samps += difi_packet.num_rx_samps;
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
  
                for (int i=0; i<difi_packet.num_rx_samps; i++ ) {
                    auto sample_i = get_abs_val((std::complex<int16_t>)buffer[difi_packet.offset+i]);
                    sum_i += sample_i;
                    if (sample_i > datatype_max*0.99)
                        clip_i++;
                }
                sum_i = sum_i/difi_packet.num_rx_samps;
                std::cout << boost::format("%.0f") % (100.0*log2(sum_i)/log2(datatype_max)) << "% I (";
                std::cout << boost::format("%.0f") % ceil(log2(sum_i)+1) << " of ";
                std::cout << (int)ceil(log2(datatype_max)+1) << " bits), ";
                std::cout << "" << boost::format("%.0f") % (100.0*clip_i/difi_packet.num_rx_samps) << "% I clip, ";
                std::cout << std::endl;

            }
        }
    }

    zmq_close(subscriber);
    zmq_ctx_destroy(context);

    return 0;

}  
