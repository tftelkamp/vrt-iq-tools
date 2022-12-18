#include <zmq.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <fftw3.h>
#include <getopt.h>
#include <time.h>
#include <sys/time.h>

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include <chrono>
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

#include "difi-tools.h"

namespace po = boost::program_options;

#define SCALE_MAX 32768

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

void usage(void)
{
  printf("rffft: FFT RF observations\n\n");
  printf("-i <file>       Input file (can be fifo) [stdin]\n");
  printf("-p <prefix>     Output prefix\n");
  printf("-o <output>     Output filename [default: YYYY-MM-DDTHH:MM:SS.sss_XXXXXX.bin]\n");
  printf("-f <frequency>  Center frequency (Hz)\n");
  printf("-s <samprate>   Sample rate (Hz)\n");
  printf("-c <chansize>   Channel size [100Hz]\n");
  printf("-t <tint>       Integration time [1s]\n");
  printf("-n <nsub>       Number of integrations per file [60]\n");
  printf("-m <use>        Use every mth integration [1]\n");
  printf("-F <format>     Input format char, int, float [int]\n");
  printf("-T <start time> YYYY-MM-DDTHH:MM:SSS.sss\n");
  printf("-R <fmin,fmax>  Frequency range to store (Hz)\n");
  printf("-S <index>      Starting index [int]\n");
  printf("-I              Invert frequencies\n");
  printf("-b              Digitize output to bytes [off]\n");
  printf("-q              Quiet mode, no output [off]\n");
  printf("-h              This help\n");

  return;
}

int main(int argc, char* argv[])
{
  int i,j,k,l,nchan,m=0,nint=1,nsub=60,flag,nuse=1,imin,imax,partial=0;
  fftwf_complex *c,*d;
  fftwf_plan fft;
  FILE *outfile;
  char outfname[128]="",prefix[32]="";
  char outformat='f';
  char *cbuf;
  float *fbuf;
  float *z,length,fchan=100.0,tint=1.0,zavg,zstd,*zw;
  char *cz;
  double freq,samp_rate,mjd,freqmin=-1,freqmax=-1;
  struct timeval start,end;
  char tbuf[30],nfd[32],header[256]="";
  int sign=1;

  // variables to be set by po
  std::string zmq_address, path, output;
  uint16_t port;
  uint32_t channel;
  int hwm;
  size_t num_requested_samples;
  double total_time;

  // setup the program options
  po::options_description desc("Allowed options");
  // clang-format off

  desc.add_options()
      ("help", "help message")
      ("nsamps", po::value<size_t>(&num_requested_samples)->default_value(0), "total number of samples to receive")
      ("duration", po::value<double>(&total_time)->default_value(0), "total number of seconds to receive")
      ("channel", po::value<uint32_t>(&channel)->default_value(0), "VRT channel")
      ("path", po::value<std::string>(&path)->default_value("."), "Output path")
      ("output", po::value<std::string>(&output), "Output filename [default: YYYY-MM-DDTHH:MM:SS.sss_XXXXXX.bin]")
      ("chan-size", po::value<float>(&fchan)->default_value(100), "Channel size [Hz]")
      ("t-int", po::value<float>(&tint)->default_value(1.0), "Integration time [sec]")
      ("n-sub", po::value<int>(&nsub)->default_value(60), "Number of integrations per file")
      ("use", po::value<int>(&nuse)->default_value(1), "Use every n-th integration")
      ("freq-min", po::value<double>(&freqmin), "Frequency range to store (Hz)")
      ("freq-max", po::value<double>(&freqmax), "Frequency range to store (Hz)")
      ("progress", "periodically display short-term bandwidth")
      ("int-second", "align start of reception to integer second")
      ("quiet", "Quiet mode, no output")
      ("continue", "don't abort on a bad packet")
      ("address", po::value<std::string>(&zmq_address)->default_value("localhost"), "VRT ZMQ address")
      ("port", po::value<uint16_t>(&port)->default_value(50100), "VRT ZMQ port")
      ("hwm", po::value<int>(&hwm)->default_value(10000), "VRT ZMQ HWM")
  ;
  // clang-format on
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  // print the help message
  if (vm.count("help")) {
      std::cout << boost::format("VRT samples to STRF rffft. %s") % desc << std::endl;
      std::cout << std::endl
                << "This application streams data from a VRT stream "
                   "to rffft.\n"
                << std::endl;
      return ~0;
  }

  bool progress               = vm.count("progress") > 0;
  bool continue_on_bad_packet = vm.count("continue") > 0;
  bool int_second             = vm.count("int-second") > 0;
  bool useoutput              = vm.count("output") > 0;
  bool quiet                  = vm.count("quiet") > 0;

  context_type difi_context;
  init_context(&difi_context);

  difi_packet_type difi_packet;

  difi_packet.channel_filt = 1<<channel;

  // ZMQ
  void *context = zmq_ctx_new();
  void *subscriber = zmq_socket(context, ZMQ_SUB);
  int rc = zmq_setsockopt (subscriber, ZMQ_RCVHWM, &hwm, sizeof hwm);
  std::string connect_string = "tcp://" + zmq_address + ":" + std::to_string(port);
  rc = zmq_connect(subscriber, connect_string.c_str());
  assert(rc == 0);
  zmq_setsockopt(subscriber, ZMQ_SUBSCRIBE, "", 0);

  // time keeping
  auto start_time = std::chrono::steady_clock::now();
  auto stop_time = start_time + std::chrono::milliseconds(int64_t(1000 * total_time));

  uint32_t buffer[ZMQ_BUFFER_SIZE];
  
  unsigned long long num_total_samps = 0;

  // Track time and samps between updating the BW summary
  auto last_update                     = start_time;
  unsigned long long last_update_samps = 0;

  bool first_frame = true;
  bool start_rx = false;
  uint64_t last_fractional_seconds_timestamp = 0;

  uint32_t signal_pointer = 0;

  // STRF
  uint32_t nint_counter = 0;
  uint32_t nsub_counter = m; // starting index

  while (not stop_signal_called
         and (num_requested_samples > num_total_samps or num_requested_samples == 0)
         and (total_time == 0.0 or std::chrono::steady_clock::now() <= stop_time)) {

      int len = zmq_recv(subscriber, buffer, ZMQ_BUFFER_SIZE, 0);

      const auto now = std::chrono::steady_clock::now();

      if (not difi_process(buffer, sizeof(buffer), &difi_context, &difi_packet)) {
          printf("Not a Vita49 packet?\n");
          continue;
      }

      if (not start_rx and difi_packet.context) {
          difi_print_context(&difi_context);
          start_rx = true;

          samp_rate = difi_context.sample_rate;
          freq = difi_context.rf_freq;

          // Ensure integer number of spectra per subintegration
          tint=ceil(fchan*tint)/fchan;
          
          // Number of channels
          nchan=(int) (samp_rate/fchan);

          // Number of integrations
          nint=(int) (tint*(float) samp_rate/(float) nchan);

          // Get channel range
          if (freqmin>0.0 && freqmax>0.0) {
            imin=(int) ((freqmin-freq+0.5*samp_rate)/fchan);
            imax=(int) ((freqmax-freq+0.5*samp_rate)/fchan);
            if (imin<0 || imin>=nchan || imax<0 || imax>=nchan || imax<=imin) {
              fprintf(stderr,"Output frequency range (%.3lf MHz -> %.3lf MHz) incompatible with\ninput settings (%.3lf MHz center frequency, %.3lf MHz sample rate)!\n",freqmin*1e-6,freqmax*1e-6,freq*1e-6,samp_rate*1e-6);
              return -1;
            }
            partial=1;
          }
          
          // Dump statistics
          printf(" Frequency: %f MHz\n",freq*1e-6);
          printf(" Bandwidth: %f MHz\n",samp_rate*1e-6);
          printf(" Sampling time: %f us\n",1e6/samp_rate);
          printf(" Number of channels: %d\n",nchan);
          printf(" Channel size: %.2f Hz\n",samp_rate/(float) nchan);
          printf(" Integration time: %.2f s\n",tint);
          printf(" Number of averaged spectra: %d\n",nint);
          printf(" Number of subints per file: %d\n",nsub);
          printf(" Starting index: %d\n",m);

          // Allocate
          c=(fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex)*nchan);
          d=(fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex)*nchan);
          cbuf=(char *) malloc(sizeof(char)*2*nchan);
          fbuf=(float *) malloc(sizeof(float)*2*nchan);
          z=(float *) malloc(sizeof(float)*nchan);
          cz=(char *) malloc(sizeof(char)*nchan);
          zw=(float *) malloc(sizeof(float)*nchan);

          // Compute window
          for (i=0;i<nchan;i++)
            zw[i]=0.54-0.46*cos(2.0*M_PI*i/(nchan-1)); 

          // Plan
          fft=fftwf_plan_dft_1d(nchan,c,d,FFTW_FORWARD,FFTW_ESTIMATE);

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
                  stop_time = start_time + std::chrono::milliseconds(int64_t(1000 * total_time));
              }
          }

          if (start_rx and first_frame) {
              std::cout << boost::format(
                               "# First frame: %u samples, %u full secs, %.09f frac secs")
                               % difi_packet.num_rx_samps
                               % difi_packet.integer_seconds_timestamp
                               % ((double)difi_packet.fractional_seconds_timestamp/1e12)
                        << std::endl;
              first_frame = false;
              // STRF Create prefix
              start.tv_sec = difi_packet.integer_seconds_timestamp;
              strftime(prefix,30,"%Y-%m-%dT%T",gmtime(&start.tv_sec));
              
              // File name
              if (not useoutput) {
                sprintf(outfname,"%s/%s_%06d.bin",path.c_str(),prefix,m);
              } else {
                sprintf(outfname,"%s/%s_%06d.bin",path.c_str(),output.c_str(),m);
              }
              outfile=fopen(outfname,"w");
          }

          // int mult = 1;
          for (uint32_t i = 0; i < difi_packet.num_rx_samps; i++) {

              if (signal_pointer==0 and nint_counter==0) {
                // gettimeofday(&start,0);
                uint64_t seconds = difi_packet.integer_seconds_timestamp;
                uint64_t frac_seconds = difi_packet.fractional_seconds_timestamp;
                frac_seconds += (i+1)*1e12/difi_context.sample_rate;
                if (frac_seconds > 1e12) {
                    frac_seconds -= 1e12;
                    seconds++;
                }
                start.tv_sec = seconds;
                start.tv_usec = frac_seconds/1e6;
              }

              int16_t re;
              memcpy(&re, (char*)&buffer[difi_packet.offset+i], 2);
              int16_t img;
              memcpy(&img, (char*)&buffer[difi_packet.offset+i]+2, 2);
              c[signal_pointer][REAL]=(float)(re/32768.0)*zw[signal_pointer];
              c[signal_pointer][IMAG]=(float)(img/32768.0)*zw[signal_pointer]*sign;
              // mult *= -1;

              signal_pointer++;

              if (signal_pointer >= nchan) { 

                  signal_pointer = 0;

                  // Execute
                  fftwf_execute(fft);
  
                  // Shift/Integrate
                  for (j=0;j<nchan;j++) {
                    if (j<nchan/2)
                      l=j+nchan/2;
                    else
                      l=j-nchan/2;

                    z[l]+=d[j][0]*d[j][0]+d[j][1]*d[j][1];
                  }
                  nint_counter++;

                  if (nint_counter >= nint) {
                    // Log end time
                    // gettimeofday(&end,0);
                    uint64_t seconds = difi_packet.integer_seconds_timestamp;
                    uint64_t frac_seconds = difi_packet.fractional_seconds_timestamp;
                    frac_seconds += (i+1)*1e12/difi_context.sample_rate;
                    if (frac_seconds > 1e12) {
                        frac_seconds -= 1e12;
                        seconds++;
                    }
                    end.tv_sec = seconds;
                    end.tv_usec = frac_seconds/1e6;

                    // Process nint block
                    // Time stats
                    length=(end.tv_sec-start.tv_sec)+(end.tv_usec-start.tv_usec)*1e-6;
                    
                    // Scale
                    for (i=0;i<nchan;i++) 
                      z[i] *= (float)nuse/(float)nchan;
                    
                    // Format start time
                    strftime(tbuf,30,"%Y-%m-%dT%T",gmtime(&start.tv_sec));
                    sprintf(nfd,"%s.%03ld",tbuf,start.tv_usec/1000);

                    // Header
                    if (partial==0) {
                      if (outformat=='f') 
                        sprintf(header,"HEADER\nUTC_START    %s\nFREQ         %lf Hz\nBW           %lf Hz\nLENGTH       %f s\nNCHAN        %d\nNSUB         %d\nEND\n",nfd,freq,samp_rate,length,nchan,nsub);
                      else if (outformat=='c')
                        sprintf(header,"HEADER\nUTC_START    %s\nFREQ         %lf Hz\nBW           %lf Hz\nLENGTH       %f s\nNCHAN        %d\nNSUB         %d\nNBITS         8\nMEAN         %e\nRMS          %e\nEND\n",nfd,freq,samp_rate,length,nchan,nsub,zavg,zstd);
                          } else if (partial==1) {
                      if (outformat=='f') 
                        sprintf(header,"HEADER\nUTC_START    %s\nFREQ         %lf Hz\nBW           %lf Hz\nLENGTH       %f s\nNCHAN        %d\nNSUB         %d\nEND\n",nfd,0.5*(freqmax+freqmin),freqmax-freqmin,length,imax-imin,nsub);
                      else if (outformat=='c')
                        sprintf(header,"HEADER\nUTC_START    %s\nFREQ         %lf Hz\nBW           %lf Hz\nLENGTH       %f s\nNCHAN        %d\nNSUB         %d\nNBITS         8\nMEAN         %e\nRMS          %e\nEND\n",nfd,0.5*(freqmax+freqmin),freqmax-freqmin,length,imax-imin,nsub,zavg,zstd);
                    }
                    // Limit output
                    if (!quiet)
                      printf("%s %s %f %d\n",outfname,nfd,length,nint_counter);

                    // Dump file
                    fwrite(header,sizeof(char),256,outfile);
                    if (partial==0) {
                      if (outformat=='f')
                        fwrite(z,sizeof(float),nchan,outfile);
                      else if (outformat=='c')
                        fwrite(cz,sizeof(char),nchan,outfile);
                    } else if (partial==1) {
                      if (outformat=='f')
                        fwrite(&z[imin],sizeof(float),imax-imin,outfile);
                      else if (outformat=='c')
                        fwrite(&cz[imin],sizeof(char),imax-imin,outfile);
                    }

                    nsub_counter++;

                    if (nsub_counter >= nsub) {
                      fclose(outfile);
                      m++;
                      if (not useoutput) {
                        sprintf(outfname,"%s/%s_%06d.bin",path.c_str(),prefix,m);
                      } else {
                        sprintf(outfname,"%s/%s_%06d.bin",path.c_str(),output.c_str(),m);
                      }
                      outfile=fopen(outfname,"w");
                      nsub_counter=0;
                    }

                    // clear z
                    for (i=0;i<nchan;i++) 
                      z[i]=0.0;

                    // reset counter
                    nint_counter = 0;
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

  // Close file
  fclose(outfile);
  
  // Destroy plan
  fftwf_destroy_plan(fft);

  // Deallocate
  free(cbuf);
  free(fbuf);
  fftwf_free(c);
  fftwf_free(d);
  free(z);
  free(cz);
  free(zw);
  
  return 0;
}

