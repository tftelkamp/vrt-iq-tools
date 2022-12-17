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

// DADA
#include <sstream>
#include <dada_def.h>
#include <dada_hdu.h>
#include <multilog.h>
#include <ipcio.h>
#include <iomanip>
#include <ascii_header.h>
// #include "interleave.h"
// END DADA

#include "difi-tools.h"

namespace po = boost::program_options;

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

    // variables to be set by po
    std::string file, type, zmq_address;
    uint16_t port;
    uint32_t channel;
    int16_t scale;
    int hwm;
    size_t num_requested_samples;
    double total_time;

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
        ("channel", po::value<uint32_t>(&channel)->default_value(0), "DIFI channel")
        ("int-second", "align start of reception to integer second")
        ("scale", po::value<int16_t>(&scale)->default_value(128), "scaling factor for 16 to 8 bit conversion (default 128)")
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
        std::cout << boost::format("DIFI samples to nothing. %s") % desc << std::endl;
        std::cout << std::endl
                  << "This application streams data from a DIFI stream "
                     "to nowhwere.\n"
                  << std::endl;
        return ~0;
    }

    bool progress               = vm.count("progress") > 0;
    bool stats                  = vm.count("stats") > 0;
    bool null                   = vm.count("null") > 0;
    bool continue_on_bad_packet = vm.count("continue") > 0;
    bool int_second             = (bool)vm.count("int-second");

    if (!vm.count("int-second")) throw std::runtime_error("Dada requires --int-second");

    context_type difi_context;
    init_context(&difi_context);

    difi_packet_type difi_packet;

    difi_packet.channel_filt = 1<<channel;

    // DADA
    dada_hdu_t *dada_hdu;
    multilog_t *dada_log;
    std::string dada_header;
    std::complex<int8_t> dadabuffer[DIFI_SAMPLES_PER_PACKET] __attribute((aligned(32)));

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
    uint64_t last_fractional_seconds_timestamp = 0;

    // set to true to process data before context
    bool start_rx = false;

    uint32_t signal_pointer = 0;

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
            // Possibly do something with context here
            // DADA
            // if (wirefmt != "sc8") throw std::runtime_error("Dada requires --wirefmt==sc8");
            // if (type != "char") throw std::runtime_error("Dada requires --type==char");
            // TODO dada multiple polarizations!
            dada_header = 
              "HEADER DADA\n"
              "HDR_VERSION 1.0\n"
              "HDR_SIZE    4096\n"
              "FREQ " + std::to_string(difi_context.rf_freq/1e6) + "\n"
              "BW " + std::to_string(int(difi_context.sample_rate/1e6)) + "\n"
              "TELESCOPE DWL\n"
              "RECEIVER DWL\n"
              "INSTRUMENT DWL\n"
              "SOURCE UNDEFINED\n"
              "NBIT " + "8\n" +
              "NDIM " + "1\n" +
              "NPOL " + "1\n" +
              "RESOLUTION 1\n"
              "OBS_OFFSET 0\n"
              "TSAMP " + std::to_string(1e6/difi_context.sample_rate) + "\n";

            // DADA hdu
            dada_log = multilog_open ("example_dada_writer", 0);
            multilog_add(dada_log, stderr);
            dada_hdu = dada_hdu_create(dada_log);
            dada_hdu_set_key(dada_hdu, 0xc2c2);
            if (dada_hdu_connect (dada_hdu) < 0)
                throw std::runtime_error("Could not connect to DADA HDU");
            if (dada_hdu_lock_write(dada_hdu) < 0 )
                throw std::runtime_error("Could not get write lock on DADA HDU");
            // END DADA
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

            // Process data here
            // Assumes ci16_le

            for (uint32_t i = 0; i < difi_packet.num_rx_samps; i++) {
                auto sample = (std::complex<int16_t>)buffer[difi_packet.offset+i];
                // Convert ci16_le to cs8
                dadabuffer[i] = sample/scale;
            }
            if (ipcio_write(dada_hdu->data_block, (char*)dadabuffer, difi_packet.num_rx_samps*sizeof(std::complex<int8_t>)) < 0)
               throw std::runtime_error("Error writing buffer to DADA");

            // data: (const char*)&buffer[difi_packet.offset]
            // size (bytes): sizeof(uint32_t)*difi_packet.num_rx_samps
             
            num_total_samps += difi_packet.num_rx_samps;

            if (start_rx and first_frame) {
                std::cout << boost::format(
                                 "  First frame: %u samples, %u full secs, %.09f frac secs")
                                 % difi_packet.num_rx_samps
                                 % difi_packet.integer_seconds_timestamp
                                 % ((double)difi_packet.fractional_seconds_timestamp/1e12)
                          << std::endl;
                first_frame = false;
                // DADA
                // Put starttime into dada_header
                std::ostringstream starttime_str;
                std::time_t starttime_time_t = difi_packet.integer_seconds_timestamp;
                std::tm *starttime_tm = std::gmtime(&starttime_time_t);
                starttime_str << std::put_time(starttime_tm, "%Y-%m-%d-%H:%M:%S");
                dada_header.append("UTC_START " + starttime_str.str() + "\n");
                dada_header.resize(4096, ' ');
                // Write dada header to dada_header.txt for debugging
                std::ofstream dada_header_txt("dada_header.txt");
                dada_header_txt << dada_header;
                dada_header_txt.close();
                {
                   // Write dada header to buffer
                   uint64_t header_size = ipcbuf_get_bufsz (dada_hdu->header_block);
                   char * ipc_header = ipcbuf_get_next_write (dada_hdu->header_block);
                   strncpy(ipc_header, dada_header.c_str(), header_size);

                   if (ipcbuf_mark_filled (dada_hdu->header_block, 4096) < 0)
                       throw std::runtime_error("Could not mark filled Header Block");
                }
            }
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

    // DADA
    // unlock write access from the HDU, performs implicit EOD
    if (dada_hdu_unlock_write (dada_hdu) < 0)
        throw std::runtime_error("dada_hdu_unlock_write failed");

    // disconnect from HDU
    if (dada_hdu_disconnect (dada_hdu) < 0)
        throw std::runtime_error("could not unlock write on DADA hdu");

    zmq_close(subscriber);
    zmq_ctx_destroy(context);

    return 0;

}  
