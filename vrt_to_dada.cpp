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
#include <csignal>
#include <fstream>
#include <iostream>
#include <thread>
#include <complex>

// VRT
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include <vrt/vrt_read.h>
#include <vrt/vrt_string.h>
#include <vrt/vrt_types.h>
#include <vrt/vrt_util.h>

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

#include "vrt-tools.h"
#include "dt-extended-context.h"

namespace po = boost::program_options;

static bool stop_signal_called = false;
void sig_int_handler(int)
{
    std::cout<<"Stop signal caught\n";
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
    std::string zmq_address, channel_list, sourcename, dadakey_str;
    uint16_t instance, main_port, port;
    uint32_t channel;
    int hwm;
    size_t num_requested_samples;
    double total_time;
    float amplitude, phase;

    // setup the program options
    po::options_description desc("Allowed options");
    // clang-format off

    desc.add_options()
        ("help", "help message")
        ("nsamps", po::value<size_t>(&num_requested_samples)->default_value(0), "total number of samples to receive")
        ("duration", po::value<double>(&total_time)->default_value(0), "total number of seconds to receive")
        ("amplitude", po::value<float>(&amplitude)->default_value(1), "amplitude correction of second channel")
        ("phase", po::value<float>(&phase)->default_value(0), "phase shift on second channel [0-1]")
        ("sourcename", po::value<std::string>(&sourcename)->default_value("undefined"), "name of tracked celestial source")
        ("key", po::value<std::string>(&dadakey_str)->default_value("c2c2"), "dada key")
        ("progress", "periodically display short-term bandwidth")
        ("channel", po::value<std::string>(&channel_list)->default_value("0"), "which VRT channel(s) to use (specify \"0\", \"1\", \"0,1\", etc)")
        ("continue", "don't abort on a bad packet")
        ("dt-trace", "add DT trace data")
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
        std::cout << boost::format("VRT samples to dada. %s") % desc << std::endl;
        std::cout << std::endl
                  << "This application streams data from a VRT stream "
                     "to a psrdada ring buffer.\n"
                  << std::endl;
        return ~0;
    }

    bool progress               = vm.count("progress") > 0;
    bool stats                  = vm.count("stats") > 0;
    bool continue_on_bad_packet = vm.count("continue") > 0;
    bool dt_trace               = vm.count("dt-trace") > 0;
    bool zmq_split              = vm.count("zmq-split") > 0;

    std::complex<float> z(0,-2*(float)M_PI*phase);
    std::complex<float> a(amplitude,0);
    std::complex<float> correction = a*exp(z);

    context_type vrt_context;
    dt_ext_context_type dt_ext_context;
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

    if (zmq_split) {
        if (channel_nums.size()>1) {
            printf("Multiple channels with --zmq-split is not supported.\n");
            exit(EXIT_FAILURE);
        }
        main_port += channel_nums[0];
        channel_nums[0] = 0;
        vrt_packet.channel_filt = 1;
    }

    // DADA
    dada_hdu_t *dada_hdu;
    multilog_t *dada_log;
    std::string dada_header;
    key_t dadakey;

    dadakey = std::stoul(dadakey_str, nullptr, 16);

    std::complex<float> dadabuffer[VRT_SAMPLES_PER_PACKET*MAX_CHANNELS] __attribute((aligned(32)));

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

    std::signal(SIGINT, &sig_int_handler);
    if (num_requested_samples == 0 and total_time == 0) {
        std::cout << "Press Ctrl + C to stop streaming..." << std::endl;
    }

    while (not stop_signal_called
           and (num_requested_samples*channel_nums.size() > num_total_samps or num_requested_samples == 0)) {

        int len = zmq_recv(subscriber, buffer, ZMQ_BUFFER_SIZE, 0);

        const auto now = std::chrono::steady_clock::now();

        if (not vrt_process(buffer, sizeof(buffer), &vrt_context, &vrt_packet)) {
            printf("Not a Vita49 packet?\n");
            continue;
        }

        uint32_t ch = 0;
        for(ch = 0; ch<channel_nums.size(); ch++)
            if (vrt_packet.stream_id & (1 << channel_nums[ch]) )
                break;

        uint32_t channel = channel_nums[ch];

        if (vrt_packet.extended_context) {
            // TODO: Find some other variable to avoid giving this warning for every extended context packet
            // This now assumes that any extended context is a DT extended context
            if (not dt_ext_context.dt_ext_context_received and not dt_trace) {
                std::cerr << "# WARNING: DT metadata is present in the stream, but it is ignored. Did you forget --dt-trace?" << std::endl;
            }
            dt_process(buffer, sizeof(buffer), &vrt_packet, &dt_ext_context);
        }

        if (not start_rx and vrt_packet.context and (dt_ext_context.dt_ext_context_received or not dt_trace)) {
            vrt_print_context(&vrt_context);
            start_rx = true;

            if (total_time > 0)
                num_requested_samples = total_time * vrt_context.sample_rate;

            // Possibly do something with context here
            // DADA
            dada_header =
              "HEADER DADA\n"
              "HDR_VERSION 1.0\n"
              "HDR_SIZE    4096\n"
              "FREQ " + std::to_string(vrt_context.rf_freq/1e6) + "\n"
              "BW " + std::to_string(int(vrt_context.sample_rate/1e6)) + "\n"
              "TELESCOPE DWL\n"
              "RECEIVER VRT\n"
              "INSTRUMENT dspsr\n"
              "SOURCE " + sourcename + "\n"
              "NBIT " + "32\n"
              "NDIM " + "2\n"
              "NPOL " + std::to_string(channel_nums.size()) + "\n"
              "RESOLUTION 1\n"
              "OBS_OFFSET 0\n"
              "TSAMP " + std::to_string(1e6/vrt_context.sample_rate) + "\n";

            if (dt_trace) {
              boost::format fmt("%02d:%02d:%06.3f");

              double ra_h = ((12.0/M_PI) * dt_ext_context.ra_current);
              int ra_hours = static_cast<int>(ra_h);
              int ra_minutes = static_cast<int>(ra_h * 60) % 60;
              double ra_seconds = fmod(ra_h * 3600.0, 60.0);
              fmt % ra_hours % ra_minutes % ra_seconds;
              dada_header += "RA " + boost::str(fmt) + "\n";

              double dec_deg = ((180.0/M_PI) * dt_ext_context.dec_current);
              std::string dec_sign = (dec_deg > 0 ? "+" : "-");
              dec_deg = abs(dec_deg);
              int dec_degrees = static_cast<int>(dec_deg);
              int dec_minutes = static_cast<int>(dec_deg * 60.0) % 60;
              double dec_seconds = fmod(dec_deg * 3600, 60.0);
              fmt % dec_degrees % dec_minutes % dec_seconds;
              dada_header += "DEC " + dec_sign + boost::str(fmt) + "\n";


            }

            // DADA hdu
            dada_log = multilog_open ("example_dada_writer", 0);
            multilog_add(dada_log, stderr);
            dada_hdu = dada_hdu_create(dada_log);
            dada_hdu_set_key(dada_hdu, dadakey);
            if (dada_hdu_connect (dada_hdu) < 0)
                throw std::runtime_error("Could not connect to DADA HDU");
            if (dada_hdu_lock_write(dada_hdu) < 0 )
                throw std::runtime_error("Could not get write lock on DADA HDU");
            // END DADA
        }

        if (start_rx and vrt_packet.data and (dt_ext_context.dt_ext_context_received or not dt_trace)) {

            if (vrt_packet.lost_frame)
               if (not continue_on_bad_packet)
                    break;

            // check if fractional second has wrapped
            if (vrt_packet.fractional_seconds_timestamp > last_fractional_seconds_timestamp ) {
                    last_fractional_seconds_timestamp = vrt_packet.fractional_seconds_timestamp;
                    continue;
            } else {
                last_update = now;
                start_time = now;
                stop_time = start_time + std::chrono::milliseconds(int64_t(1000 * total_time));
            }

            // Process data here
            // Assumes ci16_le

            for (uint32_t i = 0; i < vrt_packet.num_rx_samps; i++) {
                int16_t re;
                memcpy(&re, (char*)&buffer[vrt_packet.offset+i], 2);
                int16_t img;
                memcpy(&img, (char*)&buffer[vrt_packet.offset+i]+2, 2);
                // Convert ci16_le to float
                std::complex<float>sample(re,img);
                if (channel_nums.size() > 1) {
                    if (ch==1)
                        dadabuffer[i*channel_nums.size()+ch] = correction*sample;
                    else
                        dadabuffer[i*channel_nums.size()+ch] = sample;
                } else {
                    dadabuffer[i] = sample;
                }
            }

            // send when all channels have been received
            if (ch == channel_nums.size()-1) {
                if (ipcio_write(dada_hdu->data_block, (char*)dadabuffer, channel_nums.size()*vrt_packet.num_rx_samps*sizeof(std::complex<float>)) < 0) {
                    if (stop_signal_called) {
                        break;
                    }
                    throw std::runtime_error("Error writing buffer to DADA");
                }
            }

            // data: (const char*)&buffer[vrt_packet.offset]
            // size (bytes): sizeof(uint32_t)*vrt_packet.num_rx_samps

            num_total_samps += vrt_packet.num_rx_samps;

            if (start_rx and first_frame) {
                std::cout << boost::format(
                                 "  First frame: %u samples, %u full secs, %.09f frac secs")
                                 % vrt_packet.num_rx_samps
                                 % vrt_packet.integer_seconds_timestamp
                                 % ((double)vrt_packet.fractional_seconds_timestamp/1e12)
                          << std::endl;
                first_frame = false;
                // DADA
                // Put starttime into dada_header
                std::ostringstream starttime_str;
                std::time_t starttime_time_t = vrt_packet.integer_seconds_timestamp;
                std::tm *starttime_tm = std::gmtime(&starttime_time_t);
                starttime_str << std::put_time(starttime_tm, "%Y-%m-%d-%H:%M:%S");
                dada_header.append("UTC_START " + starttime_str.str() + "\n");
                dada_header.append("PICOSECONDS " + std::to_string(vrt_packet.fractional_seconds_timestamp) + "\n");
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

    // DADA
    // unlock write access from the HDU, performs implicit EOD
    if (dada_hdu_unlock_write (dada_hdu) < 0)
        throw std::runtime_error("dada_hdu_unlock_write failed");

    // disconnect from HDU
    if (dada_hdu_disconnect (dada_hdu) < 0)
        throw std::runtime_error("could not unlock write on DADA hdu");

    zmq_close(subscriber);
    zmq_ctx_destroy(context);
    std::cout<<"vrt_to_dada cleaned up properly after SIGINT\n";

    return 0;
}
