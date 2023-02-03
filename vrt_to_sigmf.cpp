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
#include <complex>
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

#include "vrt-tools.h"

namespace po = boost::program_options;

static bool stop_signal_called = false;
void sig_int_handler(int)
{
    stop_signal_called = true;
}

//! Change to filename, e.g. from usrp_samples.dat to usrp_samples.chan0.dat,
//  but only if multiple names are to be generated.
std::string generate_out_filename(
    const std::string& base_fn, size_t n_names, size_t this_name)
{
    if (n_names == 1) {
        return base_fn;
    }

    boost::filesystem::path base_fn_fp(base_fn);
    base_fn_fp.replace_extension(boost::filesystem::path(
        str(boost::format("chan%d%s") % this_name % base_fn_fp.extension().string())));
    return base_fn_fp.string();
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
    std::string file, type, zmq_address, channel_list;
    size_t num_requested_samples, total_time;
    uint16_t port;
    int hwm;

    // setup the program options
    po::options_description desc("Allowed options");
    // clang-format off

    desc.add_options()
        ("help", "help message")
        ("file", po::value<std::string>(&file)->default_value("usrp_samples.dat"), "name of the file to write binary samples to")
        // ("type", po::value<std::string>(&type)->default_value("short"), "sample type: double, float, or short")
        ("nsamps", po::value<size_t>(&num_requested_samples)->default_value(0), "total number of samples to receive")
        ("duration", po::value<size_t>(&total_time)->default_value(0), "total number of seconds to receive")
        ("progress", "periodically display short-term bandwidth")
        ("channel", po::value<std::string>(&channel_list)->default_value("0"), "which channel(s) to use (specify \"0\", \"1\", \"0,1\", etc)")
        // ("stats", "show average bandwidth on exit")
        ("int-second", "align start of reception to integer second")
        ("null", "run without writing to file")
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
        std::cout << boost::format("VRT samples to file %s") % desc << std::endl;
        std::cout << std::endl
                  << "This application streams data from a VRT stream "
                     "to a file.\n"
                  << std::endl;
        return ~0;
    }

    bool progress               = vm.count("progress") > 0;
    bool stats                  = vm.count("stats") > 0;
    bool null                   = vm.count("null") > 0;
    bool continue_on_bad_packet = vm.count("continue") > 0;
    bool int_second             = vm.count("int-second");

    context_type vrt_context;
    init_context(&vrt_context);

    packet_type vrt_packet;
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

    std::vector<std::shared_ptr<std::ofstream>> outfiles;
    std::vector<std::shared_ptr<std::ofstream>> metafiles;

    std::string mdfilename;
    std::ofstream mdfile;
    mdfilename = file + ".sigmf-meta";
    file = file + ".sigmf-data";

    if (not null) 
        for (size_t i = 0; i < channel_nums.size(); i++) {
            const std::string this_filename = generate_out_filename(file, channel_nums.size(), channel_nums[i]);
            outfiles.push_back(std::shared_ptr<std::ofstream>(
                new std::ofstream(this_filename.c_str(), std::ofstream::binary)));
            
            const std::string meta_filename = generate_out_filename(mdfilename, channel_nums.size(), channel_nums[i]);
            metafiles.push_back(std::shared_ptr<std::ofstream>(
                new std::ofstream(meta_filename.c_str())));
    }

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

    // ZMQ buffer
    uint32_t buffer[ZMQ_BUFFER_SIZE];

    unsigned long long num_total_samps = 0;

    // Track time and samps between updating the BW summary
    auto last_update                     = start_time;
    unsigned long long last_update_samps = 0;

    uint64_t last_fractional_seconds_timestamp = 0;

    bool first_frame = true;
    uint32_t context_recv = 0;

    if (num_requested_samples == 0) {
        std::signal(SIGINT, &sig_int_handler);
        std::cout << "Press Ctrl + C to stop receiving..." << std::endl;
    }

    while (not stop_signal_called
           and ( num_requested_samples*channel_nums.size() > num_total_samps or num_requested_samples == 0)) {

        int len = zmq_recv(subscriber, buffer, ZMQ_BUFFER_SIZE, 0);

        if (stop_signal_called)
            break;

        const auto now = std::chrono::steady_clock::now();

        if (not vrt_process(buffer, sizeof(buffer), &vrt_context, &vrt_packet)) {
            printf("Not a Vita49 packet?\n");
            continue;
        }

        if (vrt_packet.context and not first_frame and not continue_on_bad_packet and vrt_context.context_changed) {
            printf("Context changed, exiting.\n");
            break;
        }

        uint32_t ch = 0;
        for(ch = 0; ch<channel_nums.size(); ch++)
            if (vrt_packet.stream_id & (1 << channel_nums[ch]) )
                break;

        uint32_t channel = channel_nums[ch];

        if ( not (context_recv & vrt_packet.stream_id) and vrt_packet.context and not first_frame) {

            context_recv |= vrt_packet.stream_id;

            vrt_print_context(&vrt_context);

            if (not null) {
                // std::cout << "Writing SigMF metadata..." << std::endl;
                if (!metafiles[ch]) {
                    std::cout << "File not created?!";
                } else {
                    *metafiles[ch] << boost::format("{ \n"
                    "    \"global\": {\n"
                    "        \"core:version\": \"1.0.0\",\n"
                    "        \"core:recorder\": \"vrt_to_sigmf\",\n"
                    "        \"core:sample_rate\": %u,\n"
                    "        \"core:datatype\": \"ci16_le\",\n"
                    "        \"vrt:rx_gain\": %i,\n"
                    "        \"vrt:bandwidth\": %u,\n"
                    "        \"vrt:reference\": \"%s\",\n"
                    "        \"vrt:time_source\": \"%s\",\n"
                    "        \"vrt:stream_id\": %u,\n"
                    "        \"vrt:channel\": %u\n"
                    "    },\n"
                    "    \"annotations\": [],\n"
                    "    \"captures\": [\n"
                    "        {\n"
                    "            \"core:sample_start\": 0,\n"
                    "            \"core:frequency\": %u,\n"
                    "            \"core:datetime\": \"%s.%06.0f\"\n"
                    "        }\n"
                    "    ]\n"
                    "}\n")
                    % vrt_context.sample_rate
                    % vrt_context.gain
                    % vrt_context.bandwidth
                    % (vrt_context.reflock ? "external" : "internal")
                    % (vrt_context.time_cal ? "pps" : "internal")
                    % vrt_context.stream_id
                    % channel
                    % vrt_context.rf_freq
                    % (boost::posix_time::to_iso_extended_string(boost::posix_time::from_time_t(vrt_context.starttime_integer)))
                    % (double)(vrt_context.starttime_fractional/1e6);
                    *metafiles[ch] << std::endl;
                    metafiles[ch]->close();
                }
            }

            if (total_time > 0)  
                num_requested_samples = total_time * vrt_context.sample_rate;
        }

        if (vrt_packet.data) {

            if (vrt_packet.lost_frame)
               if (not continue_on_bad_packet)
                    break;

            if (int_second) {
                // check if fractional second has wrapped
                if (vrt_packet.fractional_seconds_timestamp >= last_fractional_seconds_timestamp ) {
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
                                 "# First frame: %u samples, %u full secs, %.09f frac secs (counter %i)")
                                 % vrt_packet.num_rx_samps
                                 % vrt_packet.integer_seconds_timestamp
                                 % ((double)vrt_packet.fractional_seconds_timestamp/1e12)
                                 % (int32_t)vrt_context.last_data_counter
                          << std::endl;
                first_frame = false;
                last_update = now;
                // update context starttime in case of int_second
                vrt_context.starttime_integer = vrt_packet.integer_seconds_timestamp;
                vrt_context.starttime_fractional = vrt_packet.fractional_seconds_timestamp;
            }

            // Write to file
            if (not null) {
                outfiles[ch]->write(
                    (const char*)&buffer[vrt_packet.offset], sizeof(uint32_t)*vrt_packet.num_rx_samps);
            }

            num_total_samps += vrt_packet.num_rx_samps;
        }

        if (progress and not int_second) {
            if (vrt_packet.data)
                last_update_samps += vrt_packet.num_rx_samps;
            const auto time_since_last_update = now - last_update;
            if (time_since_last_update > std::chrono::seconds(1)) {
                const double time_since_last_update_s =
                    std::chrono::duration<double>(time_since_last_update).count();
                const double rate = (double(last_update_samps) / time_since_last_update_s) / (double)channel_nums.size();
                std::cout << "\t" << (rate / 1e6) << " Msps, ";
                
                last_update_samps = 0;
                last_update       = now;
    
                float sum_i = 0;
                uint32_t clip_i = 0;

                double datatype_max = 32768.;
                // if (cpu_format == "sc8" || cpu_format == "s8")
                //     datatype_max = 128.;

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

    for (size_t i = 0; i < outfiles.size(); i++)
        outfiles[i]->close();
        
    zmq_close(subscriber);
    zmq_ctx_destroy(context);

    return 0;

}  
