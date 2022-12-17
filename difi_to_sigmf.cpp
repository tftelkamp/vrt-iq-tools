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

#include "difi-tools.h"

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
        std::cout << boost::format("DIFI samples to file %s") % desc << std::endl;
        std::cout << std::endl
                  << "This application streams data from a DIFI stream "
                     "to a file.\n"
                  << std::endl;
        return ~0;
    }

    bool progress               = vm.count("progress") > 0;
    bool stats                  = vm.count("stats") > 0;
    bool null                   = vm.count("null") > 0;
    bool continue_on_bad_packet = vm.count("continue") > 0;
    bool int_second             = vm.count("int-second");

    context_type difi_context;
    init_context(&difi_context);

    difi_packet_type difi_packet;
    difi_packet.channel_filt = 0;

    // detect which channels to use
    std::vector<std::string> channel_strings;
    std::vector<size_t> channel_nums;
    boost::split(channel_strings, channel_list, boost::is_any_of("\"',"));
    for (size_t ch = 0; ch < channel_strings.size(); ch++) {
        size_t chan = std::stoi(channel_strings[ch]);
        channel_nums.push_back(std::stoi(channel_strings[ch]));
        difi_packet.channel_filt |= 1<<std::stoi(channel_strings[ch]);
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

        const auto now = std::chrono::steady_clock::now();

        if (not difi_process(buffer, sizeof(buffer), &difi_context, &difi_packet)) {
            printf("Not a Vita49 packet?\n");
            continue;
        }

        uint32_t ch = 0;
        while(not (difi_packet.stream_id & (1 << channel_nums[ch]) ) )
            ch++;

        uint32_t channel = channel_nums[ch];

        if ( not (context_recv & difi_packet.stream_id) and difi_packet.context and not first_frame) {

            context_recv |= difi_packet.stream_id;

            difi_print_context(&difi_context);

            if (not null) {
                // std::cout << "Writing SigMF metadata..." << std::endl;
                if (!mdfile) {
                    std::cout << "File not created?!";
                } else {
                    *metafiles[ch] << boost::format("{ \n"
                    "    \"global\": {\n"
                    "        \"core:version\": \"1.0.0\",\n"
                    "        \"core:recorder\": \"difi_to_sigmf\",\n"
                    "        \"core:sample_rate\": %u,\n"
                    "        \"core:datatype\": \"ci16_le\",\n"
                    "        \"difi:rx_gain\": %i,\n"
                    "        \"difi:bandwidth\": %u,\n"
                    "        \"difi:reference\": \"%s\",\n"
                    "        \"difi:time_source\": \"%s\",\n"
                    "        \"difi:stream_id\": %u,\n"
                    "        \"difi:channel\": %u\n"
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
                    % difi_context.sample_rate
                    % difi_context.gain
                    % difi_context.bandwidth
                    % (difi_context.reflock ? "external" : "internal")
                    % (difi_context.time_cal ? "pps" : "internal")
                    % difi_context.stream_id
                    % channel
                    % difi_context.rf_freq
                    % (boost::posix_time::to_iso_extended_string(boost::posix_time::from_time_t(difi_context.starttime_integer)))
                    % (double)(difi_context.starttime_fractional/1e6);
                    *metafiles[ch] << std::endl;
                    metafiles[ch]->close();
                }
            }

            if (total_time > 0)  
                num_requested_samples = total_time * difi_context.sample_rate;
        }

        if (difi_packet.data) {

            if (difi_packet.lost_frame)
               if (not continue_on_bad_packet)
                    break;

            if (int_second) {
                // check if fractional second has wrapped
                if (difi_packet.fractional_seconds_timestamp >= last_fractional_seconds_timestamp ) {
                        last_fractional_seconds_timestamp = difi_packet.fractional_seconds_timestamp;
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
                                 % difi_packet.num_rx_samps
                                 % difi_packet.integer_seconds_timestamp
                                 % ((double)difi_packet.fractional_seconds_timestamp/1e12)
                                 % (int32_t)difi_context.last_data_counter
                          << std::endl;
                first_frame = false;
            }

            // Write to file
            outfiles[ch]->write(
                (const char*)&buffer[difi_packet.offset], sizeof(uint32_t)*difi_packet.num_rx_samps);

            num_total_samps += difi_packet.num_rx_samps;
        }

        if (progress and not int_second) {
            if (difi_packet.data)
                last_update_samps += difi_packet.num_rx_samps;
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

    for (size_t i = 0; i < outfiles.size(); i++)
        outfiles[i]->close();
        
    zmq_close(subscriber);
    zmq_ctx_destroy(context);

    return 0;

}  
