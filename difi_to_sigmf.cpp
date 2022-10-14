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

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

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

namespace po = boost::program_options;

using boost::property_tree::ptree;

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
    std::string file, type, zmq_address;
    size_t num_requested_samples;
    double total_time;

    // setup the program options
    po::options_description desc("Allowed options");
    // clang-format off

    desc.add_options()
        ("help", "help message")
        ("file", po::value<std::string>(&file)->default_value("usrp_samples.dat"), "name of the file to write binary samples to")
        // ("type", po::value<std::string>(&type)->default_value("short"), "sample type: double, float, or short")
        ("nsamps", po::value<size_t>(&num_requested_samples)->default_value(0), "total number of samples to receive")
        ("duration", po::value<double>(&total_time)->default_value(0), "total number of seconds to receive")
        ("progress", "periodically display short-term bandwidth")
        // ("stats", "show average bandwidth on exit")
        ("int-second", "align start of reception to integer second")
        ("null", "run without writing to file")
        ("continue", "don't abort on a bad packet")
        ("address", po::value<std::string>(&zmq_address)->default_value("localhost"), "DIFI ZMQ address")
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

    int64_t rf_freq = 0;
    uint32_t sample_rate = 0;
    uint32_t gain = 0;
    uint32_t bandwidth = 0;
    bool reflock = false;
    bool time_cal = false;
    uint32_t stream_id;

    uint64_t starttime_integer;
    uint64_t starttime_fractional;

    std::vector<std::shared_ptr<std::ofstream>> outfiles;
    std::vector<size_t> channel_nums = {0}; // single channel (0)

    std::string mdfilename;
    mdfilename = file + ".sigmf-meta";
    file = file + ".sigmf-data";

    if (not null)
        for (size_t i = 0; i < channel_nums.size(); i++) {
            const std::string this_filename = generate_out_filename(file, channel_nums.size(), channel_nums[i]);
            outfiles.push_back(std::shared_ptr<std::ofstream>(
                new std::ofstream(this_filename.c_str(), std::ofstream::binary)));
        }

    // ZMQ
    void *context = zmq_ctx_new();
    void *subscriber = zmq_socket(context, ZMQ_SUB);
    int hwm = 10000;
    int rc = zmq_setsockopt (subscriber, ZMQ_RCVHWM, &hwm, sizeof hwm);
    std::string connect_string = "tcp://" + zmq_address + ":50100";
    rc = zmq_connect(subscriber, connect_string.c_str());
    assert(rc == 0);
    zmq_setsockopt(subscriber, ZMQ_SUBSCRIBE, "", 0);

    bool first_frame = true;

    // time keeping
    auto start_time = std::chrono::steady_clock::now();
    auto stop_time =
        start_time + std::chrono::milliseconds(int64_t(1000 * total_time));

    int8_t packet_count = -1;

    uint32_t buffer[100000];

    struct vrt_header h;
    struct vrt_fields f;
    
    unsigned long long num_total_samps = 0;

    // Track time and samps between updating the BW summary
    auto last_update                     = start_time;
    unsigned long long last_update_samps = 0;

    bool start_rx = false;
    uint64_t last_fractional_seconds_timestamp = 0;

    while (not stop_signal_called
           and (num_requested_samples > num_total_samps or num_requested_samples == 0)
           and (total_time == 0.0 or std::chrono::steady_clock::now() <= stop_time)) {

        int len = zmq_recv(subscriber, buffer, 100000, 0);

        const auto now = std::chrono::steady_clock::now();

        int32_t offset = 0;
        int32_t rv = vrt_read_header(buffer + offset, 100000 - offset, &h, true);

        /* Parse header */
        if (rv < 0) {
            fprintf(stderr, "Failed to parse header: %s\n", vrt_string_error(rv));
            return EXIT_FAILURE;
        }
        offset += rv;

        if (not start_rx and (h.packet_type == VRT_PT_IF_CONTEXT)) {
            start_rx = true;
            printf("Packet type: %s\n", vrt_string_packet_type(h.packet_type));

            /* Parse fields */
            rv = vrt_read_fields(&h, buffer + offset, 100000 - offset, &f, true);
            if (rv < 0) {
                fprintf(stderr, "Failed to parse fields section: %s\n", vrt_string_error(rv));
                return EXIT_FAILURE;
            }
            offset += rv;

            stream_id = f.stream_id;

            struct vrt_if_context c;
            rv = vrt_read_if_context(buffer + offset, 100000 - offset, &c, true);
            if (rv < 0) {
                fprintf(stderr, "Failed to parse IF context section: %s\n", vrt_string_error(rv));
                return EXIT_FAILURE;
            }
            if (c.has.sample_rate) {
                printf("Sample Rate [samples per second]: %.0f\n", c.sample_rate);
                sample_rate = (uint32_t)c.sample_rate;
            } else {
                printf("No Rate\n");
            }
            if (c.has.rf_reference_frequency) {
                printf("RF Freq [Hz]: %.0f\n", c.rf_reference_frequency);
                rf_freq = (int64_t)round(c.rf_reference_frequency);
            } else {
                printf("No Freq\n");
            }
            if (c.has.bandwidth) {
                printf("Bandwidth [Hz]: %.0f\n", c.bandwidth);
                bandwidth = c.bandwidth;
            } else {
                printf("No Bandwidth\n");
            }
            if (c.has.gain) {
                printf("Gain [dB]: %.0f\n", c.gain.stage1);
                gain = c.gain.stage1;
            } else {
                printf("No Gain\n");
            }
            if (c.state_and_event_indicators.has.reference_lock) {
                printf("Ref lock: %i\n", c.state_and_event_indicators.reference_lock);
                reflock = c.state_and_event_indicators.reference_lock;
            } else {
                printf("No Ref lock.\n");
            }
            if (c.state_and_event_indicators.has.calibrated_time) {
                printf("Time cal: %i\n", c.state_and_event_indicators.calibrated_time);
                time_cal = c.state_and_event_indicators.calibrated_time;
            } else {
                printf("No Time cal.\n");
            }
            start_time = now;
            stop_time = start_time + std::chrono::milliseconds(int64_t(1000 * total_time));
        }

        if (start_rx and (h.packet_type == VRT_PT_IF_DATA_WITH_STREAM_ID)) {

            if (not first_frame and (h.packet_count != (packet_count+1)%16) ) {
                printf("Error: lost frame (expected %i, received %i)\n", packet_count, h.packet_count);
                if (not continue_on_bad_packet)
                    break;
                else
                    packet_count = h.packet_count;
            } else {
                packet_count = h.packet_count;
            }

            /* Parse fields */
            rv = vrt_read_fields(&h, buffer + offset, 100000 - offset, &f, true);
            if (rv < 0) {
                fprintf(stderr, "Failed to parse fields section: %s\n", vrt_string_error(rv));
                return EXIT_FAILURE;
            }
            offset += rv;

            if (int_second) {
                // check if fractional second has wrapped
                if (f.fractional_seconds_timestamp > last_fractional_seconds_timestamp ) {
                        last_fractional_seconds_timestamp = f.fractional_seconds_timestamp;
                        continue;
                } else {
                    int_second = false;
                    last_update = now; 
                    start_time = now;
                    stop_time = start_time + std::chrono::milliseconds(int64_t(1000 * total_time));
                }
            }

            uint32_t num_rx_samps = (h.packet_size-offset);

            if (first_frame) {
                std::cout << boost::format(
                                 "First frame: %u samples, %u full secs, %.09f frac secs")
                                 % num_rx_samps
                                 % f.integer_seconds_timestamp
                                 % ((double)f.fractional_seconds_timestamp/1e12)
                          << std::endl;
                starttime_integer = f.integer_seconds_timestamp;
                starttime_fractional = f.fractional_seconds_timestamp;
                first_frame = false;
            }

            for (size_t i = 0; i < outfiles.size(); i++) {
                 outfiles[i]->write(
                    (const char*)&buffer[offset], sizeof(uint32_t)*h.packet_size-sizeof(uint32_t)*offset);
               }

            num_total_samps += num_rx_samps;
        }

        if (progress and not int_second) {
            last_update_samps += (h.packet_size-offset);
            const auto time_since_last_update = now - last_update;
            if (time_since_last_update > std::chrono::seconds(1)) {
                const double time_since_last_update_s =
                    std::chrono::duration<double>(time_since_last_update).count();
                const double rate = double(last_update_samps) / time_since_last_update_s;
                std::cout << "\t" << (rate / 1e6) << " Msps, ";
                
                last_update_samps = 0;
                last_update       = now;

                uint32_t num_rx_samps = (h.packet_size-offset);
    
                float sum_i = 0;
                uint32_t clip_i = 0;

                double datatype_max = 32768.;
                // if (cpu_format == "sc8" || cpu_format == "s8")
                //     datatype_max = 128.;

                for (int i=0; i<num_rx_samps; i++ ) {
                    auto sample_i = get_abs_val((std::complex<int16_t>)buffer[offset+i]);
                    sum_i += sample_i;
                    if (sample_i > datatype_max*0.99)
                        clip_i++;
                }
                sum_i = sum_i/num_rx_samps;
                std::cout << boost::format("%.0f") % (100.0*log2(sum_i)/log2(datatype_max)) << "% I (";
                std::cout << boost::format("%.0f") % ceil(log2(sum_i)+1) << " of ";
                std::cout << (int)ceil(log2(datatype_max)+1) << " bits), ";
                std::cout << "" << boost::format("%.0f") % (100.0*clip_i/num_rx_samps) << "% I clip, ";
                std::cout << std::endl;

            }
        }
    }

    for (size_t i = 0; i < outfiles.size(); i++) {
        outfiles[i]->close();
    }

    zmq_close(subscriber);
    zmq_ctx_destroy(context);

    std::cout << "Writing SigMF metadata..." << std::endl;

    ptree pt;
    ptree global;
    ptree captures;
    ptree global_item;
    ptree captures_item;

    // global
    global_item.put("core:version", "1.0.0"); 
    global_item.put("core:recorder", "difi_to_sigmf"); 
    global_item.put("core:sample_rate", sample_rate); 
    global_item.put("core:datatype", "ci16_le"); 
    // global_item.put("core:dataset", file);
    global_item.put("camras:usrp:rx_gain", gain); 
    // global_item.put("camras:usrp:channel", 0); 
    global_item.put("camras:usrp:bandwidth", bandwidth); 
    global_item.put("camras:usrp:reference", (reflock ? "external" : "internal")); 
    global_item.put("camras:usrp:time_source", (time_cal ? "pps" : "internal")); 
    global_item.put("camras:usrp:rx_gain", gain); 
    global_item.put("difi:stream_id", stream_id); 
    global.push_back(std::make_pair("", global_item));

    // captures
    captures_item.put("core:sample_start", 0); 
    captures_item.put("core:frequency", rf_freq); 
    captures_item.put("core:datetime", boost::format("%d.%06.0f") % starttime_integer % (double)(starttime_fractional/1e6)); 
    captures_item.put("core:datetime", boost::format("%s.%06.0f")
        % (boost::posix_time::to_iso_extended_string(boost::posix_time::from_time_t(starttime_integer)))
        % (double)(starttime_fractional/1e6)
    );
    captures.push_back(std::make_pair("", captures_item));
   
    pt.add_child("global", global);
    pt.add_child("captures", captures);

    write_json(mdfilename, pt);

    return 0;

}  
