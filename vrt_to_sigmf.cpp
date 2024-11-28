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
#include "dt-extended-context.h"

namespace po = boost::program_options;

static bool stop_signal_called = false;
void sig_int_handler(int)
{
    stop_signal_called = true;
}

//! Change to filename, e.g. from usrp_samples.dat to usrp_samples.chan0.dat,
//  but only if multiple names are to be generated.
std::string generate_out_filename(
    const std::string& base_fn, size_t n_names, size_t this_name, bool vrt=false)
{
    if (n_names == 1 or vrt) {
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
    std::string file, auto_file, type, zmq_address, channel_list, author, description;
    size_t num_requested_samples, total_time;
    uint16_t instance, main_port, port;
    int hwm;

    bool stream_has_dt_extended_context = false;

    // DT trace data
    float azimuth = NAN, elevation = NAN;

    // setup the program options
    po::options_description desc("Allowed options");
    // clang-format off

    desc.add_options()
        ("help", "help message")
        ("file", po::value<std::string>(&file)->default_value("vrt_samples"), "name of the file to write binary samples to")
        ("auto-file", po::value<std::string>(&auto_file), "prefix of the auto generated filename to write binary samples to")
        // ("type", po::value<std::string>(&type)->default_value("short"), "sample type: double, float, or short")
        ("nsamps", po::value<size_t>(&num_requested_samples)->default_value(0), "total number of samples to receive")
        ("duration", po::value<size_t>(&total_time)->default_value(0), "total number of seconds to receive")
        ("progress", "periodically display short-term bandwidth")
        ("channel", po::value<std::string>(&channel_list)->default_value("0"), "which channel(s) to use (specify \"0\", \"1\", \"0,1\", etc)")
        ("author", po::value<std::string>(&author), "core:author in sigmf-meta")
        ("description", po::value<std::string>(&description), "core:description in sigmf-meta")
        ("int-second", "align start of reception to integer second")
        ("null", "run without writing to file")
        ("continue", "don't abort on a bad packet")
        ("meta-only", "only create sigmf-meta file")
        ("dt-trace", "add DT trace data")
        ("vrt", "write VRT stream to file")
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
        std::cout << boost::format("VRT samples to file. %s") % desc << std::endl;
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
    bool meta_only              = vm.count("meta-only") > 0;
    bool dt_trace               = vm.count("dt-trace") > 0;
    bool do_auto_file           = vm.count("auto-file") > 0;
    bool int_second             = vm.count("int-second") > 0;
    bool has_author             = vm.count("author") > 0;
    bool has_desc               = vm.count("description") > 0;
    bool vrt                    = vm.count("vrt") > 0;
    bool zmq_split              = vm.count("zmq-split") > 0;

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

    std::vector<std::shared_ptr<std::ofstream>> outfiles;
    std::vector<std::shared_ptr<std::ofstream>> metafiles;

    std::string mdfilename;
    std::ofstream mdfile;
    mdfilename = file + ".sigmf-meta";
    if (vrt)
        file = file + ".sigmf-vrt";
    else
        file = file + ".sigmf-data";

    if (not null)
        for (size_t i = 0; i < channel_nums.size(); i++) {
            if (vrt and i > 0) break; // Only one data and metadata file for VRT
            const std::string meta_filename = generate_out_filename(mdfilename, channel_nums.size(), channel_nums[i], vrt);
            metafiles.push_back(std::shared_ptr<std::ofstream>(
                new std::ofstream(meta_filename.c_str())));

            if (not meta_only) {
                const std::string this_filename = generate_out_filename(file, channel_nums.size(), channel_nums[i], vrt);
                outfiles.push_back(std::shared_ptr<std::ofstream>(
                    new std::ofstream(this_filename.c_str(), std::ofstream::binary)));
            }
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

    // ZMQ buffer
    uint32_t buffer[ZMQ_BUFFER_SIZE];

    unsigned long long num_total_samps = 0;

    // Track time and samps between updating the BW summary
    auto last_update                     = start_time;
    unsigned long long last_update_samps = 0;

    uint64_t last_fractional_seconds_timestamp = 0;

    bool first_frame = true;
    bool dt_trace_received = false;
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

        std::string channel = std::to_string(channel_nums[ch]);
        if (vrt) {
            channel = channel_list;
        }

        if (vrt and not null and not meta_only)
            outfiles[0]->write((const char*)&buffer, len);

        if ( not (context_recv & vrt_packet.stream_id) and vrt_packet.context
             and not first_frame and not (dt_trace and not dt_ext_context.dt_ext_context_received)) {

            context_recv |= vrt_packet.stream_id;

            vrt_print_context(&vrt_context);

            if (not null) {
                // std::cout << "Writing SigMF metadata..." << std::endl;
                if (metafiles.size() < ch + 1) {
                    if (!vrt)
                        std::cout << "File not created?!";
                } else {
                    std::string json = str(boost::format("{ \n"
                    "    \"global\": {\n"
                    "        \"core:version\": \"1.0.0\",\n"
                    "        \"core:recorder\": \"vrt_to_sigmf\",\n"
                    "        \"core:sample_rate\": %u,\n") % vrt_context.sample_rate);
                    if (vrt) {
                        json += str(boost::format(
                        "        \"core:datatype\": \"vrt\",\n"));
                    } else {
                        json += str(boost::format(
                        "        \"core:datatype\": \"ci16_le\",\n"));
                    }
                    if (vrt and not do_auto_file)
                        json += str(boost::format(
                        "        \"core:dataset\": \"%s\",\n")
                        % file);
                    if (has_author) {
                        json += str(boost::format(
                        "        \"core:author\": \"%s\",\n")
                        % author);
                    }
                    if (has_desc) {
                        json += str(boost::format(
                        "        \"core:description\": \"%s\",\n")
                        % description);
                    }
                    if (dt_trace) {
                        char const *trackerStrings[] = {"idle", "azel", "j2000tracker", "moontracker", "suntracker", "sattracker", "manual"};
                        json += str(boost::format(
                        "        \"dt:datetime\": \"%s.%06.0f\",\n"
                        "        \"dt:pointing:active_tracker\": \"%s\",\n"
                        "        \"dt:pointing:tracking_enabled\": \"%s\",\n"
                        "        \"dt:pointing:refraction\": \"%s\",\n"
                        "        \"dt:pointing:dt_model\": \"%s\",\n"
                        "        \"dt:pointing:refraction_j2000\": \"%s\",\n"
                        "        \"dt:pointing:dt_model_j2000\": \"%s\",\n"
                        "        \"dt:pointing:current:az_deg\": %.3f,\n"
                        "        \"dt:pointing:current:el_deg\": %.3f,\n"
                        "        \"dt:pointing:error:az_deg\": %.3f,\n"
                        "        \"dt:pointing:error:el_deg\": %.3f,\n"
                        "        \"dt:pointing:offset:az_deg\": %.3f,\n"
                        "        \"dt:pointing:offset:el_deg\": %.3f,\n"
                        "        \"dt:pointing:az_speed_deg_s\": %.6f,\n"
                        "        \"dt:pointing:el_speed_deg_s\": %.6f,\n"
                        "        \"dt:pointing:setpoint:ra_h\": %.3f,\n"
                        "        \"dt:pointing:setpoint:dec_deg\": %.3f,\n"
                        "        \"dt:pointing:current:ra_h\": %.3f,\n"
                        "        \"dt:pointing:current:dec_deg\": %.3f,\n"
                        "        \"dt:pointing:model:a0\": %.6f,\n"
                        "        \"dt:pointing:model:c1\": %.6f,\n"
                        "        \"dt:pointing:model:c2\": %.6f,\n"
                        "        \"dt:pointing:model:e0\": %.6f,\n"
                        "        \"dt:pointing:model:b\": %.6f,\n"
                        "        \"dt:pointing:model:za\": %.6f,\n"
                        "        \"dt:pointing:model:aa\": %.6f,\n"
                        "        \"dt:focusbox_position_mm\": %.0f,\n" )
                        % (boost::posix_time::to_iso_extended_string(boost::posix_time::from_time_t(dt_ext_context.integer_seconds_timestamp)))
                        % ((double)(dt_ext_context.fractional_seconds_timestamp/1e6))
                        % (trackerStrings[dt_ext_context.active_tracker])
                        % (dt_ext_context.tracking_enabled ? "true" : "false")
                        % (dt_ext_context.refraction ? "true" : "false")
                        % (dt_ext_context.dt_model ? "true" : "false")
                        % (dt_ext_context.refraction_j2000 ? "true" : "false")
                        % (dt_ext_context.dt_model_j2000 ? "true" : "false")
                        % ((180.0/M_PI)*dt_ext_context.azimuth)
                        % ((180.0/M_PI)*dt_ext_context.elevation)
                        % ((180.0/M_PI)*dt_ext_context.azimuth_error)
                        % ((180.0/M_PI)*dt_ext_context.elevation_error)
                        % ((180.0/M_PI)*dt_ext_context.azimuth_offset)
                        % ((180.0/M_PI)*dt_ext_context.elevation_offset)
                        % ((180.0/M_PI)*dt_ext_context.azimuth_speed)
                        % ((180.0/M_PI)*dt_ext_context.elevation_speed)
                        % ((12.0/M_PI)*dt_ext_context.ra_setpoint)
                        % ((180.0/M_PI)*dt_ext_context.dec_setpoint)
                        % ((12.0/M_PI)*dt_ext_context.ra_current)
                        % ((180.0/M_PI)*dt_ext_context.dec_current)
                        % dt_ext_context.model_a0
                        % dt_ext_context.model_c1
                        % dt_ext_context.model_c2
                        % dt_ext_context.model_e0
                        % dt_ext_context.model_b
                        % dt_ext_context.model_za
                        % dt_ext_context.model_aa
                        % dt_ext_context.focusbox );
                    }
                    if (vrt_context.timestamp_calibration_time != 0) {
                        json += str(boost::format("        \"vrt:cal_time\": %u,\n") % vrt_context.timestamp_calibration_time);
                    }
                    json += str(boost::format(
                    "        \"vrt:rx_gain\": %i,\n"
                    "        \"vrt:bandwidth\": %u,\n"
                    "        \"vrt:reference\": \"%s\",\n"
                    "        \"vrt:time_source\": \"%s\",\n"
                    "        \"vrt:stream_id\": %u,\n"
                    "        \"vrt:channel\": %s\n"
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
                    % vrt_context.gain
                    % vrt_context.bandwidth
                    % (vrt_context.reflock ? "external" : "internal")
                    % (vrt_context.time_cal ? "pps" : "internal")
                    % vrt_context.stream_id
                    % channel
                    % vrt_context.rf_freq
                    % (boost::posix_time::to_iso_extended_string(boost::posix_time::from_time_t(vrt_context.starttime_integer)))
                    % (double)(vrt_context.starttime_fractional/1e6) );
                    *metafiles[ch] << json;
                    *metafiles[ch] << std::endl;
                    metafiles[ch]->close();
                    if (meta_only and ch==(channel_nums.size()-1))
                        break;
                }
            }

            if (total_time > 0)
                num_requested_samples = total_time * vrt_context.sample_rate;
        }


        if (vrt_packet.extended_context) {
            if (stream_has_dt_extended_context and not dt_trace) {
                std::cerr << "WARNING: DT metadata is present in the stream, but it is ignored. Did you forget --dt-trace?" << std::endl;
            }
            stream_has_dt_extended_context = dt_process(buffer, sizeof(buffer), &vrt_packet, &dt_ext_context);
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
            if (not vrt and not null and not meta_only) {
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

    // Auto file
    if (do_auto_file) {
        boost::format auto_format;
        boost::posix_time::ptime starttime = boost::posix_time::from_time_t(vrt_context.starttime_integer);
        std::string timestring = boost::posix_time::to_iso_extended_string(starttime);
        std::replace( timestring.begin(), timestring.end(), ':', '_');
        std::replace( timestring.begin(), timestring.end(), '-', '_');
        std::replace( timestring.begin(), timestring.end(), 'T', '_');
        auto_format = boost::format("%s_%s_%.3fMHz_%.2fMsps_ci16_le")
                    % (auto_file)
                    % (timestring)
                    % (vrt_context.rf_freq/1e6)
                    % (vrt_context.sample_rate/1e6);

        std::string auto_mdfilename = auto_format.str() + ".sigmf-meta";
        std::string auto_bin_file;
        if (vrt)
            auto_bin_file = auto_format.str() + ".sigmf-vrt";
        else
            auto_bin_file = auto_format.str() + ".sigmf-data";

        if (not null)
            for (size_t i = 0; i < channel_nums.size(); i++) {
                if (vrt and i > 0) break;
                const std::string meta_filename = generate_out_filename(mdfilename, channel_nums.size(), channel_nums[i], vrt);
                const std::string auto_meta_filename = generate_out_filename(auto_mdfilename, channel_nums.size(), channel_nums[i], vrt);
                boost::filesystem::rename(meta_filename, auto_meta_filename);

                if (not meta_only) {
                    if (vrt and i > 0) break; // Only one data file for VRT
                    const std::string this_filename = generate_out_filename(file, channel_nums.size(), channel_nums[i], vrt);
                    const std::string this_auto_filename = generate_out_filename(auto_bin_file, channel_nums.size(), channel_nums[i], vrt);
                    boost::filesystem::rename( this_filename, this_auto_filename );
                }

        }
    }

    zmq_close(subscriber);
    zmq_ctx_destroy(context);

    return 0;

}
