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
    // variables to be set by po
    std::string zmq_address;
    size_t num_requested_samples;
    float update_time;
    double total_time;
    uint16_t port;
    uint32_t channel;
    int hwm;

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
        // ("int-interval", "align start of reception to integer number of integration intervals (implies --int-second)")
        ("update-time", po::value<float>(&update_time)->default_value(1.0), "update time (seconds)")
        // ("ecsv", "output in ECSV format (Astropy)")
        ("temperature", "output temperature")
        ("null", "run without writing to file")
        ("continue", "don't abort on a bad packet")
        ("dt-trace", "use DT trace data in VRT stream")
        ("address", po::value<std::string>(&zmq_address)->default_value("localhost"), "VRT ZMQ address")
        ("port", po::value<uint16_t>(&port)->default_value(50100), "VRT ZMQ port")
        ("hwm", po::value<int>(&hwm)->default_value(10000), "VRT ZMQ HWM")
    ;
    // clang-format on
    po::variables_map vm;
    // po::store(po::parse_command_line(argc, argv, desc), vm);
    po::store(po::command_line_parser(argc, argv).options(desc).style(po::command_line_style::unix_style ^ po::command_line_style::allow_short).run(), vm);
    po::notify(vm);

    // print the help message
    if (vm.count("help")) {
        std::cout << boost::format("VRT samples to gnuplot %s") % desc << std::endl;
        std::cout << std::endl
                  << "This application streams data from a VRT stream "
                     "to gnuplot.\n"
                  << std::endl;
        return ~0;
    }

    bool progress               = vm.count("progress") > 0;
    bool stats                  = vm.count("stats") > 0;
    bool null                   = vm.count("null") > 0;
    bool continue_on_bad_packet = vm.count("continue") > 0;
    bool int_interval           = (bool)vm.count("int-interval");
    bool int_second             = int_interval || (bool)vm.count("int-second");
    bool dt_trace               = vm.count("dt-trace") > 0;
    bool log_temp               = vm.count("temperature") > 0;
   
    context_type vrt_context;
    dt_ext_context_type dt_ext_context;
    init_context(&vrt_context);

    packet_type vrt_packet;

    vrt_packet.channel_filt = 1<<channel;

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
    int32_t samples_per_update = 0;
    int32_t samples_last_update = 0;

    // Track time and samps between updating the BW summary
    auto last_update                     = start_time;
    unsigned long long last_update_samps = 0;

    bool start_rx = false;
    uint64_t last_fractional_seconds_timestamp = 0;

    while (not stop_signal_called
           and (num_requested_samples > num_total_samps or num_requested_samples == 0)
           and (total_time == 0.0 or std::chrono::steady_clock::now() <= stop_time)) {

        int len = zmq_recv(subscriber, buffer, ZMQ_BUFFER_SIZE, 0);

        const auto now = std::chrono::steady_clock::now();

        if (not vrt_process(buffer, sizeof(buffer), &vrt_context, &vrt_packet)) {
            printf("Not a Vita49 packet?\n");
            continue;
        }

        if (vrt_packet.context) {
            samples_per_update = ((double)vrt_context.sample_rate * update_time);
            if (total_time > 0)  
                num_requested_samples = total_time * vrt_context.sample_rate;
        }

        if (not start_rx and vrt_packet.context) {
            // vrt_print_context(&vrt_context);
            start_rx = true;

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
            printf("# - metadata: !!omap\n");
            printf("#   - {update_time: %.2f}\n", (double)update_time);
            
            printf("# datatype:\n");
            printf("# - {name: data_timestamp, datatype: float64}\n");
            printf("# - {name: context_timestamp, datatype: float64}\n");
            printf("# - {name: center_freq_hz, unit: Hz, datatype: float64}\n");
            printf("# - {name: sample_rate, datatype: float64}\n");
            printf("# - {name: rx_gain, unit: dB, datatype: float64}\n");
            if (log_temp) {
                printf("# - {name: temperature_deg_c, datatype: float64}\n");
            }
            if (dt_trace) {
                printf("# - {name: ext_context_timestamp, datatype: float64}\n");
                printf("# - {name: current_az_deg, unit: deg, datatype: float64}\n");
                printf("# - {name: current_el_deg, unit: deg, datatype: float64}\n");
                printf("# - {name: current_az_error_deg, unit: deg, datatype: float64}\n");
                printf("# - {name: current_el_error_deg, unit: deg, datatype: float64}\n");
                printf("# - {name: current_ra_h, unit: h, datatype: float64}\n");
                printf("# - {name: current_dec_deg, unit: deg, datatype: float64}\n");
                printf("# - {name: setpoint_ra_h, unit: deg, datatype: float64}\n");
                printf("# - {name: setpoint_dec_deg, unit: deg, datatype: float64}\n");
                printf("# - {name: radec_error_angle_deg, unit: deg, datatype: float64}\n");
                printf("# - {name: radec_error_bearing_deg, unit: deg, datatype: float64}\n");
                printf("# - {name: focusbox_mm, unit: mm, datatype: float64}\n");
            }
            printf("# schema: astropy-2.0\n");
        
            // Header
            printf("data_timestamp, context_timestamp, center_freq_hz, sample_rate, rx_gain");
            if (log_temp)
                printf(", temperature_deg_c");
            if (dt_trace) 
                printf(", ext_context_timestamp, current_az_deg, current_el_deg, current_az_error_deg, current_el_error_deg, current_ra_h, current_dec_deg, setpoint_ra_h, setpoint_dec_deg, radec_error_angle_deg, radec_error_bearing_deg, focusbox_mm");
            printf("\n");
            fflush(stdout);
        }

        if (vrt_packet.extended_context) {
            // TODO: Find some other variable to avoid giving this warning for every extended context packet
            // This now assumes that any extended context is a DT extended context
            if (not dt_ext_context.dt_ext_context_received and not dt_trace) {
                std::cerr << "# WARNING: DT metadata is present in the stream, but it is ignored. Did you forget --dt-trace?" << std::endl;
            }
            dt_process(buffer, sizeof(buffer), &vrt_packet, &dt_ext_context);
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
                    // if (int_interval && int(vrt_packet.integer_seconds_timestamp) % int(update_time) != 0) {
                    //     continue;
                    // }
                    int_second = false;
                    last_update = now; 
                    start_time = now;
                }
            }

            num_total_samps += vrt_packet.num_rx_samps;

            if ( num_total_samps-samples_last_update >= samples_per_update || samples_last_update==0) {

                uint64_t data_seconds = vrt_packet.integer_seconds_timestamp;
                uint64_t data_frac_seconds = vrt_packet.fractional_seconds_timestamp;

                printf("%lu.%09li", static_cast<unsigned long>(data_seconds), static_cast<long>(data_frac_seconds/1e3));
                printf(", %lu.%09li", static_cast<unsigned long>(vrt_context.integer_seconds_timestamp), static_cast<long>(vrt_context.fractional_seconds_timestamp/1e3));
                printf(", %li", static_cast<long>(vrt_context.rf_freq));
                printf(", %li", static_cast<long>(vrt_context.sample_rate));
                printf(", %li", static_cast<long>(vrt_context.gain));

                if (log_temp)
                    printf(", %.2f", vrt_context.temperature); 
                if (dt_trace) {
                    printf(", %lu.%09li", static_cast<unsigned long>(dt_ext_context.integer_seconds_timestamp), static_cast<long>(dt_ext_context.fractional_seconds_timestamp/1e3));
                    printf(", %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f",
                        ((180.0/M_PI)*dt_ext_context.azimuth),
                        ((180.0/M_PI)*dt_ext_context.elevation),
                        ((180.0/M_PI)*dt_ext_context.azimuth_error),
                        ((180.0/M_PI)*dt_ext_context.elevation_error),
                        ((12.0/M_PI)*dt_ext_context.ra_current),
                        ((180.0/M_PI)*dt_ext_context.dec_current),
                        ((12.0/M_PI)*dt_ext_context.ra_setpoint),
                        ((180.0/M_PI)*dt_ext_context.dec_setpoint),
                        ((180.0/M_PI)*haversine(dt_ext_context.dec_setpoint, dt_ext_context.dec_current, dt_ext_context.ra_setpoint, dt_ext_context.ra_current)),
                        ((180.0/M_PI)*bearing(dt_ext_context.dec_setpoint, dt_ext_context.dec_current, dt_ext_context.ra_setpoint, dt_ext_context.ra_current)),
                        dt_ext_context.focusbox);
                }

                printf("\n");
                fflush(stdout);

                samples_last_update = num_total_samps;
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
                std::cout << "# " << (rate / 1e6) << " Msps, ";
                
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
