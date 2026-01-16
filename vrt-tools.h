/* VRT tools helper functions */

#ifndef _VRTTOOLS_H
#define _VRTTOOLS_H

#define VRT_SAMPLES_PER_PACKET 10000

#define SIZE (VRT_SAMPLES_PER_PACKET+7)
#define VRT_DATA_PACKET_SIZE (VRT_SAMPLES_PER_PACKET+7)

#define ZMQ_BUFFER_SIZE 100000

#define MAX_CHANNELS    10

#define DEFAULT_MAIN_PORT       50100
#define DEFAULT_GNURADIO_PORT   (DEFAULT_MAIN_PORT+100)
#define DEFAULT_CONTROL_PORT    (DEFAULT_MAIN_PORT+200)
#define DEFAULT_TX_PORT         (DEFAULT_MAIN_PORT+400)

// Context update interval in ms
#define VRT_CONTEXT_INTERVAL 200

// VRT
#include <vrt/vrt_init.h>
#include <vrt/vrt_string.h>
#include <vrt/vrt_types.h>
#include <vrt/vrt_util.h>
#include <vrt/vrt_write.h>
#include <vrt/vrt_read.h>

struct context_type {
    bool context_received;
    bool context_changed;
    int64_t rf_freq;
    double rf_frac_freq;
    uint32_t sample_rate;
    int32_t gain;
    float temperature;
    uint32_t bandwidth;
    bool reflock;
    bool time_cal;
    uint32_t stream_id;
    uint64_t starttime_integer;
    uint64_t starttime_fractional;
    int32_t last_data_counter;
    uint64_t fractional_seconds_timestamp;
    uint64_t integer_seconds_timestamp;
    uint32_t timestamp_calibration_time;
};

struct packet_type {
    bool context;
    bool data;
    bool extended_context;
    bool lost_frame;
    bool first_frame;
    uint32_t oui;
    uint16_t information_class_code;
    uint16_t packet_class_code;
    uint32_t stream_id;
    uint32_t channel_filt;
    uint32_t num_rx_samps;
    uint32_t offset;
    uint64_t fractional_seconds_timestamp;
    uint64_t integer_seconds_timestamp;
};

void init_context(context_type* context) {
    context->context_received = false;
    context->context_changed = false;
    context->last_data_counter = -1;
    context->rf_freq = 0;
    context->sample_rate = 0;
    context->gain = 0;
    context->bandwidth = 0;
    context->stream_id = 0;
    context->starttime_integer = 0;
    context->starttime_fractional = 0;
    context->reflock = false;
    context->time_cal = false;
    context->timestamp_calibration_time = 0;
}

bool check_packet_count(int8_t counter, context_type* vrt_context) {
    if ( (vrt_context->last_data_counter > 0) and
            ( (counter != (vrt_context->last_data_counter+1)%16) and
              (counter != (vrt_context->last_data_counter  )%16) ) ) {
        printf("# Error: lost frame (expected %i, received %i)\n", vrt_context->last_data_counter, counter);
        vrt_context->last_data_counter = counter;
        return false;
    } else {
        vrt_context->last_data_counter = counter;
        return true;
    }
}

void vrt_print_context(context_type* vrt_context) {

    uint32_t ch=0;
    while(not (vrt_context->stream_id & (1 << ch) ) )
            ch++;

    printf("# VRT Context:\n");
    printf("#    Stream ID (channel): %u (%u)\n", vrt_context->stream_id, ch);
    printf("#    Sample Rate [samples per second]: %i\n", vrt_context->sample_rate);
    printf("#    RF Freq [Hz]: %lld\n", (long long int)vrt_context->rf_freq);
    printf("#    RF frac. Freq [Hz]: %e\n", vrt_context->rf_frac_freq);
    printf("#    Bandwidth [Hz]: %i\n", vrt_context->bandwidth);
    printf("#    Gain [dB]: %i\n", vrt_context->gain);
    printf("#    Ref lock: %s\n", vrt_context->reflock == 1 ? "external" : "internal");
    printf("#    Time cal: %s\n", vrt_context->time_cal == 1? "pps" : "internal");
    if (vrt_context->timestamp_calibration_time != 0)
        printf("#    Cal time: %u\n", vrt_context->timestamp_calibration_time);

}

bool vrt_process(uint32_t* buffer, uint32_t size, context_type* vrt_context, packet_type* vrt_packet) {

    struct vrt_header h;
    struct vrt_fields f;

    int32_t offset = 0;
    int32_t rv = vrt_read_header(buffer + offset, size - offset, &h, true);

    vrt_packet->context = false;
    vrt_packet->data = false;
    vrt_packet->extended_context = false;

    /* Parse header */
    if (rv < 0) {
        fprintf(stderr, "Failed to parse header: %s\n", vrt_string_error(rv));
        return false;
    }
    offset += rv;

    if (h.packet_type == VRT_PT_IF_CONTEXT) {
        // Context

        /* Parse fields */
        rv = vrt_read_fields(&h, buffer + offset, size - offset, &f, true);
        if (rv < 0) {
            fprintf(stderr, "Failed to parse fields section: %s\n", vrt_string_error(rv));
            return false;
        }
        offset += rv;

        vrt_context->stream_id = f.stream_id;

        if (f.stream_id & vrt_packet->channel_filt) {
            struct vrt_if_context c;
            rv = vrt_read_if_context(buffer + offset, ZMQ_BUFFER_SIZE - offset, &c, true);
            if (rv < 0) {
                fprintf(stderr, "Failed to parse IF context section: %s\n", vrt_string_error(rv));
                return false;
            }

            vrt_context->integer_seconds_timestamp = f.integer_seconds_timestamp;
            vrt_context->fractional_seconds_timestamp = f.fractional_seconds_timestamp;
            if (c.has.sample_rate)
                vrt_context->sample_rate = (uint32_t)round(c.sample_rate);

            if (c.has.rf_reference_frequency) {
                vrt_context->rf_freq = (int64_t)round(c.rf_reference_frequency);
                vrt_context->rf_frac_freq = c.rf_reference_frequency - (double)vrt_context->rf_freq;
            }

            if (c.has.bandwidth)
                vrt_context->bandwidth = c.bandwidth;

            if (c.has.gain)
                vrt_context->gain = c.gain.stage1;

            if (c.state_and_event_indicators.has.reference_lock)
                vrt_context->reflock = c.state_and_event_indicators.reference_lock;

            if (c.state_and_event_indicators.has.calibrated_time)
                vrt_context->time_cal = c.state_and_event_indicators.calibrated_time;

            if (c.has.temperature)
                vrt_context->temperature = c.temperature;

            if (c.has.timestamp_calibration_time)
                vrt_context->timestamp_calibration_time = c.timestamp_calibration_time;

            vrt_context->context_changed = c.context_field_change_indicator;
            vrt_packet->context = true;
            vrt_context->context_received = true;
            vrt_packet->stream_id = f.stream_id;

            vrt_packet->oui = f.class_id.oui;
            vrt_packet->information_class_code = f.class_id.information_class_code;
            vrt_packet->packet_class_code = f.class_id.packet_class_code;

        }
    } else if (h.packet_type == VRT_PT_IF_DATA_WITH_STREAM_ID) {
        // Data
        /* Parse fields */
        rv = vrt_read_fields(&h, buffer + offset, ZMQ_BUFFER_SIZE - offset, &f, true);
        if (rv < 0) {
            fprintf(stderr, "Failed to parse fields section: %s\n", vrt_string_error(rv));
            return false;
        }
        offset += rv;
        if (f.stream_id & vrt_packet->channel_filt) {

            if (not check_packet_count(h.packet_count, vrt_context))
                vrt_packet->lost_frame = true;
            else
                vrt_packet->lost_frame = false;

            vrt_packet->integer_seconds_timestamp = f.integer_seconds_timestamp;
            vrt_packet->fractional_seconds_timestamp = f.fractional_seconds_timestamp;
            vrt_packet->num_rx_samps = (h.packet_size-offset);
            vrt_packet->offset = offset;
            vrt_packet->stream_id = f.stream_id;
            vrt_packet->data = true;

            vrt_packet->oui = f.class_id.oui;
            vrt_packet->information_class_code = f.class_id.information_class_code;
            vrt_packet->packet_class_code = f.class_id.packet_class_code;

            if (vrt_packet->first_frame) {
                vrt_context->starttime_integer = f.integer_seconds_timestamp;
                vrt_context->starttime_fractional = f.fractional_seconds_timestamp;
                vrt_packet->first_frame = false;
            }
        }
    } else if (h.packet_type == VRT_PT_EXT_CONTEXT) {

         /* Parse fields */
        rv = vrt_read_fields(&h, buffer + offset, size - offset, &f, true);
        if (rv < 0) {
            fprintf(stderr, "Failed to parse fields section: %s\n", vrt_string_error(rv));
            return false;
        }
        offset += rv;

        vrt_packet->integer_seconds_timestamp = f.integer_seconds_timestamp;
        vrt_packet->fractional_seconds_timestamp = f.fractional_seconds_timestamp;
        vrt_packet->num_rx_samps = (h.packet_size-offset);
        vrt_packet->offset = offset;
        vrt_packet->stream_id = f.stream_id;

        vrt_packet->oui = f.class_id.oui;
        vrt_packet->information_class_code = f.class_id.information_class_code;
        vrt_packet->packet_class_code = f.class_id.packet_class_code;

        vrt_packet->extended_context = true;
    }

    return true;
}

void vrt_init_data_packet(struct vrt_packet* p) {

    p->header.packet_type         = VRT_PT_IF_DATA_WITH_STREAM_ID;

    p->header.packet_size         = SIZE;
    p->header.tsm                 = VRT_TSM_FINE;
    p->header.tsi                 = VRT_TSI_OTHER; // unix time
    p->header.tsf                 = VRT_TSF_REAL_TIME;
    p->fields.stream_id           = 0;
    p->words_body                 = VRT_SAMPLES_PER_PACKET;

    p->header.has.class_id        = true;
    p->fields.class_id.oui        = 0xFF5454;
    p->fields.class_id.information_class_code = 0;
    p->fields.class_id.packet_class_code = 0;

    p->header.has.trailer         = false;
}

void vrt_init_context_packet(struct vrt_packet* pc) {

    pc->header.packet_type = VRT_PT_IF_CONTEXT;
    pc->header.has.class_id = true;

    pc->fields.class_id.oui        = 0xFF5454;
    pc->fields.class_id.information_class_code = 0;
    pc->fields.class_id.packet_class_code = 0;

    pc->if_context.has.bandwidth   = true;
    pc->if_context.has.sample_rate = true;
    pc->if_context.has.reference_point_identifier = true;
    pc->if_context.has.if_reference_frequency = true;
    pc->if_context.has.rf_reference_frequency = true;
    pc->if_context.has.if_band_offset = true;
    pc->if_context.has.reference_level = true;
    pc->if_context.has.gain = true;
    pc->if_context.has.timestamp_adjustment = true;
    pc->if_context.has.timestamp_calibration_time = true;
    pc->if_context.has.state_and_event_indicators = true;
    pc->if_context.has.data_packet_payload_format = true;

    pc->if_context.data_packet_payload_format.packing_method = VRT_PM_LINK_EFFICIENT;
    pc->if_context.data_packet_payload_format.real_or_complex = VRT_ROC_COMPLEX_CARTESIAN;
    pc->if_context.data_packet_payload_format.data_item_format = VRT_DIF_SIGNED_FIXED_POINT;
    pc->if_context.data_packet_payload_format.sample_component_repeat = false;
    pc->if_context.data_packet_payload_format.item_packing_field_size = 31;
    pc->if_context.data_packet_payload_format.data_item_size = 15;

    pc->header.tsm                 = VRT_TSM_COARSE;
    pc->header.tsi                 = VRT_TSI_OTHER; // unix time
    pc->header.tsf                 = VRT_TSF_REAL_TIME;

    pc->if_context.state_and_event_indicators.has.reference_lock = true;
    pc->if_context.state_and_event_indicators.has.calibrated_time = true;

}

void show_progress_stats(
    std::chrono::time_point<std::chrono::steady_clock> now,
    std::chrono::time_point<std::chrono::steady_clock> *last_update,
    uint64_t *last_update_samps,
    uint32_t *buffer,
    size_t num_rx_samps,
    uint32_t channel) {

    *last_update_samps += num_rx_samps;

    const auto time_since_last_update = now - *last_update;
    if (time_since_last_update > std::chrono::seconds(1)) {
        const double time_since_last_update_s =
            std::chrono::duration<double>(time_since_last_update).count();
        const double rate = double(*last_update_samps) / time_since_last_update_s;
        *last_update_samps = 0;
        *last_update       = now;

        double max_iq = 0;
        uint32_t clip_iq = 0;

        double datatype_max = 32767.;

        for (int i=0; i < num_rx_samps; i++ ) {
            std::complex<int16_t> sample = (std::complex<int16_t>)buffer[i];
            max_iq = fmax(max_iq, fmax(fabs(sample.real()), fabs(sample.imag())));
            if (fabs(sample.real()) > datatype_max*0.99 || fabs(sample.imag()) > datatype_max*0.99)
                clip_iq++;
        }
        std::cout << "\t" << boost::format("%.6f") % (rate / 1e6) << " Msps, ";
        std::cout << "CH" << boost::format("%u") % channel << ": ";
        std::cout << boost::format("%3.0f") % (20*log10(max_iq/datatype_max)) << " dBFS (";
        std::cout << boost::format("%2.0f") % ceil(log2(max_iq)+1) << "/";
        std::cout << (int)ceil(log2(datatype_max)+1) << " bits), ";
        std::cout << "" << boost::format("%2.0f") % (100.0*clip_iq/num_rx_samps) << "% clip. ";
        std::cout << std::endl;
    }

}

#endif
