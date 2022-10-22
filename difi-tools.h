/* DIFI tools helper functions */

#ifndef _DIFITOOLS_H
#define _DIFITOOLS_H

#define DIFI_SAMPLES_PER_PACKET 10000

#define SIZE (DIFI_SAMPLES_PER_PACKET+7)
#define DIFI_DATA_PACKET_SIZE (DIFI_SAMPLES_PER_PACKET+7)

#define ZMQ_BUFFER_SIZE 100000

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
    uint32_t sample_rate;
    int32_t gain;
    uint32_t bandwidth;
    bool reflock;
    bool time_cal;
    uint32_t stream_id;
    uint64_t starttime_integer;
    uint64_t starttime_fractional;
    int8_t last_data_counter;
};

struct difi_packet_type {
    bool context;
    bool data;
    uint32_t data_offset;
    bool lost_frame;
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
}

bool check_packet_count(int8_t counter, context_type* difi_context) {
    if ( (difi_context->last_data_counter > 0) and (counter != (difi_context->last_data_counter+1)%16) ) {
        printf("Error: lost frame (expected %i, received %i)\n", difi_context->last_data_counter, counter);
        difi_context->last_data_counter = counter;
        return false;
    } else {
        difi_context->last_data_counter = counter;
        return true;
    }
}

void difi_print_context(context_type* difi_context) {

    printf("  DIFI Context:\n");
    printf("    Sample Rate [samples per second]: %i\n", difi_context->sample_rate);
    printf("    RF Freq [Hz]: %li\n", difi_context->rf_freq);
    printf("    Bandwidth [Hz]: %i\n", difi_context->bandwidth);
    printf("    Gain [dB]: %i\n", difi_context->gain);
    printf("    Ref lock: %i\n", difi_context->reflock);
    printf("    Time cal: %i\n", difi_context->time_cal);

}

bool difi_process(uint32_t* buffer, uint32_t size, context_type* difi_context, difi_packet_type* difi_packet) {

    struct vrt_header h;
    struct vrt_fields f;

    int32_t offset = 0;
    int32_t rv = vrt_read_header(buffer + offset, size - offset, &h, true);

    difi_packet->context = false;
    difi_packet->data = false;

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

        difi_context->stream_id = f.stream_id;

        struct vrt_if_context c;
        rv = vrt_read_if_context(buffer + offset, 100000 - offset, &c, true);
        if (rv < 0) {
            fprintf(stderr, "Failed to parse IF context section: %s\n", vrt_string_error(rv));
            return false;
        }
        if (c.has.sample_rate) {
            difi_context->sample_rate = (uint32_t)c.sample_rate;
        } else {
            printf("No Rate\n");
        }
        if (c.has.rf_reference_frequency) {
            difi_context->rf_freq = (int64_t)round(c.rf_reference_frequency);
        } else {
            printf("No Freq\n");
        }
        if (c.has.bandwidth) {
            difi_context->bandwidth = c.bandwidth;
        } else {
            printf("No Bandwidth\n");
        }
        if (c.has.gain) {
            difi_context->gain = c.gain.stage1;
        } else {
            printf("No Gain\n");
        }
        if (c.state_and_event_indicators.has.reference_lock) {
            difi_context->reflock = c.state_and_event_indicators.reference_lock;
        } else {
            printf("No Ref lock.\n");
        }
        if (c.state_and_event_indicators.has.calibrated_time) {
            difi_context->time_cal = c.state_and_event_indicators.calibrated_time;
        } else {
            printf("No Time cal.\n");
        }
        difi_packet->context = true;
        difi_context->context_received = true;

    } else if (h.packet_type == VRT_PT_IF_DATA_WITH_STREAM_ID) {
        // Data
        if (not check_packet_count(h.packet_count, difi_context))
            difi_packet->lost_frame = true;
        else 
            difi_packet->lost_frame = false;
       /* Parse fields */
        rv = vrt_read_fields(&h, buffer + offset, 100000 - offset, &f, true);
        if (rv < 0) {
            fprintf(stderr, "Failed to parse fields section: %s\n", vrt_string_error(rv));
            return false;
        }
        offset += rv;

        difi_packet->integer_seconds_timestamp = f.integer_seconds_timestamp;
        difi_packet->fractional_seconds_timestamp = f.fractional_seconds_timestamp;
        difi_packet->num_rx_samps = (h.packet_size-offset);
        difi_packet->offset = offset;
        difi_packet->data = true;
    }

    return true;
}

void difi_init_data_packet(struct vrt_packet* p) {

	p->header.packet_type         = VRT_PT_IF_DATA_WITH_STREAM_ID;

    p->header.packet_size         = SIZE;
    p->header.tsm                 = VRT_TSM_FINE;
    p->header.tsi                 = VRT_TSI_OTHER; // unix time
    p->header.tsf                 = VRT_TSF_REAL_TIME;
    p->fields.stream_id           = 0;
    p->words_body                 = 10000;

    p->header.has.class_id        = true;
    p->fields.class_id.oui        = 0x6A621E; // DIFI OUI
    p->fields.class_id.information_class_code = 0;
    p->fields.class_id.packet_class_code = 0;

    p->header.has.trailer         = false;
}

void difi_init_context_packet(struct vrt_packet* pc) {

	pc->header.packet_type = VRT_PT_IF_CONTEXT;
    pc->header.has.class_id = true;
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

#endif