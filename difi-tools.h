/* DIFI tools helper functions */

#ifndef _DIFITOOLS_H
#define _DIFITOOLS_H

#define SIZE 10007
#define DIFI_DATA_PACKET_SIZE 10007

// VRT
#include <vrt/vrt_init.h>
#include <vrt/vrt_string.h>
#include <vrt/vrt_types.h>
#include <vrt/vrt_util.h>
#include <vrt/vrt_write.h>
#include <vrt/vrt_read.h>

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