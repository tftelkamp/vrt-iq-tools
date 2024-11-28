#ifndef _TRACKEREXTENDEDCONTEXT_H
#define _TRACKEREXTENDEDCONTEXT_H

struct tracker_ext_context_type {
    bool tracker_ext_context_received = false;

    char object_name[32] = "";
    char tracking_source[32] = "";
    int32_t object_id = NAN;
    float azimuth = NAN;
    float elevation = NAN;
    float ra = NAN;
    float dec = NAN;
    double distance = NAN;
    double speed = NAN;
    double frequency = NAN;
    double doppler = NAN;
    double doppler_rate = NAN;

    uint32_t stream_id;
    uint64_t fractional_seconds_timestamp;
    uint64_t integer_seconds_timestamp;
};

bool tracker_process(uint32_t* buffer, uint32_t size, packet_type* vrt_packet, tracker_ext_context_type* tracker_ext_context) {

    if (vrt_packet->oui == 0xFF0043) { // add information and packet class checks

        tracker_ext_context->stream_id = vrt_packet->stream_id;
        tracker_ext_context->fractional_seconds_timestamp = vrt_packet->fractional_seconds_timestamp;
        tracker_ext_context->integer_seconds_timestamp = vrt_packet->integer_seconds_timestamp;

        memcpy(tracker_ext_context->object_name, (char*)&buffer[vrt_packet->offset], sizeof(tracker_ext_context->object_name));
        memcpy(tracker_ext_context->tracking_source, (char*)&buffer[vrt_packet->offset+8], sizeof(tracker_ext_context->tracking_source));
        memcpy(&tracker_ext_context->object_id, (char*)&buffer[vrt_packet->offset+16], sizeof(int32_t));

        memcpy(&tracker_ext_context->azimuth, (char*)&buffer[vrt_packet->offset+17], sizeof(float));
        memcpy(&tracker_ext_context->elevation, (char*)&buffer[vrt_packet->offset+18], sizeof(float));
        memcpy(&tracker_ext_context->ra, (char*)&buffer[vrt_packet->offset+19], sizeof(float));
        memcpy(&tracker_ext_context->dec, (char*)&buffer[vrt_packet->offset+20], sizeof(float));

        memcpy(&tracker_ext_context->distance, (char*)&buffer[vrt_packet->offset+21], sizeof(double));
        memcpy(&tracker_ext_context->speed, (char*)&buffer[vrt_packet->offset+23], sizeof(double));
        memcpy(&tracker_ext_context->frequency, (char*)&buffer[vrt_packet->offset+25], sizeof(double));
        memcpy(&tracker_ext_context->doppler, (char*)&buffer[vrt_packet->offset+27], sizeof(double));
        memcpy(&tracker_ext_context->doppler_rate, (char*)&buffer[vrt_packet->offset+29], sizeof(double));

        tracker_ext_context->tracker_ext_context_received = true;
        return true;
    } else
        return false;
}

#endif
