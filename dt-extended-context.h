#ifndef _DTEXTENDEDCONTEXT_H
#define _DTEXTENDEDCONTEXT_H

struct dt_ext_context_type {
	bool dt_ext_context_received = false;
	float azimuth = NAN;
	float elevation = NAN;
	float azimuth_error = NAN;
	float elevation_error = NAN;
	float azimuth_speed = NAN;
	float elevation_speed = NAN;
	float focusbox = NAN;
	uint32_t stream_id;
	uint64_t fractional_seconds_timestamp;
    uint64_t integer_seconds_timestamp;
};

bool dt_process(uint32_t* buffer, uint32_t size, packet_type* vrt_packet, dt_ext_context_type* dt_ext_context) {

	if (vrt_packet->oui == 0xFF0042) { // add information and packet class checks

		dt_ext_context->stream_id = vrt_packet->stream_id;
		dt_ext_context->fractional_seconds_timestamp = vrt_packet->fractional_seconds_timestamp;
		dt_ext_context->integer_seconds_timestamp = vrt_packet->integer_seconds_timestamp;

		memcpy(&dt_ext_context->azimuth, (char*)&buffer[vrt_packet->offset], sizeof(float));
		memcpy(&dt_ext_context->elevation, (char*)&buffer[vrt_packet->offset+1], sizeof(float));
		memcpy(&dt_ext_context->azimuth_error, (char*)&buffer[vrt_packet->offset+2], sizeof(float));
		memcpy(&dt_ext_context->elevation_error, (char*)&buffer[vrt_packet->offset+3], sizeof(float));
		memcpy(&dt_ext_context->azimuth_speed, (char*)&buffer[vrt_packet->offset+4], sizeof(float));
		memcpy(&dt_ext_context->elevation_speed, (char*)&buffer[vrt_packet->offset+5], sizeof(float));
		memcpy(&dt_ext_context->focusbox, (char*)&buffer[vrt_packet->offset+6], sizeof(float));

		dt_ext_context->dt_ext_context_received = true;
		return true;
	} else
		return false;
}


#endif