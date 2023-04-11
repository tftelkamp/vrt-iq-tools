#ifndef _DTEXTENDEDCONTEXT_H
#define _DTEXTENDEDCONTEXT_H

struct dt_ext_context_type {
	bool dt_ext_context_received = false;

	uint32_t active_tracker;
	bool tracking_enabled = false;
	bool refraction = false;
	bool dt_model = false;
	bool refraction_j2000 = false;
	bool dt_model_j2000 = false;

	float azimuth = NAN;
	float elevation = NAN;
	float azimuth_error = NAN;
	float elevation_error = NAN;
	float azimuth_speed = NAN;
	float elevation_speed = NAN;
	float elevation_offset = NAN;
	float azimuth_offset = NAN;
	float focusbox = NAN;
	float ra_setpoint = NAN;
	float dec_setpoint = NAN;
	float ra_current = NAN;
	float dec_current = NAN;
	float model_a0 = NAN;
	float model_c1 = NAN;
	float model_c2 = NAN;
	float model_e0 = NAN;
	float model_b = NAN;
	float model_za = NAN;
	float model_aa = NAN;

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
		memcpy(&dt_ext_context->azimuth_offset, (char*)&buffer[vrt_packet->offset+8], sizeof(float));
		memcpy(&dt_ext_context->elevation_offset, (char*)&buffer[vrt_packet->offset+9], sizeof(float));
		memcpy(&dt_ext_context->ra_setpoint, (char*)&buffer[vrt_packet->offset+10], sizeof(float));
		memcpy(&dt_ext_context->dec_setpoint, (char*)&buffer[vrt_packet->offset+11], sizeof(float));
		memcpy(&dt_ext_context->ra_current, (char*)&buffer[vrt_packet->offset+12], sizeof(float));
		memcpy(&dt_ext_context->dec_current, (char*)&buffer[vrt_packet->offset+13], sizeof(float));
		memcpy(&dt_ext_context->model_a0, (char*)&buffer[vrt_packet->offset+14], sizeof(float));
		memcpy(&dt_ext_context->model_c1, (char*)&buffer[vrt_packet->offset+15], sizeof(float));
		memcpy(&dt_ext_context->model_c2, (char*)&buffer[vrt_packet->offset+16], sizeof(float));
		memcpy(&dt_ext_context->model_e0, (char*)&buffer[vrt_packet->offset+17], sizeof(float));
		memcpy(&dt_ext_context->model_b, (char*)&buffer[vrt_packet->offset+18], sizeof(float));
		memcpy(&dt_ext_context->model_za, (char*)&buffer[vrt_packet->offset+19], sizeof(float));
		memcpy(&dt_ext_context->model_aa, (char*)&buffer[vrt_packet->offset+20], sizeof(float));

		dt_ext_context->active_tracker = (buffer[vrt_packet->offset+7]) & 0x0F;
		dt_ext_context->tracking_enabled = (buffer[vrt_packet->offset+7] >> 8) & (1<<7);
		dt_ext_context->refraction = (buffer[vrt_packet->offset+7] >> 8) & (1<<6);
		dt_ext_context->dt_model = (buffer[vrt_packet->offset+7] >> 8) & (1<<5);
		dt_ext_context->refraction_j2000 = (buffer[vrt_packet->offset+7] >> 8) & (1<<1);
		dt_ext_context->dt_model_j2000 = (buffer[vrt_packet->offset+7] >> 8) & (1<<0);

		dt_ext_context->dt_ext_context_received = true;
		return true;
	} else
		return false;
}

#endif