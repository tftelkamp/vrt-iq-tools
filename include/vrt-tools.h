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

#include <boost/format.hpp>
#include <chrono>

// VRT
#include <vrt/vrt_init.h>
#include <vrt/vrt_string.h>
#include <vrt/vrt_types.h>
#include <vrt/vrt_util.h>
#include <vrt/vrt_write.h>
#include <vrt/vrt_read.h>

#include <cctype>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <sys/time.h>
#include <time.h>

/* Timestamps
 *
 * VRT real-time fractional timestamps (VRT_TSF_REAL_TIME) count picoseconds
 * since the last second shift, and must stay below 1e12. Keep that integer
 * picosecond count as the internal currency everywhere and treat the ISO 8601
 * text form as a presentation layer only.
 *
 */
#define VRT_PS_PER_SECOND 1000000000000ULL

/* Seconds since the epoch plus a picosecond fraction. */
struct vrt_time_ps {
    int64_t  seconds;
    uint64_t frac_ps;
};

/* Carry any whole seconds out of the fraction, so frac_ps < 1e12 and the value
 * is within the bounds vrt_write_packet() enforces for VRT_TSF_REAL_TIME. */
inline void vrt_time_normalize(struct vrt_time_ps* t) {
    if (t->frac_ps >= VRT_PS_PER_SECOND) {
        t->seconds += (int64_t)(t->frac_ps / VRT_PS_PER_SECOND);
        t->frac_ps %= VRT_PS_PER_SECOND;
    }
}

/* Format as an ISO 8601 extended timestamp with nanosecond resolution, e.g.
 * "2026-08-25T12:00:00.123456789". Times are UTC; no zone designator is
 * appended, matching what earlier versions wrote.
 *
 * The fraction is truncated rather than rounded. Rounding here used to be able
 * to carry into a tenth digit: the old "%s.%06.0f" formatting of a fraction in
 * the last half microsecond of a second produced ".1000000", a malformed
 * timestamp a second adrift. */
inline std::string vrt_iso_datetime_ns(uint64_t seconds, uint64_t frac_ps) {
    struct vrt_time_ps t = {(int64_t)seconds, frac_ps};
    vrt_time_normalize(&t);

    const time_t secs = (time_t)t.seconds;
    struct tm     utc  = {};
    gmtime_r(&secs, &utc);

    char buf[48];
    snprintf(buf, sizeof(buf), "%04d-%02d-%02dT%02d:%02d:%02d.%09llu",
             utc.tm_year + 1900, utc.tm_mon + 1, utc.tm_mday,
             utc.tm_hour, utc.tm_min, utc.tm_sec,
             (unsigned long long)(t.frac_ps / 1000));
    return std::string(buf);
}

/* Days since 1970-01-01 for a civil date, proleptic Gregorian. Used instead of
 * timegm(), which is not standard C, and unlike mktime() needs no time zone. */
inline int64_t vrt_days_from_civil(int64_t y, unsigned m, unsigned d) {
    y -= m <= 2;
    const int64_t  era = (y >= 0 ? y : y - 399) / 400;
    const unsigned yoe = (unsigned)(y - era * 400);                       /* [0, 399]    */
    const unsigned doy = (153 * (m + (m > 2 ? -3 : 9)) + 2) / 5 + d - 1;  /* [0, 365]    */
    const unsigned doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;           /* [0, 146096] */
    return era * 146097 + (int64_t)doe - 719468;
}

/* Days in a month of a proleptic Gregorian year. */
inline unsigned vrt_days_in_month(int64_t y, unsigned m) {
    static const unsigned len[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    if (m == 2 and ((y % 4 == 0 and y % 100 != 0) or y % 400 == 0))
        return 29;
    return len[m - 1];
}

/* Parse "YYYY-MM-DDTHH:MM:SS[.frac]" into whole seconds and a picosecond
 * fraction. A space is accepted in place of the 'T' and a trailing 'Z' is
 * ignored. Any number of fractional digits is taken: the fraction is padded or
 * truncated to picoseconds, so the microsecond timestamps written by earlier
 * versions and the nanosecond ones written now both read back exactly.
 * Returns false if the date and time part cannot be parsed. */
inline bool vrt_parse_iso_datetime_ns(const std::string& text, struct vrt_time_ps* t) {
    std::string s = text;

    /* Times are UTC; drop a zone designator if one is present. */
    if (not s.empty() and (s.back() == 'Z' or s.back() == 'z'))
        s.pop_back();

    /* Split the fraction off: the civil-time part below carries whole seconds
     * only, and the fraction is kept to picoseconds. */
    std::string frac_str;
    const size_t dot = s.find('.');
    if (dot != std::string::npos) {
        frac_str = s.substr(dot + 1);
        s.erase(dot);
    }

    /* "YYYY-MM-DDTHH:MM:SS", or the same with the seconds left off. */
    int  y = 0, mo = 0, d = 0, h = 0, mi = 0, sec = 0, used = 0;
    char sep = 0;
    if (sscanf(s.c_str(), "%d-%d-%d%c%d:%d:%d%n", &y, &mo, &d, &sep, &h, &mi, &sec, &used) < 7) {
        used = 0;
        if (sscanf(s.c_str(), "%d-%d-%d%c%d:%d%n", &y, &mo, &d, &sep, &h, &mi, &used) < 6)
            return false;
        sec = 0;
    }
    if (sep != 'T' and sep != ' ')
        return false;

    /* Reject trailing text rather than silently accepting a malformed time. */
    for (size_t i = (size_t)used; i < s.size(); i++)
        if (not isspace((unsigned char)s[i]))
            return false;

    if (mo < 1 or mo > 12 or d < 1 or (unsigned)d > vrt_days_in_month(y, (unsigned)mo) or
        h < 0 or h > 23 or mi < 0 or mi > 59 or sec < 0 or sec > 60)
        return false;

    t->seconds = vrt_days_from_civil(y, (unsigned)mo, (unsigned)d) * 86400
                 + (int64_t)h * 3600 + (int64_t)mi * 60 + (int64_t)sec;

    /* Keep the leading digits only, then pad or truncate to 12 of them. */
    size_t digits = 0;
    while (digits < frac_str.size() and isdigit((unsigned char)frac_str[digits]))
        digits++;

    uint64_t ps = 0;
    for (size_t i = 0; i < 12; i++) {
        ps *= 10;
        if (i < digits)
            ps += (uint64_t)(frac_str[i] - '0');
    }
    t->frac_ps = ps;

    return true;
}

/* Parse either a decimal Unix time or an ISO 8601 timestamp, as the
 * --start-time option of several of these tools accepts.
 *
 * Note that a Unix time given as a decimal number cannot carry nanoseconds: a
 * double holds about 0.2 us of resolution at present-day epochs. Use the ISO
 * form when the fraction matters. */
inline bool vrt_parse_time_arg(const std::string& text, struct vrt_time_ps* t) {
    /* A bare decimal number is a Unix time; anything else is tried as ISO. */
    const char* begin = text.c_str();
    char*       end   = NULL;
    const double unix_start = strtod(begin, &end);
    if (end != begin and *end == '\0' and std::isfinite(unix_start)) {
        const double whole = std::floor(unix_start);
        t->seconds         = (int64_t)whole;
        t->frac_ps         = (uint64_t)std::llround((unix_start - whole) * (double)VRT_PS_PER_SECOND);
        vrt_time_normalize(t);
        return true;
    }
    return vrt_parse_iso_datetime_ns(text, t);
}

/* Order two picosecond timestamps. */
inline bool vrt_time_before(const struct vrt_time_ps& a, const struct vrt_time_ps& b) {
    if (a.seconds != b.seconds)
        return a.seconds < b.seconds;
    return a.frac_ps < b.frac_ps;
}

/* Signed difference a - b in picoseconds, for comparing two streams that are
 * meant to carry the same instant. Do not compute such a difference from
 * timestamps converted to a double holding seconds since the epoch: at
 * present-day epochs a double resolves about 240 ns, which is coarser than a
 * sample period above ~4 Msps.
 *
 * Differences of more than a second are clamped to +/- 2 seconds rather than
 * overflowing the multiply, so the result is only exact for nearby times. The
 * clamp keeps the sign, and its magnitude stays safe to negate. */
inline int64_t vrt_time_diff_ps(const struct vrt_time_ps& a, const struct vrt_time_ps& b) {
    const int64_t ds = a.seconds - b.seconds;
    if (ds > 1 or ds < -1)
        return ds > 0 ? 2 * (int64_t)VRT_PS_PER_SECOND : -2 * (int64_t)VRT_PS_PER_SECOND;
    return ds * (int64_t)VRT_PS_PER_SECOND + ((int64_t)a.frac_ps - (int64_t)b.frac_ps);
}

/* Time of an absolute sample index, given the time of sample 0.
 *
 * Integer arithmetic, so an exact sample rate gives an exact answer however
 * long the stream runs. Deriving each packet's fraction from a running double
 * seconds count, as these tools used to, quantized every timestamp to a whole
 * microsecond. */
inline struct vrt_time_ps vrt_time_add_samples(struct vrt_time_ps t0, uint64_t sample, double rate) {
    struct vrt_time_ps t   = t0;
    const uint64_t rate_i  = (uint64_t)std::llround(rate);

    if (rate_i != 0 and (double)rate_i == rate) {
        t.seconds += (int64_t)(sample / rate_i);
        /* 128-bit intermediate: (sample % rate_i) * 1e12 overflows 64 bits
         * for rates above about 18 Msps. */
        t.frac_ps += (uint64_t)(((unsigned __int128)(sample % rate_i) * VRT_PS_PER_SECOND) / rate_i);
    } else {
        /* Non-integer sample rate: fall back to floating point. */
        const double offset = (double)sample / rate;
        const double whole  = std::floor(offset);
        t.seconds += (int64_t)whole;
        t.frac_ps += (uint64_t)std::llround((offset - whole) * (double)VRT_PS_PER_SECOND);
    }

    vrt_time_normalize(&t);
    return t;
}

/* Current wall-clock time in picoseconds. clock_gettime() reports
 * nanoseconds, where gettimeofday() stops at microseconds. */
inline struct vrt_time_ps vrt_time_now(void) {
    struct timespec ts {};
    clock_gettime(CLOCK_REALTIME, &ts);
    struct vrt_time_ps t = {(int64_t)ts.tv_sec, (uint64_t)ts.tv_nsec * 1000ULL};
    vrt_time_normalize(&t);
    return t;
}

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
    int64_t timestamp_adjustment;
    bool has_payload_format;
    uint8_t data_item_size;
    uint8_t item_packing_field_size;
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
    context->timestamp_adjustment = 0;
    context->has_payload_format = false;
    context->data_item_size = 16;
    context->item_packing_field_size = 32;
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
    if (vrt_context->has_payload_format)
        printf("#    Sample size [bits per component]: %u\n", vrt_context->data_item_size);
    printf("#    Ref lock: %s\n", vrt_context->reflock == 1 ? "external" : "internal");
    printf("#    Time cal: %s\n", vrt_context->time_cal == 1? "pps" : "internal");
    if (vrt_context->timestamp_calibration_time != 0)
        printf("#    Cal time: %u\n", vrt_context->timestamp_calibration_time);
    if (vrt_context->timestamp_adjustment != 0)
        printf("#    Timestamp adjust: %.9f\n", (double)vrt_context->timestamp_adjustment/1e12);

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

            if (c.has.timestamp_adjustment)
                vrt_context->timestamp_adjustment = c.timestamp_adjustment;

            if (c.has.data_packet_payload_format) {
                vrt_context->has_payload_format = true;
                vrt_context->data_item_size = c.data_packet_payload_format.data_item_size + 1;
                vrt_context->item_packing_field_size = c.data_packet_payload_format.item_packing_field_size + 1;
            }

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
