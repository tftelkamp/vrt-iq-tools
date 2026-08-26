/* VDIF helper functions
 *
 * Frame header layout, Data Array format and sample representation follow the
 * VLBI Data Interchange Format (VDIF) specification release 1.1.1 (June 2014),
 * sections 6, 9 and 10. Section and note numbers below refer to that document.
 */

#ifndef _VDIF_H
#define _VDIF_H

#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

#include <string>

/* Header of a Data Frame. Legacy (16 byte) headers are discouraged by the
 * specification (Note 1) and are not written here. */
#define VDIF_HEADER_BYTES 32
#define VDIF_HEADER_WORDS 8

/* VDIF version number in Word 2, zero for this release (Note 3). */
#define VDIF_VERSION 0

/* A Data Frame, header included, is at most 2^27 bytes (Note 5). */
#define VDIF_MAX_FRAME_BYTES (1 << 27)

/* Data Array size aimed for when it is not given on the command line. Frames of
 * 8000 byte payload plus header are the common choice in VLBI recordings. */
#define VDIF_DEFAULT_PAYLOAD_BYTES 8000

struct vdif_header {
    uint32_t word[VDIF_HEADER_WORDS];
};

/* Words 2 and 3 hold information that is static for a Data Thread, so they are
 * written once. Words 4-7 are the extended user data, all zero for EDV 0
 * (Note 9). For complex data, bits is the number of bits in one component, not
 * in the complete sample (Note 7). */
static inline void vdif_init_header(struct vdif_header* h,
                                    uint32_t            payload_bytes,
                                    uint32_t            bits,
                                    bool                is_complex,
                                    uint32_t            log2_channels,
                                    uint16_t            station_id,
                                    uint16_t            thread_id)
{
    memset(h, 0, sizeof(*h));

    h->word[2] = ((uint32_t)VDIF_VERSION << 29) | ((log2_channels & 0x1F) << 24) |
                 (((payload_bytes + VDIF_HEADER_BYTES) / 8) & 0x00FFFFFF);

    h->word[3] = ((uint32_t)(is_complex ? 1 : 0) << 31) | (((bits - 1) & 0x1F) << 26) |
                 (((uint32_t)thread_id & 0x3FF) << 16) | station_id;
}

/* The reference epoch in Word 1 is static as well, but is only known once the
 * first timestamp has been seen. */
static inline void vdif_set_epoch(struct vdif_header* h, uint8_t epoch)
{
    h->word[1] = (h->word[1] & 0x00FFFFFF) | (((uint32_t)epoch & 0x3F) << 24);
}

/* Words 0 and 1 are the only ones that change per frame. Bits 31-30 of Word 1
 * are unassigned and are left at zero. */
static inline void vdif_set_time(struct vdif_header* h, uint32_t seconds, uint32_t frame, bool invalid)
{
    h->word[0] = (invalid ? (1u << 31) : 0u) | (seconds & 0x3FFFFFFF);
    h->word[1] = (h->word[1] & 0x3F000000) | (frame & 0x00FFFFFF);
}

/* The reference epoch is the 6-month period in which the clock was set, with
 * zero being the first half of 2000 (Note 2a). Returns the epoch number and the
 * Unix time at which it starts, from which the second count of Word 0 follows
 * as a difference. Note that the seconds field counts all seconds including
 * leap seconds, where Unix time repeats one, so the two only agree as long as
 * no leap second falls between the epoch and the recording. */
static inline bool vdif_epoch_from_unix(int64_t t, uint8_t* epoch, int64_t* epoch_start)
{
    const time_t tt = (time_t)t;
    struct tm    g;

    if (gmtime_r(&tt, &g) == NULL)
        return false;

    const int year = g.tm_year + 1900;
    const int half = (g.tm_mon >= 6) ? 1 : 0;

    if (year < 2000)
        return false;

    struct tm e;
    memset(&e, 0, sizeof(e));
    e.tm_year = g.tm_year;
    e.tm_mon  = half ? 6 : 0;
    e.tm_mday = 1;

    const time_t start = timegm(&e);
    if (start == (time_t)-1)
        return false;

    /* Wraps every 32 years */
    *epoch       = (uint8_t)(((year - 2000) * 2 + half) & 0x3F);
    *epoch_start = (int64_t)start;

    return true;
}

static inline uint32_t vdif_gcd(uint32_t a, uint32_t b)
{
    while (b != 0) {
        const uint32_t t = a % b;
        a                = b;
        b                = t;
    }
    return a;
}

/* Pick the size of the Data Array. The number of samples in a frame has to
 * divide the sample rate, so that there is an integral number of Data Frames
 * per second (Word 1, bits 23-0), and has to fill an even number of 32-bit
 * words (section 9.1, rule 1). Of the sizes that qualify, the largest payload
 * not above the target is taken, or else the smallest one above it.
 *
 * bits_per_sample is the length of a complete sample, so twice the component
 * size for complex data. */
static inline bool vdif_solve_geometry(uint32_t  sample_rate,
                                       uint32_t  bits_per_sample,
                                       uint32_t  target_payload_bytes,
                                       uint32_t* payload_bytes,
                                       uint32_t* samples_per_frame,
                                       uint32_t* frames_per_second)
{
    if (sample_rate == 0 or bits_per_sample == 0)
        return false;

    /* Samples per frame must be a multiple of this to end on an even word */
    const uint32_t granularity = 64 / vdif_gcd(64, bits_per_sample);

    uint32_t best     = 0;
    uint32_t best_alt = 0;

    for (uint32_t i = 1; (uint64_t)i * i <= (uint64_t)sample_rate; i++) {

        if (sample_rate % i != 0)
            continue;

        const uint32_t divisors[2] = {i, sample_rate / i};

        for (int d = 0; d < 2; d++) {

            const uint32_t spf = divisors[d];

            if (spf % granularity != 0)
                continue;

            const uint64_t bytes = (uint64_t)spf * bits_per_sample / 8;

            if (bytes < 8 or bytes > (uint64_t)(VDIF_MAX_FRAME_BYTES - VDIF_HEADER_BYTES))
                continue;

            if (bytes <= target_payload_bytes) {
                if (spf > best)
                    best = spf;
            } else {
                if (best_alt == 0 or spf < best_alt)
                    best_alt = spf;
            }
        }
    }

    if (best == 0)
        best = best_alt;

    if (best == 0)
        return false;

    *samples_per_frame = best;
    *payload_bytes     = (uint32_t)((uint64_t)best * bits_per_sample / 8);
    *frames_per_second = sample_rate / best;

    return true;
}

/* Read n bits (at most 32) at an arbitrary bit position of a packed stream. The
 * second word is only touched when the field really extends into it. */
static inline uint32_t vdif_get_bits(const uint32_t* src, size_t bit, unsigned n)
{
    const size_t   w  = bit / 32;
    const unsigned sh = (unsigned)(bit % 32);

    uint32_t v = src[w] >> sh;

    if (sh + n > 32)
        v |= src[w + 1] << (32 - sh);

    return (n == 32) ? v : (v & ((1u << n) - 1));
}

/* Insert n bits (at most 32) at an arbitrary bit position. The destination is
 * expected to be zeroed. */
static inline void vdif_put_bits(uint32_t* dst, size_t bit, uint32_t v, unsigned n)
{
    const size_t   w  = bit / 32;
    const unsigned sh = (unsigned)(bit % 32);

    dst[w] |= v << sh;

    if (sh + n > 32)
        dst[w + 1] |= v >> (32 - sh);
}

/* Copy a run of bits between two packed streams. Both VRT and VDIF hold the
 * oldest sample in the least significant bits of the first word (section 9.1,
 * rule 3), so this moves samples without touching their values. Whole words are
 * copied directly when both ends happen to be word aligned, which is the common
 * case. The destination is expected to be zeroed. */
static inline void vdif_copy_bits(uint32_t*       dst,
                                  size_t          dst_bit,
                                  const uint32_t* src,
                                  size_t          src_bit,
                                  size_t          nbits)
{
    if (dst_bit % 32 == 0 and src_bit % 32 == 0) {
        const size_t words = nbits / 32;
        memcpy(&dst[dst_bit / 32], &src[src_bit / 32], words * sizeof(uint32_t));
        dst_bit += words * 32;
        src_bit += words * 32;
        nbits -= words * 32;
    }

    while (nbits > 0) {
        const unsigned n = (nbits >= 32) ? 32 : (unsigned)nbits;
        vdif_put_bits(dst, dst_bit, vdif_get_bits(src, src_bit, n), n);
        dst_bit += n;
        src_bit += n;
        nbits -= n;
    }
}

/* Station ID is either the globally assigned 2-character ASCII identifier or an
 * unsigned 16-bit number, the two being told apart by the first character being
 * below the ASCII value of '0' (Note 8). Something that parses as a number is
 * taken as one, so that a numeric identifier is never turned into two
 * characters by accident. */
static inline bool vdif_parse_station_id(const std::string& s, uint16_t* station_id, bool* is_ascii)
{
    try {
        size_t              end = 0;
        const unsigned long v   = std::stoul(s, &end, 0);
        if (end == s.size()) {
            if (v > 0xFFFF)
                return false;
            *station_id = (uint16_t)v;
            *is_ascii   = false;
            return true;
        }
    } catch (...) {
        /* Not a number, try the ASCII form below */
    }

    if (s.size() == 2 and (unsigned char)s[0] >= '0') {
        *station_id = (uint16_t)(((unsigned char)s[0] << 8) | (unsigned char)s[1]);
        *is_ascii   = true;
        return true;
    }

    return false;
}

#endif
