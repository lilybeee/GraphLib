/*
 * png_writer.c — see png_writer.h.
 *
 * PNG = signature + chunks (IHDR, IDAT, IEND), each with a CRC-32.  IDAT holds
 * a zlib stream: 2 header bytes, a DEFLATE stream, and an Adler-32 checksum of
 * the uncompressed data.  We emit one DEFLATE block that uses the fixed
 * Huffman code (no tables to build) and a simple greedy LZ77 matcher.
 *
 * C note: DEFLATE packs bits starting at the LEAST significant bit of each
 * byte, but Huffman codes are defined MSB-first, so codes are bit-reversed
 * before being written.
 */

#include "png_writer.h"

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* ---- checksums ---- */

static uint32_t crc_table[256];
static bool     crc_ready = false;

static void crc_init(void)
{
    for (uint32_t n = 0; n < 256; ++n) {
        uint32_t c = n;
        for (int k = 0; k < 8; ++k) {
            c = (c & 1u) ? 0xEDB88320u ^ (c >> 1) : c >> 1;
        }
        crc_table[n] = c;
    }
    crc_ready = true;
}

static uint32_t crc_update(uint32_t crc, const unsigned char *p, size_t n)
{
    if (!crc_ready) {
        crc_init();
    }
    uint32_t c = crc ^ 0xFFFFFFFFu;
    for (size_t i = 0; i < n; ++i) {
        c = crc_table[(c ^ p[i]) & 0xFFu] ^ (c >> 8);
    }
    return c ^ 0xFFFFFFFFu;
}

static uint32_t adler32(const unsigned char *p, size_t n)
{
    uint32_t a = 1, b = 0;
    while (n > 0) {
        const size_t chunk = n < 5552 ? n : 5552;      /* largest run before the mod can overflow */
        for (size_t i = 0; i < chunk; ++i) {
            a += p[i];
            b += a;
        }
        a %= 65521u;
        b %= 65521u;
        p += chunk;
        n -= chunk;
    }
    return (b << 16) | a;
}

/* ---- growable byte buffer + LSB-first bit writer ---- */

typedef struct {
    unsigned char *data;
    size_t         len, cap;
    uint32_t       bits;      /* pending bits, LSB first */
    int            nbits;
    bool           failed;
} Out;

static void out_byte(Out *o, unsigned char b)
{
    if (o->len == o->cap) {
        const size_t cap = o->cap ? o->cap * 2 : 1 << 16;
        unsigned char *p = realloc(o->data, cap);
        if (p == NULL) {
            o->failed = true;
            return;
        }
        o->data = p;
        o->cap  = cap;
    }
    o->data[o->len++] = b;
}

static void put_bits(Out *o, uint32_t value, int count)       /* value's low bits, LSB first */
{
    o->bits |= value << o->nbits;
    o->nbits += count;
    while (o->nbits >= 8) {
        out_byte(o, (unsigned char)(o->bits & 0xFFu));
        o->bits >>= 8;
        o->nbits -= 8;
    }
}

static uint32_t reverse_bits(uint32_t v, int n)
{
    uint32_t r = 0;
    for (int i = 0; i < n; ++i) {
        r = (r << 1) | ((v >> i) & 1u);
    }
    return r;
}

static void put_literal_code(Out *o, unsigned sym)            /* fixed Huffman, symbols 0..287 */
{
    if (sym < 144)      { put_bits(o, reverse_bits(0x30u + sym, 8), 8); }
    else if (sym < 256) { put_bits(o, reverse_bits(0x190u + (sym - 144), 9), 9); }
    else if (sym < 280) { put_bits(o, reverse_bits(sym - 256, 7), 7); }
    else                { put_bits(o, reverse_bits(0xC0u + (sym - 280), 8), 8); }
}

static const uint16_t LEN_BASE[29]  = { 3,4,5,6,7,8,9,10,11,13,15,17,19,23,27,31,35,43,51,59,67,83,99,115,131,163,195,227,258 };
static const uint8_t  LEN_EXTRA[29] = { 0,0,0,0,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,0 };
static const uint16_t DIST_BASE[30] = { 1,2,3,4,5,7,9,13,17,25,33,49,65,97,129,193,257,385,513,769,1025,1537,2049,3073,4097,6145,8193,12289,16385,24577 };
static const uint8_t  DIST_EXTRA[30]= { 0,0,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13 };

static void put_match(Out *o, unsigned len, unsigned dist)
{
    int lc = 28;
    while (LEN_BASE[lc] > len) { --lc; }
    put_literal_code(o, 257u + (unsigned)lc);
    put_bits(o, len - LEN_BASE[lc], LEN_EXTRA[lc]);
    int dc = 29;
    while (DIST_BASE[dc] > dist) { --dc; }
    put_bits(o, reverse_bits((uint32_t)dc, 5), 5);
    put_bits(o, dist - DIST_BASE[dc], DIST_EXTRA[dc]);
}

#define HASH_BITS 15
#define WINDOW    32768u

static void deflate_fixed(Out *o, const unsigned char *in, size_t n)
{
    int32_t *head = malloc(sizeof(int32_t) << HASH_BITS);
    if (head == NULL) {
        o->failed = true;
        return;
    }
    for (size_t i = 0; i < (1u << HASH_BITS); ++i) {
        head[i] = -1;
    }
    put_bits(o, 1, 1);          /* BFINAL */
    put_bits(o, 1, 2);          /* BTYPE = 01: fixed Huffman */
    size_t i = 0;
    while (i < n && !o->failed) {
        unsigned best_len = 0, best_dist = 0;
        if (i + 3 <= n) {
            const uint32_t h = (((uint32_t)in[i] << 10) ^ ((uint32_t)in[i + 1] << 5) ^ in[i + 2])
                               & ((1u << HASH_BITS) - 1u);
            const int32_t cand = head[h];
            head[h] = (int32_t)i;
            if (cand >= 0 && i - (size_t)cand <= WINDOW) {
                unsigned l = 0;
                const size_t max = n - i < 258 ? n - i : 258;
                while (l < max && in[(size_t)cand + l] == in[i + l]) {
                    ++l;
                }
                if (l >= 3) {
                    best_len  = l;
                    best_dist = (unsigned)(i - (size_t)cand);
                }
            }
        }
        if (best_len >= 3) {
            put_match(o, best_len, best_dist);
            /* Index the skipped positions too (cheaply) so later matches can find them. */
            for (size_t k = 1; k < best_len && i + k + 3 <= n; ++k) {
                const size_t p = i + k;
                head[(((uint32_t)in[p] << 10) ^ ((uint32_t)in[p + 1] << 5) ^ in[p + 2])
                     & ((1u << HASH_BITS) - 1u)] = (int32_t)p;
            }
            i += best_len;
        } else {
            put_literal_code(o, in[i]);
            ++i;
        }
    }
    put_literal_code(o, 256);   /* end of block */
    if (o->nbits > 0) {
        out_byte(o, (unsigned char)(o->bits & 0xFFu));
        o->bits = 0;
        o->nbits = 0;
    }
    free(head);
}

/* ---- PNG container ---- */

static bool write_be32(FILE *f, uint32_t v)
{
    const unsigned char b[4] = { (unsigned char)(v >> 24), (unsigned char)(v >> 16),
                                 (unsigned char)(v >> 8), (unsigned char)v };
    return fwrite(b, 1, 4, f) == 4;
}

static bool write_chunk(FILE *f, const char type[4], const unsigned char *data, size_t len)
{
    uint32_t crc = crc_update(0, (const unsigned char *)type, 4);
    crc = crc_update(crc, data, len);
    return write_be32(f, (uint32_t)len) && fwrite(type, 1, 4, f) == 4 &&
           (len == 0 || fwrite(data, 1, len, f) == len) && write_be32(f, crc);
}

GraphStatus png_write_rgb(FILE *out, const unsigned char *rgb, unsigned w, unsigned h)
{
    if (out == NULL || rgb == NULL || w == 0 || h == 0) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t stride = (size_t)w * 3;
    if (h > (SIZE_MAX / 2) / (stride + 1) || (stride + 1) * h > 0xFFFFFFF0u) {
        return GRAPH_ERR_TOO_LARGE;
    }
    /* Filter every row with "Sub" (type 1): each byte minus the byte one pixel to the left. */
    const size_t raw_len = (stride + 1) * h;
    unsigned char *raw = malloc(raw_len);
    if (raw == NULL) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    for (unsigned y = 0; y < h; ++y) {
        unsigned char       *dst = raw + (size_t)y * (stride + 1);
        const unsigned char *src = rgb + (size_t)y * stride;
        dst[0] = 1;
        for (size_t x = 0; x < stride; ++x) {
            dst[1 + x] = (unsigned char)(src[x] - (x >= 3 ? src[x - 3] : 0));
        }
    }

    Out z = { 0 };
    out_byte(&z, 0x78);         /* zlib header: deflate, 32K window */
    out_byte(&z, 0x01);
    deflate_fixed(&z, raw, raw_len);
    const uint32_t ad = adler32(raw, raw_len);
    free(raw);
    out_byte(&z, (unsigned char)(ad >> 24));
    out_byte(&z, (unsigned char)(ad >> 16));
    out_byte(&z, (unsigned char)(ad >> 8));
    out_byte(&z, (unsigned char)ad);
    if (z.failed) {
        free(z.data);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }

    static const unsigned char sig[8] = { 0x89, 'P', 'N', 'G', '\r', '\n', 0x1A, '\n' };
    unsigned char ihdr[13] = { (unsigned char)(w >> 24), (unsigned char)(w >> 16),
                               (unsigned char)(w >> 8), (unsigned char)w,
                               (unsigned char)(h >> 24), (unsigned char)(h >> 16),
                               (unsigned char)(h >> 8), (unsigned char)h,
                               8, 2, 0, 0, 0 };   /* 8-bit, colour type 2 (RGB) */
    const bool ok = fwrite(sig, 1, 8, out) == 8 && write_chunk(out, "IHDR", ihdr, 13) &&
                    write_chunk(out, "IDAT", z.data, z.len) && write_chunk(out, "IEND", NULL, 0);
    free(z.data);
    return ok ? GRAPH_OK : GRAPH_ERR_IO;
}
