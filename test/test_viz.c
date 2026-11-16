/*
 * test_viz.c — visualization: PNG encoder round trip, SVG content, layout
 * properties, large graphs, error codes.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "graph_algorithms.h"
#include "graph_viz.h"
#include "png_writer.h"

#define CHECK(cond)                                                                  \
    do {                                                                             \
        if (!(cond)) {                                                               \
            fprintf(stderr, "%s:%d: CHECK failed: %s\n", __FILE__, __LINE__, #cond); \
            exit(1);                                                                 \
        }                                                                            \
    } while (0)

/* ------------------------------------------------------------------ */
/*  A tiny PNG decoder (only what our encoder emits: fixed Huffman)     */
/* ------------------------------------------------------------------ */

typedef struct {
    const unsigned char *p;
    size_t               len, pos;
    unsigned             bits;
    int                  nbits;
} BitIn;

static unsigned get_bits(BitIn *b, int n)
{
    while (b->nbits < n) {
        CHECK(b->pos < b->len);
        b->bits |= (unsigned)b->p[b->pos++] << b->nbits;
        b->nbits += 8;
    }
    const unsigned v = b->bits & ((1u << n) - 1u);
    b->bits >>= n;
    b->nbits -= n;
    return v;
}

static unsigned get_code_msb(BitIn *b, int n)             /* Huffman codes are read MSB-first */
{
    unsigned v = 0;
    for (int i = 0; i < n; ++i) { v = (v << 1) | get_bits(b, 1); }
    return v;
}

static unsigned decode_literal_symbol(BitIn *b)
{
    unsigned c = get_code_msb(b, 7);
    if (c <= 0x17) { return 256 + c; }                       /* 256..279 */
    c = (c << 1) | get_bits(b, 1);
    if (c >= 0x30 && c <= 0xBF) { return c - 0x30; }          /* 0..143 */
    if (c >= 0xC0 && c <= 0xC7) { return 280 + (c - 0xC0); }  /* 280..287 */
    c = (c << 1) | get_bits(b, 1);
    CHECK(c >= 0x190 && c <= 0x1FF);
    return 144 + (c - 0x190);                                 /* 144..255 */
}

static const unsigned LB[29] = { 3,4,5,6,7,8,9,10,11,13,15,17,19,23,27,31,35,43,51,59,67,83,99,115,131,163,195,227,258 };
static const unsigned LE[29] = { 0,0,0,0,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,0 };
static const unsigned DB[30] = { 1,2,3,4,5,7,9,13,17,25,33,49,65,97,129,193,257,385,513,769,1025,1537,2049,3073,4097,6145,8193,12289,16385,24577 };
static const unsigned DE[30] = { 0,0,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13 };

static unsigned char *inflate_fixed(const unsigned char *z, size_t zlen, size_t expect)
{
    CHECK(zlen > 6 && z[0] == 0x78);
    BitIn b = { z + 2, zlen - 6, 0, 0, 0 };
    unsigned char *out = malloc(expect + 1);
    size_t n = 0;
    CHECK(out != NULL);
    CHECK(get_bits(&b, 1) == 1);          /* final block */
    CHECK(get_bits(&b, 2) == 1);          /* fixed Huffman */
    for (;;) {
        const unsigned sym = decode_literal_symbol(&b);
        if (sym == 256) { break; }
        if (sym < 256) {
            CHECK(n < expect);
            out[n++] = (unsigned char)sym;
            continue;
        }
        const unsigned li = sym - 257;
        CHECK(li < 29);
        const unsigned len = LB[li] + get_bits(&b, (int)LE[li]);
        const unsigned di = get_code_msb(&b, 5);
        CHECK(di < 30);
        const unsigned dist = DB[di] + get_bits(&b, (int)DE[di]);
        CHECK(dist <= n && n + len <= expect);
        for (unsigned i = 0; i < len; ++i) { out[n] = out[n - dist]; ++n; }
    }
    CHECK(n == expect);
    return out;
}

static uint32_t be32(const unsigned char *p)
{
    return ((uint32_t)p[0] << 24) | ((uint32_t)p[1] << 16) | ((uint32_t)p[2] << 8) | p[3];
}

static uint32_t crc32_ref(const unsigned char *p, size_t n)
{
    uint32_t c = 0xFFFFFFFFu;
    for (size_t i = 0; i < n; ++i) {
        c ^= p[i];
        for (int k = 0; k < 8; ++k) { c = (c & 1u) ? 0xEDB88320u ^ (c >> 1) : c >> 1; }
    }
    return c ^ 0xFFFFFFFFu;
}

/** Reads a whole stream into memory. */
static unsigned char *slurp(FILE *f, size_t *len)
{
    fflush(f);
    fseek(f, 0, SEEK_END);
    const long n = ftell(f);
    rewind(f);
    unsigned char *buf = malloc((size_t)n + 1);
    CHECK(buf != NULL && fread(buf, 1, (size_t)n, f) == (size_t)n);
    buf[n] = 0;
    *len = (size_t)n;
    return buf;
}

/** Parses and fully verifies a PNG; returns the RGB pixels (caller frees). */
static unsigned char *decode_png(const unsigned char *d, size_t len, unsigned *w, unsigned *h)
{
    static const unsigned char sig[8] = { 0x89, 'P', 'N', 'G', '\r', '\n', 0x1A, '\n' };
    CHECK(len > 8 && memcmp(d, sig, 8) == 0);
    size_t pos = 8, idat_len = 0;
    unsigned char *idat = malloc(len);
    CHECK(idat != NULL);
    bool seen_ihdr = false, seen_iend = false;
    while (pos + 12 <= len) {
        const uint32_t clen = be32(d + pos);
        CHECK(pos + 12 + clen <= len);
        CHECK(crc32_ref(d + pos + 4, 4 + clen) == be32(d + pos + 8 + clen));   /* chunk CRC */
        if (memcmp(d + pos + 4, "IHDR", 4) == 0) {
            *w = be32(d + pos + 8);
            *h = be32(d + pos + 12);
            CHECK(d[pos + 16] == 8 && d[pos + 17] == 2);
            seen_ihdr = true;
        } else if (memcmp(d + pos + 4, "IDAT", 4) == 0) {
            memcpy(idat + idat_len, d + pos + 8, clen);
            idat_len += clen;
        } else if (memcmp(d + pos + 4, "IEND", 4) == 0) {
            seen_iend = true;
        }
        pos += 12 + clen;
    }
    CHECK(seen_ihdr && seen_iend && pos == len);
    const size_t stride = (size_t)*w * 3, raw_len = (stride + 1) * *h;
    unsigned char *raw = inflate_fixed(idat, idat_len, raw_len);
    /* Adler-32 of the raw data sits in the last four bytes of the zlib stream. */
    uint32_t a = 1, bsum = 0;
    for (size_t i = 0; i < raw_len; ++i) { a = (a + raw[i]) % 65521u; bsum = (bsum + a) % 65521u; }
    CHECK(((bsum << 16) | a) == be32(idat + idat_len - 4));
    free(idat);
    unsigned char *rgb = malloc(stride * *h);
    CHECK(rgb != NULL);
    for (unsigned y = 0; y < *h; ++y) {                       /* undo the Sub filter */
        const unsigned char *row = raw + (size_t)y * (stride + 1);
        CHECK(row[0] == 1);
        for (size_t x = 0; x < stride; ++x) {
            rgb[(size_t)y * stride + x] =
                (unsigned char)(row[1 + x] + (x >= 3 ? rgb[(size_t)y * stride + x - 3] : 0));
        }
    }
    free(raw);
    return rgb;
}

static void test_png_encoder(void)
{
    /* Mixed content: flat runs, gradients, noise, and a size that is not a multiple of anything. */
    const unsigned sizes[][2] = { { 1, 1 }, { 7, 5 }, { 64, 3 }, { 301, 200 } };
    unsigned seed = 12345;
    for (int t = 0; t < 4; ++t) {
        const unsigned w = sizes[t][0], h = sizes[t][1];
        unsigned char *img = malloc((size_t)w * h * 3);
        CHECK(img != NULL);
        for (size_t i = 0; i < (size_t)w * h * 3; ++i) {
            seed = seed * 1103515245u + 12345u;
            const size_t pix = i / 3;
            img[i] = (pix % w < w / 3) ? 255 : (pix % w < 2 * w / 3) ? (unsigned char)(pix * 7 + i % 3)
                                                                      : (unsigned char)(seed >> 16);
        }
        FILE *f = tmpfile();
        CHECK(f != NULL);
        CHECK(png_write_rgb(f, img, w, h) == GRAPH_OK);
        size_t len;
        unsigned char *file = slurp(f, &len);
        unsigned rw, rh;
        unsigned char *back = decode_png(file, len, &rw, &rh);
        CHECK(rw == w && rh == h && memcmp(back, img, (size_t)w * h * 3) == 0);
        free(file); free(back); free(img);
        fclose(f);
    }
    /* A flat image must compress to a tiny file (the point of Sub filter + LZ77). */
    unsigned char *flat = malloc(500u * 500u * 3);
    memset(flat, 255, 500u * 500u * 3);
    FILE *f = tmpfile();
    CHECK(png_write_rgb(f, flat, 500, 500) == GRAPH_OK);
    size_t len;
    unsigned char *file = slurp(f, &len);
    CHECK(len < 12000);                    /* 750 KB of raw data */
    free(file); fclose(f);
    CHECK(png_write_rgb(NULL, flat, 1, 1) == GRAPH_ERR_INVALID_ARGUMENT);
    free(flat);
}

/* ------------------------------------------------------------------ */
/*  SVG                                                                 */
/* ------------------------------------------------------------------ */

static size_t count_substr(const char *hay, const char *needle)
{
    size_t c = 0;
    for (const char *p = hay; (p = strstr(p, needle)) != NULL; p += strlen(needle)) { ++c; }
    return c;
}

static char *render_svg(Graph *g, VizOptions o)
{
    o.format = VIZ_FORMAT_SVG;
    o.legend = false;                    /* the element counts below must be exact */
    FILE *f = tmpfile();
    CHECK(f != NULL && graph_visualize_to(g, f, &o) == GRAPH_OK);
    size_t len;
    char *text = (char *)slurp(f, &len);
    fclose(f);
    return text;
}

static void test_svg(void)
{
    Graph *g = graph_gen_cycle(GRAPH_LIST, 12, false);
    CHECK(g != NULL);
    VizOptions o = graph_viz_default_options();
    char *svg = render_svg(g, o);
    CHECK(strncmp(svg, "<svg ", 5) == 0 && strstr(svg, "</svg>") != NULL);
    CHECK(count_substr(svg, "<circle") == 12);
    CHECK(count_substr(svg, "<line") == 12);
    CHECK(count_substr(svg, "<text") == 12);               /* auto labels: small graph */
    CHECK(strstr(svg, "marker-end") == NULL);               /* undirected: no arrows */
    free(svg);

    o.show_labels = 0;
    svg = render_svg(g, o);
    CHECK(count_substr(svg, "<text") == 0);
    free(svg);

    /* Custom labels are XML-escaped; weights are optional. */
    const char *names[12] = { "a<b", "x&y", "\"q\"", "d", "e", "f", "g", "h", "i", "j", "k", "l" };
    o.show_labels = 1;
    o.labels = names;
    o.title = "A & B";
    svg = render_svg(g, o);
    CHECK(strstr(svg, "a&lt;b") && strstr(svg, "x&amp;y") && strstr(svg, "&quot;q&quot;"));
    CHECK(strstr(svg, "A &amp; B") != NULL);
    CHECK(strstr(svg, "a<b") == NULL);
    free(svg);

    o = graph_viz_default_options();
    o.show_weights = true;
    o.show_labels = 0;
    svg = render_svg(g, o);
    CHECK(count_substr(svg, "<text") == 12);                /* one weight per edge */
    free(svg);
    graph_destroy(g);

    /* Directed graph: arrowheads, one line per arc; removed vertices are not drawn. */
    Graph *d = graph_gen_path(GRAPH_MATRIX, 6, true);
    CHECK(d != NULL && graph_remove_vertex(d, 5));
    svg = render_svg(d, graph_viz_default_options());
    CHECK(strstr(svg, "marker-end") != NULL);
    CHECK(count_substr(svg, "<circle") == 5 && count_substr(svg, "<line") == 4);
    free(svg);
    graph_destroy(d);

    /* An empty graph and a single vertex still produce valid output. */
    Graph *e = graph_create(false);
    svg = render_svg(e, graph_viz_default_options());
    CHECK(strstr(svg, "</svg>") && count_substr(svg, "<circle") == 0);
    free(svg);
    CHECK(graph_add_vertex(e) == 0);
    svg = render_svg(e, graph_viz_default_options());
    CHECK(count_substr(svg, "<circle") == 1);
    free(svg);
    graph_destroy(e);
}

/* ------------------------------------------------------------------ */
/*  Layout                                                              */
/* ------------------------------------------------------------------ */

static void test_layout(void)
{
    Graph *g = graph_gen_erdos_renyi(GRAPH_LIST, 150, 0.04, 3, false, 1.0, 1.0);
    CHECK(g != NULL && graph_remove_vertex(g, 7));
    const size_t n = graph_id_bound(g);
    double x[150], y[150], x2[150], y2[150], x3[150], y3[150];
    VizOptions o = graph_viz_default_options();
    CHECK(graph_layout(g, x, y, &o) == GRAPH_OK);
    for (size_t v = 0; v < n; ++v) {
        if (!graph_is_alive(g, (int)v)) {
            CHECK(isnan(x[v]) && isnan(y[v]));
        } else {
            CHECK(x[v] >= 0.0 && x[v] <= 1.0 && y[v] >= 0.0 && y[v] <= 1.0);
        }
    }
    /* Deterministic for a seed, different for another. */
    CHECK(graph_layout(g, x2, y2, &o) == GRAPH_OK);
    for (size_t v = 0; v < n; ++v) {
        if (graph_is_alive(g, (int)v)) { CHECK(x[v] == x2[v] && y[v] == y2[v]); }
    }
    o.seed = 99;
    CHECK(graph_layout(g, x3, y3, &o) == GRAPH_OK);
    bool differs = false;
    for (size_t v = 0; v < n; ++v) {
        differs = differs || (graph_is_alive(g, (int)v) && x3[v] != x[v]);
    }
    CHECK(differs);

    /* The layout is meaningful: connected pairs sit closer than random pairs. */
    double edge_sum = 0, pair_sum = 0;
    size_t edge_n = 0, pair_n = 0;
    for (int u = 0; u < (int)n; ++u) {
        if (!graph_is_alive(g, u)) { continue; }
        for (int v = u + 1; v < (int)n; ++v) {
            if (!graph_is_alive(g, v)) { continue; }
            const double d = hypot(x[u] - x[v], y[u] - y[v]);
            pair_sum += d; ++pair_n;
            if (graph_has_edge(g, u, v, NULL)) { edge_sum += d; ++edge_n; }
        }
    }
    CHECK(edge_n > 0 && edge_sum / (double)edge_n < 0.6 * pair_sum / (double)pair_n);

    /* Dense clusters: members end up much closer to each other than to other clusters. */
    Graph *c = graph_create_sized(120, false);
    for (int u = 0; u < 120; ++u) {
        for (int v = u + 1; v < 120; ++v) {
            if (u / 30 == v / 30 && (u * 7 + v * 3) % 4 != 0) { CHECK(graph_add_edge(c, u, v, 1.0)); }
        }
    }
    for (int k = 0; k < 3; ++k) { CHECK(graph_add_edge(c, k * 30, k * 30 + 30, 1.0)); }
    double cx[120], cy[120];
    CHECK(graph_layout(c, cx, cy, NULL) == GRAPH_OK);
    double within = 0, between = 0;
    size_t wn = 0, bn = 0;
    for (int u = 0; u < 120; ++u) {
        for (int v = u + 1; v < 120; ++v) {
            const double d = hypot(cx[u] - cx[v], cy[u] - cy[v]);
            if (u / 30 == v / 30) { within += d; ++wn; } else { between += d; ++bn; }
        }
    }
    CHECK(within / (double)wn < 0.5 * between / (double)bn);
    graph_destroy(c);

    /* Tiny cases and bad arguments. */
    Graph *e = graph_create(false);
    double ex[1], ey[1];
    CHECK(graph_layout(e, ex, ey, NULL) == GRAPH_OK);
    CHECK(graph_add_vertex(e) == 0 && graph_layout(e, ex, ey, NULL) == GRAPH_OK && ex[0] == 0.5);
    graph_destroy(e);
    CHECK(graph_layout(NULL, x, y, NULL) == GRAPH_ERR_INVALID_ARGUMENT);
    CHECK(graph_layout(g, NULL, y, NULL) == GRAPH_ERR_INVALID_ARGUMENT);
    graph_destroy(g);
}

/* ------------------------------------------------------------------ */
/*  Raster output, big graphs, options, errors                          */
/* ------------------------------------------------------------------ */

static void test_raster_and_scale(void)
{
    /* Same graph, either representation: the picture is identical byte for byte. */
    Graph *a = graph_gen_barabasi_albert(GRAPH_LIST, 400, 2, 5);
    Graph *b = graph_copy(a);
    CHECK(a != NULL && b != NULL && graph_convert(b, GRAPH_MATRIX));
    VizOptions o = graph_viz_default_options();
    o.format = VIZ_FORMAT_PNG;
    o.width = o.height = 300;
    FILE *fa = tmpfile(), *fb = tmpfile();
    CHECK(graph_visualize_to(a, fa, &o) == GRAPH_OK && graph_visualize_to(b, fb, &o) == GRAPH_OK);
    size_t la, lb;
    unsigned char *da = slurp(fa, &la), *db = slurp(fb, &lb);
    CHECK(la == lb && memcmp(da, db, la) == 0);
    unsigned w, h;
    unsigned char *rgb = decode_png(da, la, &w, &h);
    CHECK(w == 300 && h == 300);
    size_t non_white = 0;
    for (size_t i = 0; i < (size_t)w * h * 3; i += 3) {
        non_white += !(rgb[i] == 255 && rgb[i + 1] == 255 && rgb[i + 2] == 255);
    }
    CHECK(non_white > 2000);                          /* something was drawn */
    free(rgb); free(da); free(db); fclose(fa); fclose(fb);

    /* Frozen handles can be drawn. */
    graph_freeze(b);
    FILE *f = tmpfile();
    CHECK(graph_visualize_to(b, f, &o) == GRAPH_OK);
    fclose(f);
    graph_destroy(a); graph_destroy(b);

    /* A big graph (never SVG under AUTO) renders to a valid PNG. */
    Graph *big = graph_gen_barabasi_albert(GRAPH_LIST, 12000, 2, 9);
    CHECK(big != NULL);
    o = graph_viz_default_options();
    o.width = o.height = 500;
    f = tmpfile();
    CHECK(graph_visualize_to(big, f, &o) == GRAPH_OK);
    size_t len;
    unsigned char *file = slurp(f, &len);
    CHECK(len > 8 && memcmp(file + 1, "PNG", 3) == 0);   /* AUTO chose PNG */
    rgb = decode_png(file, len, &w, &h);
    CHECK(w == 500 && h == 500);
    free(rgb); free(file); fclose(f);

    /* Matrix heat-map: PNG of the requested size, blocks on the diagonal for clusters. */
    Graph *c = graph_create_sized(100, false);
    for (int u = 0; u < 100; ++u) {
        for (int v = u + 1; v < 100; ++v) {
            if (u / 25 == v / 25) { CHECK(graph_add_edge(c, u, v, 1.0)); }
        }
    }
    o = graph_viz_default_options();
    o.layout = VIZ_LAYOUT_MATRIX;
    o.width = o.height = 200;
    f = tmpfile();
    CHECK(graph_visualize_to(c, f, &o) == GRAPH_OK);
    file = slurp(f, &len);
    rgb = decode_png(file, len, &w, &h);
    CHECK(w == 200 && h == 200);
    #define PIXEL_DARK(px, py) (rgb[((size_t)(py) * 200 + (size_t)(px)) * 3] < 128)
    CHECK(PIXEL_DARK(10, 20) && PIXEL_DARK(60, 70));        /* inside block 0 and block 1 */
    CHECK(!PIXEL_DARK(10, 120) && !PIXEL_DARK(180, 20));    /* off-diagonal: empty */
    free(rgb); free(file); fclose(f);
    graph_destroy(c);

    /* More vertices than pixels: each pixel covers a bucket of the matrix (the density path). */
    Graph *dense = graph_create_sized(600, false);
    for (int u = 0; u < 600; ++u) {
        for (int v = u + 1; v < 600; ++v) {
            if (u / 150 == v / 150 && (u + v) % 2 == 0) { CHECK(graph_add_edge(dense, u, v, 1.0)); }
        }
    }
    f = tmpfile();
    CHECK(graph_visualize_to(dense, f, &o) == GRAPH_OK);
    file = slurp(f, &len);
    rgb = decode_png(file, len, &w, &h);
    CHECK(w == 200 && h == 200);
    size_t in_block = 0, off_block = 0;
    for (unsigned py = 0; py < 200; ++py) {
        for (unsigned px = 0; px < 200; ++px) {
            const bool dark = rgb[((size_t)py * 200 + px) * 3] < 200;
            if (px / 50 == py / 50) { in_block += dark; } else { off_block += dark; }
        }
    }
    CHECK(in_block > 3000 && off_block == 0);       /* strictly block-diagonal, rows not mixed up */
    free(rgb); free(file); fclose(f);
    graph_destroy(dense);
    graph_destroy(big);
}

/* ------------------------------------------------------------------ */
/*  Legend, community tags, temperature scale, hub names                 */
/* ------------------------------------------------------------------ */

static void test_legend_and_labels(void)
{
    /* Four dense groups of 40 joined by single bridges: Louvain must find them. */
    Graph *g = graph_create_sized(160, false);
    char *names[160];
    const char *const_names[160];
    for (int v = 0; v < 160; ++v) {
        names[v] = malloc(16);
        snprintf(names[v], 16, "p<%d>", v);                 /* the < and > must be escaped in SVG */
        const_names[v] = names[v];
    }
    for (int u = 0; u < 160; ++u) {
        for (int v = u + 1; v < 160; ++v) {
            if (u / 40 == v / 40 && (u * 13 + v * 7) % 3 != 0) { CHECK(graph_add_edge(g, u, v, 1.0)); }
        }
    }
    for (int k = 0; k < 3; ++k) { CHECK(graph_add_edge(g, k * 40, k * 40 + 40, 1.0)); }
    for (int v = 1; v < 30; ++v) { (void)graph_add_edge(g, 7, v + 100, 1.0); }      /* vertex 7: the hub */

    VizOptions o = graph_viz_default_options();
    o.format = VIZ_FORMAT_SVG;
    o.color = VIZ_COLOR_COMMUNITY;
    o.group_name = "Community";
    o.labels = const_names;
    o.show_labels = 0;
    o.label_top = 3;
    FILE *f = tmpfile();
    CHECK(graph_visualize_to(g, f, &o) == GRAPH_OK);
    size_t len;
    char *svg = (char *)slurp(f, &len);
    fclose(f);
    CHECK(strstr(svg, "Communities") != NULL);                          /* legend heading */
    CHECK(strstr(svg, "Community 1: 4") != NULL || strstr(svg, "Community 1:") != NULL);
    CHECK(count_substr(svg, "Community 1") >= 2 && count_substr(svg, "Community 4") >= 2);   /* tag + legend row */
    CHECK(strstr(svg, "hub: ") != NULL);                                /* each group names its hub */
    CHECK(strstr(svg, "p&lt;7&gt;") != NULL);                           /* the hub is labelled, escaped */
    CHECK(count_substr(svg, "<rect") >= 5);                             /* background + 4 swatches */
    free(svg);

    /* Without the legend: no panel, no tags, no heading. */
    o.legend = false;
    f = tmpfile();
    CHECK(graph_visualize_to(g, f, &o) == GRAPH_OK);
    svg = (char *)slurp(f, &len);
    fclose(f);
    CHECK(strstr(svg, "Communities") == NULL && strstr(svg, "Community 1") == NULL);
    free(svg);

    /* Temperature mode: colour bar, tick labels, popularity list, and the hub is the first entry. */
    o = graph_viz_default_options();
    o.format = VIZ_FORMAT_SVG;
    o.color = VIZ_COLOR_DEGREE;
    o.labels = const_names;
    o.show_labels = 0;
    o.label_top = 3;
    f = tmpfile();
    CHECK(graph_visualize_to(g, f, &o) == GRAPH_OK);
    svg = (char *)slurp(f, &len);
    fclose(f);
    CHECK(strstr(svg, "linearGradient") != NULL && strstr(svg, "Popularity") != NULL);
    CHECK(strstr(svg, "friends (log scale)") != NULL && strstr(svg, "Most popular") != NULL);
    CHECK(strstr(svg, "1. p&lt;7&gt;") != NULL);                        /* the hub tops the list */
    CHECK(count_substr(svg, "hot = very popular") == 1);
    /* The hub's circle is the hottest colour in the picture; ordinary people are cool. */
    CHECK(strstr(svg, "fill=\"#a50026\"") != NULL || strstr(svg, "fill=\"#a6") != NULL || strstr(svg, "fill=\"#a5") != NULL);
    free(svg);

    /* PNG: legend panel and text are drawn; with the legend off the right-hand panel stays empty. */
    o = graph_viz_default_options();
    o.format = VIZ_FORMAT_PNG;
    o.color = VIZ_COLOR_COMMUNITY;
    o.labels = const_names;
    o.width = 1000;
    o.height = 800;
    o.title = "Test";
    unsigned w, h;
    size_t ink_on = 0, ink_off = 0;
    for (int pass = 0; pass < 2; ++pass) {
        o.legend = pass == 0;
        f = tmpfile();
        CHECK(graph_visualize_to(g, f, &o) == GRAPH_OK);
        unsigned char *file = slurp(f, &len);
        unsigned char *rgb = decode_png(file, len, &w, &h);
        CHECK(w == 1000 && h == 800);
        size_t ink = 0;
        for (unsigned y = 0; y < h; ++y) {
            for (unsigned x = 800; x < w; ++x) {              /* the right-hand panel */
                const unsigned char *p = rgb + ((size_t)y * w + x) * 3;
                ink += !(p[0] > 250 && p[1] > 250 && p[2] > 250);
            }
        }
        (pass == 0 ? &ink_on : &ink_off)[0] = ink;
        free(rgb); free(file); fclose(f);
    }
    CHECK(ink_on > 1500);                    /* swatches and text */
    CHECK(ink_off < ink_on / 3);             /* off: the area carries (at most) a few plot pixels */

    for (int v = 0; v < 160; ++v) { free(names[v]); }
    graph_destroy(g);
}

static void test_files_and_errors(void)
{
    Graph *g = graph_gen_grid(GRAPH_LIST, 5, 5);
    CHECK(g != NULL);
    const char *svg_path = "viz_test_out.svg", *png_path = "viz_test_out.png";
    CHECK(graph_visualize(g, svg_path, NULL) == GRAPH_OK);           /* extension decides */
    CHECK(graph_visualize(g, png_path, NULL) == GRAPH_OK);
    FILE *f = fopen(svg_path, "rb");
    char head[8] = { 0 };
    CHECK(f != NULL && fread(head, 1, 5, f) == 5 && strcmp(head, "<svg ") == 0);
    fclose(f);
    f = fopen(png_path, "rb");
    unsigned char sig[4] = { 0 };
    CHECK(f != NULL && fread(sig, 1, 4, f) == 4 && sig[1] == 'P');
    fclose(f);
    remove(svg_path);
    remove(png_path);

    graph_clear_error();
    CHECK(graph_visualize(g, "/nonexistent_dir/x.png", NULL) == GRAPH_ERR_IO);
    CHECK(graph_last_error() == GRAPH_ERR_IO);
    CHECK(graph_visualize(NULL, "x.png", NULL) == GRAPH_ERR_INVALID_ARGUMENT);
    CHECK(graph_visualize(g, NULL, NULL) == GRAPH_ERR_INVALID_ARGUMENT);
    CHECK(graph_visualize_to(g, NULL, NULL) == GRAPH_ERR_INVALID_ARGUMENT);

    VizOptions o = graph_viz_default_options();
    o.width = 10;
    CHECK(graph_visualize(g, "viz_test_bad.png", &o) == GRAPH_ERR_INVALID_ARGUMENT);
    CHECK(fopen("viz_test_bad.png", "rb") == NULL);                   /* no half-written file left */
    o = graph_viz_default_options();
    o.width = o.height = 8192;                                         /* 67M pixels: refused, not crashed */
    o.format = VIZ_FORMAT_PNG;
    FILE *t = tmpfile();
    CHECK(graph_visualize_to(g, t, &o) == GRAPH_ERR_TOO_LARGE);
    fclose(t);
    graph_destroy(g);
}

int main(void)
{
    test_png_encoder();
    test_svg();
    test_layout();
    test_raster_and_scale();
    test_legend_and_labels();
    test_files_and_errors();
    puts("All visualization tests passed.");
    return 0;
}
