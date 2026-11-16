/*
 * million_users.c — a social network of ONE MILLION users, analysed and drawn.
 *
 * It builds a made-up but realistic network (communities, a few very popular hubs, and the 20% of users
 * who never connect to anyone), runs the library's analyses on it, prints what it finds, and draws a
 * picture for each: the degree distribution, the friendships themselves, connected groups, communities,
 * PageRank, k-core, clustering, how far information travels, what the memory costs, and what the
 * sharding engine does on a small device.
 *
 *   make build/examples/million_users
 *   ./build/examples/million_users docs/images/million        # about 4-5 minutes, < 2 GB
 *   ./build/examples/million_users docs/images/million 100000 # the same for 100 000 users (about 20 s)
 *
 * It is also the library's BENCHMARK.  Every analysis is timed (best / median / worst of several runs),
 * its answer is checked by independent code in this file (never by the library itself), the peak memory is
 * recorded, and a report is written next to the pictures: REPORT.md (also printed) and results.csv.
 * Build it with optimisation or the numbers mean nothing:   make bench-million
 * Compare with an earlier run:   make bench-million BASELINE=path/to/old/results.csv
 *
 * The network drawings come from the library's visualiser (graph_viz.h).  The charts (histograms, bar
 * charts, scatter plots) are not part of the library: this file draws them itself with a few helpers
 * at the top, on top of the library's internal PNG encoder and bitmap font (hence -Isrc).
 */

#define _POSIX_C_SOURCE 200809L

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/utsname.h>
#include <time.h>
#include <unistd.h>

#include "graph_algorithms.h"
#include "graph_shard.h"
#include "graph_social.h"
#include "graph_viz.h"
#include "graph_internal.h"      /* the backend's memory_bytes(), to report what a graph costs */
#include "png_writer.h"
#include "viz_font.h"

/* ================================================================== */
/*  Tiny chart toolkit (example-only)                                   */
/* ================================================================== */

typedef struct { int w, h; unsigned char *px; } Img;
typedef struct { unsigned char r, g, b; } Rgb;

static const Rgb INK = {30, 34, 48}, GREY = {110, 116, 130}, GRID = {228, 230, 236};
static const Rgb BLUE = {45, 105, 205}, ORANGE = {232, 122, 30};

static Img *img_new(int w, int h)
{
    Img *im = malloc(sizeof *im);
    im->w = w; im->h = h;
    im->px = malloc((size_t)w * (size_t)h * 3);
    memset(im->px, 255, (size_t)w * (size_t)h * 3);
    return im;
}

static void put(Img *im, int x, int y, Rgb c)
{
    if (x < 0 || y < 0 || x >= im->w || y >= im->h) return;
    unsigned char *p = im->px + ((size_t)y * (size_t)im->w + (size_t)x) * 3;
    p[0] = c.r; p[1] = c.g; p[2] = c.b;
}

static void rect(Img *im, int x0, int y0, int x1, int y1, Rgb c)      /* inclusive corners */
{
    if (x0 > x1) { int t = x0; x0 = x1; x1 = t; }
    if (y0 > y1) { int t = y0; y0 = y1; y1 = t; }
    for (int y = y0; y <= y1; y++) for (int x = x0; x <= x1; x++) put(im, x, y, c);
}

static void line(Img *im, int x0, int y0, int x1, int y1, Rgb c)
{
    const int dx = abs(x1 - x0), dy = -abs(y1 - y0), sx = x0 < x1 ? 1 : -1, sy = y0 < y1 ? 1 : -1;
    int err = dx + dy;
    for (;;) {
        put(im, x0, y0, c);
        if (x0 == x1 && y0 == y1) break;
        const int e2 = 2 * err;
        if (e2 >= dy) { err += dy; x0 += sx; }
        if (e2 <= dx) { err += dx; y0 += sy; }
    }
}

static int text_w(const char *s, int scale) { return (int)strlen(s) * VIZ_FONT_W * scale; }

static void text(Img *im, int x, int y, const char *s, int scale, Rgb c)
{
    for (; *s; s++, x += VIZ_FONT_W * scale) {
        if (*s < 32 || *s > 126) continue;
        const unsigned char *g = VIZ_FONT[*s - 32];
        for (int r = 0; r < VIZ_FONT_H; r++)
            for (int b = 0; b < VIZ_FONT_W; b++)
                if (g[r] & (0x80 >> b)) rect(im, x + b * scale, y + r * scale, x + b * scale + scale - 1, y + r * scale + scale - 1, c);
    }
}

/* Text that wraps onto a second line (at a space) if it would run past max_w pixels. */
static void text_wrapped(Img *im, int x, int y, const char *s, int scale, Rgb c, int max_w)
{
    if (text_w(s, scale) <= max_w) { text(im, x, y, s, scale, c); return; }
    char buf[420];
    snprintf(buf, sizeof buf, "%s", s);
    int cut = max_w / (VIZ_FONT_W * scale);
    if ((size_t)cut >= strlen(buf)) cut = (int)strlen(buf) - 1;
    while (cut > 0 && buf[cut] != ' ') cut--;
    if (cut == 0) cut = max_w / (VIZ_FONT_W * scale);
    char first[420];
    snprintf(first, sizeof first, "%.*s", cut, buf);
    text(im, x, y, first, scale, c);
    char second[420];
    snprintf(second, sizeof second, "%s", buf + cut + (buf[cut] == ' ' ? 1 : 0));
    if (text_w(second, scale) > max_w) { second[max_w / (VIZ_FONT_W * scale) - 3] = 0; strcat(second, "..."); }
    text(im, x, y + VIZ_FONT_H * scale + 2, second, scale, c);
}

static void save(Img *im, const char *dir, const char *name)
{
    char path[600];
    snprintf(path, sizeof path, "%s/%s", dir, name);
    FILE *f = fopen(path, "wb");
    if (f == NULL || png_write_rgb(f, im->px, (unsigned)im->w, (unsigned)im->h) != GRAPH_OK) {
        fprintf(stderr, "cannot write %s\n", path);
        exit(1);
    }
    fclose(f);
    printf("  drew %s\n", name);
    free(im->px); free(im);
}

/* A plot area with a title, a subtitle and axes.  Coordinates are mapped by map_x / map_y. */
typedef struct {
    Img *im;
    int x0, y0, x1, y1;                       /* the plot rectangle in pixels */
    double xmin, xmax, ymin, ymax;
    bool logx, logy;
} Plot;

static double tf(double v, bool lg) { return lg ? log10(v < 1e-300 ? 1e-300 : v) : v; }
static int map_x(const Plot *p, double v) { return p->x0 + (int)((tf(v, p->logx) - tf(p->xmin, p->logx)) / (tf(p->xmax, p->logx) - tf(p->xmin, p->logx)) * (p->x1 - p->x0)); }
static int map_y(const Plot *p, double v) { return p->y1 - (int)((tf(v, p->logy) - tf(p->ymin, p->logy)) / (tf(p->ymax, p->logy) - tf(p->ymin, p->logy)) * (p->y1 - p->y0)); }

static void fmt_num(char *out, size_t cap, double v)
{
    if (v >= 1e6 && fabs(v - round(v / 1e6) * 1e6) < 1) snprintf(out, cap, "%.0fM", v / 1e6);
    else if (v >= 1e3 && fabs(v - round(v / 1e3) * 1e3) < 1) snprintf(out, cap, "%.0fk", v / 1e3);
    else if (v >= 100 || v == floor(v)) snprintf(out, cap, "%.0f", v);
    else snprintf(out, cap, "%.3g", v);
}

static Plot plot_begin(int w, int h, const char *title, const char *subtitle, double xmin, double xmax, double ymin,
                       double ymax, bool logx, bool logy, const char *xlabel, const char *ylabel, bool xticks)
{
    Plot p = {img_new(w, h), 110, 124, w - 40, h - 80, xmin, xmax, ymin, ymax, logx, logy};
    text(p.im, 30, 22, title, 2, INK);
    text_wrapped(p.im, 30, 58, subtitle, 1, GREY, w - 60);
    /* grid and tick labels */
    char buf[32];
    for (int axis = xticks ? 0 : 1; axis < 2; axis++) {      /* a category chart skips the numeric x axis */
        const bool lg = axis == 0 ? logx : logy;
        const double lo = axis == 0 ? xmin : ymin, hi = axis == 0 ? xmax : ymax;
        if (lg) {
            for (double v = pow(10, floor(log10(lo))); v <= hi * 1.0001; v *= 10) {
                if (v < lo * 0.9999) continue;
                fmt_num(buf, sizeof buf, v);
                if (axis == 0) { const int x = map_x(&p, v); line(p.im, x, p.y0, x, p.y1, GRID); text(p.im, x - text_w(buf, 1) / 2, p.y1 + 8, buf, 1, GREY); }
                else { const int y = map_y(&p, v); line(p.im, p.x0, y, p.x1, y, GRID); text(p.im, p.x0 - 8 - text_w(buf, 1), y - 8, buf, 1, GREY); }
            }
        } else {
            const double step = pow(10, floor(log10((hi - lo) / 5.0))) * ((hi - lo) / 5.0 / pow(10, floor(log10((hi - lo) / 5.0))) > 5 ? 10 : (hi - lo) / 5.0 / pow(10, floor(log10((hi - lo) / 5.0))) > 2 ? 5 : (hi - lo) / 5.0 / pow(10, floor(log10((hi - lo) / 5.0))) > 1 ? 2 : 1);
            for (double v = ceil(lo / step) * step; v <= hi + step * 1e-6; v += step) {
                fmt_num(buf, sizeof buf, v);
                if (axis == 0) { const int x = map_x(&p, v); line(p.im, x, p.y0, x, p.y1, GRID); text(p.im, x - text_w(buf, 1) / 2, p.y1 + 8, buf, 1, GREY); }
                else { const int y = map_y(&p, v); line(p.im, p.x0, y, p.x1, y, GRID); text(p.im, p.x0 - 8 - text_w(buf, 1), y - 8, buf, 1, GREY); }
            }
        }
    }
    line(p.im, p.x0, p.y1, p.x1, p.y1, GREY);
    line(p.im, p.x0, p.y0, p.x0, p.y1, GREY);
    text(p.im, (p.x0 + p.x1) / 2 - text_w(xlabel, 1) / 2, h - 36, xlabel, 1, INK);
    text(p.im, 30, 100, ylabel, 1, INK);
    return p;
}

/* Vertical bars; labels (optional) under each bar. */
static void chart_bars(const char *dir, const char *file, const char *title, const char *subtitle, const char *xlabel,
                       const char *ylabel, const char *const *labels, const double *v, size_t n, bool logy, int highlight)
{
    double mx = 1;
    for (size_t i = 0; i < n; i++) if (v[i] > mx) mx = v[i];
    Plot p = plot_begin(1100, 640, title, subtitle, 0, (double)n, logy ? 1 : 0, logy ? mx * 2 : mx * 1.1, false, logy, xlabel, ylabel, false);
    const double bw = (double)(p.x1 - p.x0) / (double)n;
    for (size_t i = 0; i < n; i++) {
        const int xa = p.x0 + (int)(bw * (double)i) + 2, xb = p.x0 + (int)(bw * (double)(i + 1)) - 2;
        const double val = v[i] < (logy ? 1 : 0) ? (logy ? 1 : 0) : v[i];
        rect(p.im, xa, map_y(&p, val), xb, p.y1 - 1, (int)i == highlight ? ORANGE : BLUE);
        char buf[32];
        fmt_num(buf, sizeof buf, v[i]);
        if (n <= 30 && v[i] > 0) text(p.im, (xa + xb) / 2 - text_w(buf, 1) / 2, map_y(&p, val) - 18, buf, 1, INK);
        if (labels != NULL && labels[i] != NULL && (n <= 30 || i % (n / 20 + 1) == 0))
            text(p.im, (xa + xb) / 2 - text_w(labels[i], 1) / 2, p.y1 + 28, labels[i], 1, GREY);
    }
    save(p.im, dir, file);
}

/* Horizontal bars with a label to the left of each. */
static void chart_hbars(const char *dir, const char *file, const char *title, const char *subtitle, const char *unit,
                        const char *const *labels, const double *v, size_t n, bool logx, int highlight)
{
    double mx = 1e-300;
    for (size_t i = 0; i < n; i++) if (v[i] > mx) mx = v[i];
    int labw = 0;
    for (size_t i = 0; i < n; i++) if (text_w(labels[i], 1) > labw) labw = text_w(labels[i], 1);
    const int w = 1100, h = 170 + (int)n * 44;
    Img *im = img_new(w, h);
    text(im, 30, 22, title, 2, INK);
    text_wrapped(im, 30, 58, subtitle, 1, GREY, w - 60);
    const int x0 = 30 + labw + 20, x1 = w - 190, y0 = 124;
    const double lo = logx ? 1e-3 : 0, hi = logx ? pow(10, ceil(log10(mx)) + 0.0) : mx * 1.05;
    for (size_t i = 0; i < n; i++) {
        const int y = y0 + (int)i * 44;
        text(im, x0 - 12 - text_w(labels[i], 1), y + 8, labels[i], 1, INK);
        const double frac = logx ? (log10(v[i] < lo ? lo : v[i]) - log10(lo)) / (log10(hi) - log10(lo)) : v[i] / hi;
        rect(im, x0, y, x0 + (int)(frac * (x1 - x0)), y + 28, (int)i == highlight ? ORANGE : BLUE);
        char buf[48], num[32];
        fmt_num(num, sizeof num, v[i]);
        snprintf(buf, sizeof buf, "%s %s", num, unit);
        text(im, x0 + (int)(frac * (x1 - x0)) + 10, y + 8, buf, 1, INK);
    }
    line(im, x0, y0 - 6, x0, y0 + (int)n * 44, GREY);
    save(im, dir, file);
}

/* Points on a plot (optionally joined by lines). */
static void chart_points(const char *dir, const char *file, const char *title, const char *subtitle, const char *xlabel,
                         const char *ylabel, const double *xs, const double *ys, size_t n, bool logx, bool logy, bool join)
{
    double xmin = 1e300, xmax = -1e300, ymin = 1e300, ymax = -1e300;
    for (size_t i = 0; i < n; i++) {
        if (logx && xs[i] <= 0) continue;
        if (logy && ys[i] <= 0) continue;
        if (xs[i] < xmin) { xmin = xs[i]; }
        if (xs[i] > xmax) { xmax = xs[i]; }
        if (ys[i] < ymin) { ymin = ys[i]; }
        if (ys[i] > ymax) { ymax = ys[i]; }
    }
    if (logx) { xmin = pow(10, floor(log10(xmin))); xmax = pow(10, ceil(log10(xmax))); } else { xmin = 0; xmax *= 1.05; }
    if (logy) { ymin = pow(10, floor(log10(ymin))); ymax = pow(10, ceil(log10(ymax))); } else { ymin = 0; ymax *= 1.1; }
    Plot p = plot_begin(1100, 700, title, subtitle, xmin, xmax, ymin, ymax, logx, logy, xlabel, ylabel, true);
    int px = -1, py = -1;
    for (size_t i = 0; i < n; i++) {
        if ((logx && xs[i] <= 0) || (logy && ys[i] <= 0)) continue;
        const int x = map_x(&p, xs[i]), y = map_y(&p, ys[i]);
        if (join && px >= 0) line(p.im, px, py, x, y, BLUE);
        rect(p.im, x - 2, y - 2, x + 2, y + 2, join ? BLUE : (Rgb){45, 105, 205});
        px = x; py = y;
    }
    save(p.im, dir, file);
}

/* ================================================================== */
/*  The network                                                         */
/* ================================================================== */

static uint64_t rs = 88172645463325252ULL;
static uint64_t rnd64(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static unsigned rnd(unsigned bound) { return (unsigned)(rnd64() % bound); }
static double u01(void) { return (double)(rnd64() >> 11) / 9007199254740992.0; }

static double now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec + (double)t.tv_nsec * 1e-9;
}

typedef struct { int *a; size_t n, cap; } Vec;

static void push(Vec *v, int x)
{
    if (v->n == v->cap) {
        v->cap = v->cap ? v->cap * 2 : 16;
        v->a = realloc(v->a, v->cap * sizeof *v->a);
        if (v->a == NULL) { fprintf(stderr, "out of memory\n"); exit(1); }
    }
    v->a[v->n++] = x;
}

/*
 * The model.  Users belong to communities (a few large ones, many small); each active user makes
 * 1-80 friendships when they join.  85% of the time the friend is in their own community, and most
 * of those picks favour people who already have many friends (so each community grows hubs); 15% of
 * the time the friend comes from anywhere, again favouring the popular.  Every fifth user is silent:
 * they exist but never befriend anyone and nobody befriends them.  Result: a heavy-tailed degree
 * distribution (hubs, but not a clean power law), clear communities, and a giant connected group.
 */
static Graph *make_network(size_t users, size_t ncomm, int *comm, Vec *eu, Vec *ev)
{
    GraphConfig cfg = graph_default_config(false);
    cfg.unweighted = true;                      /* friendships have no weight: 23 bytes each, not 47 */
    cfg.adaptive   = false;
    cfg.cache_entries = 0;                      /* no result cache: the core calls (components, Dijkstra, ...) would
                                                   otherwise answer repeat runs from memory, and a benchmark must
                                                   time the computation, not the cache */
    Graph *g = graph_create_with(&cfg, users);
    if (g == NULL) { fprintf(stderr, "cannot create the graph\n"); exit(1); }
    Vec *ends = calloc(ncomm, sizeof *ends), *members = calloc(ncomm, sizeof *members);
    Vec gends = {0}, active = {0};
    for (size_t u = 0; u < users; u++) comm[u] = (int)((double)ncomm * pow(u01(), 1.7));
    for (size_t v = 0; v < users; v++) {
        if (v % 5 == 4) continue;                                       /* the silent fifth */
        const int c = comm[v];
        const int links = 1 + (int)fmin(79.0, -log(1.0 - u01()) * 3.5);
        for (int k = 0; k < links; k++) {
            int f = -1;
            if (rnd(100) < 85 && members[c].n > 0)
                f = (rnd(100) < 70 && ends[c].n > 0) ? ends[c].a[rnd((unsigned)ends[c].n)] : members[c].a[rnd((unsigned)members[c].n)];
            else if (active.n > 0)
                f = (rnd(100) < 70 && gends.n > 0) ? gends.a[rnd((unsigned)gends.n)] : active.a[rnd((unsigned)active.n)];
            if (f < 0 || f == (int)v || !graph_add_edge(g, (int)v, f, 1.0)) continue;
            push(eu, (int)v); push(ev, f);
            push(&ends[c], (int)v); push(&ends[comm[f]], f);
            push(&gends, (int)v); push(&gends, f);
        }
        push(&members[c], (int)v);
        push(&active, (int)v);
    }
    for (size_t c = 0; c < ncomm; c++) { free(ends[c].a); free(members[c].a); }
    free(ends); free(members); free(gends.a); free(active.a);
    return g;
}

/* ================================================================== */
/*  Benchmark harness: stages, repeats, memory, independent checks      */
/* ================================================================== */

#ifndef BUILD_CFLAGS
#define BUILD_CFLAGS "(unknown: compile with the Makefile)"
#endif

#define MAXREP 9
#define MAXSTAGE 40

typedef struct {
    char   name[64];
    int    reps;
    double min, med, max;
    double rate;            /* million friendships per second at the BEST run; 0 = not applicable */
    bool   scored;          /* counted in the headline score */
    long   rss_mb;          /* peak resident memory of the whole process so far */
    bool   ok;              /* the independent check passed */
    char   check[200];      /* what was checked */
    double baseline;        /* median of the same stage in the baseline run, or 0 */
} Stage;

static Stage g_stage[MAXSTAGE];
static int   g_nstage;
static int   g_repeats = 5;

static long rss_peak_mb(void)               /* the process's high-water mark of resident memory */
{
    FILE *f = fopen("/proc/self/status", "r");
    char line[256];
    long kb = 0;
    while (f != NULL && fgets(line, sizeof line, f)) if (strncmp(line, "VmHWM:", 6) == 0) sscanf(line + 6, "%ld", &kb);
    if (f != NULL) fclose(f);
    return kb / 1024;
}

static int cmp_d(const void *a, const void *b)
{
    const double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

/* Records a finished stage: `ts` are the wall times of `reps` runs. */
static Stage *stage_done(const char *name, double *ts, int reps, size_t work_edges, bool scored, bool ok, const char *check)
{
    Stage *s = &g_stage[g_nstage++];
    snprintf(s->name, sizeof s->name, "%s", name);
    qsort(ts, (size_t)reps, sizeof *ts, cmp_d);
    s->reps = reps; s->min = ts[0]; s->max = ts[reps - 1]; s->med = ts[reps / 2];
    s->rate = work_edges > 0 && s->min > 0 ? (double)work_edges / s->min / 1e6 : 0;      /* best run: least disturbed */
    s->scored = scored; s->ok = ok; s->rss_mb = rss_peak_mb();
    snprintf(s->check, sizeof s->check, "%s", check);
    printf("  %-34s %8.3f s (best %.3f, worst %.3f)  %s%s\n", name, s->med, s->min, s->max, ok ? "check passed: " : "CHECK FAILED: ", check);
    fflush(stdout);
    return s;
}

#define REPEAT(ts, reps, work) do { for (int r_ = 0; r_ < (reps); r_++) { const double t0_ = now(); work; (ts)[r_] = now() - t0_; } } while (0)

/* ---- independent code: the answers are checked against THESE, built only from the generator's own edge list ---- */

static int uf_find(int *p, int x) { while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; } return x; }

static bool verify_components(size_t users, const Vec *eu, const Vec *ev, const int *comp, size_t ngroups, char *msg, size_t cap)
{
    int *p = malloc(users * sizeof *p), *to = malloc(users * sizeof *to), *back = malloc(users * sizeof *back);
    for (size_t v = 0; v < users; v++) { p[v] = (int)v; to[v] = back[v] = -1; }
    size_t groups = users;
    for (size_t i = 0; i < eu->n; i++) {
        const int a = uf_find(p, eu->a[i]), b = uf_find(p, ev->a[i]);
        if (a != b) { p[a] = b; groups--; }
    }
    bool ok = groups == ngroups;
    for (size_t v = 0; v < users && ok; v++) {              /* the labels must be the same partition */
        const int r = uf_find(p, (int)v);
        if (to[comp[v]] == -1 && back[r] == -1) { to[comp[v]] = r; back[r] = comp[v]; }
        else if (to[comp[v]] != r || back[r] != comp[v]) ok = false;
    }
    snprintf(msg, cap, "union-find over the %zu edges gives %zu groups and the same partition", eu->n, groups);
    free(p); free(to); free(back);
    return ok;
}

static bool verify_core(Graph *g, size_t users, const int *core, char *msg, size_t cap)
{
    bool ok = true;
    size_t maxd = 0;
    for (size_t v = 0; v < users; v++) if (graph_degree(g, (int)v) > maxd) maxd = graph_degree(g, (int)v);
    int *nb = malloc((maxd + 1) * sizeof *nb);
    for (size_t v = 0; v < users && ok; v++) {              /* a k-core member keeps >= k friends inside the k-core */
        const size_t d = graph_neighbors(g, (int)v, nb, NULL, maxd);
        size_t inside = 0;
        for (size_t i = 0; i < d; i++) inside += core[nb[i]] >= core[v];
        if (inside < (size_t)core[v] || (size_t)core[v] > d) ok = false;
    }
    snprintf(msg, cap, "every user keeps at least k friends inside their own k-core, and k <= degree");
    free(nb);
    return ok;
}

static bool verify_clustering(Graph *g, size_t users, const double *local, long long triangles, char *msg, size_t cap)
{
    bool ok = true;
    double sum_tri = 0;
    for (size_t v = 0; v < users; v++) {
        const double d = (double)graph_degree(g, (int)v);
        sum_tri += local[v] * d * (d - 1) / 2.0;            /* a user's triangles = clustering x pairs of friends */
    }
    ok = fabs(sum_tri / 3.0 - (double)triangles) <= 1e-6 * (double)triangles + 1.0;
    size_t checked = 0;
    int nb[128];
    for (size_t k = 0; k < 4000 && ok; k++) {               /* brute force on a sample of users */
        const size_t v = (size_t)rnd((unsigned)users);
        const size_t d = graph_degree(g, (int)v);
        if (d < 2 || d > 120) continue;
        graph_neighbors(g, (int)v, nb, NULL, 128);
        size_t links = 0;
        for (size_t a = 0; a < d; a++) for (size_t b = a + 1; b < d; b++) links += graph_has_edge(g, nb[a], nb[b], NULL);
        if (fabs(local[v] - 2.0 * (double)links / ((double)d * (double)(d - 1))) > 1e-9) ok = false;
        checked++;
    }
    snprintf(msg, cap, "triangles = sum of local clustering / 3, and %zu users recomputed by brute force", checked);
    return ok;
}

static bool verify_pagerank(size_t users, const Vec *eu, const Vec *ev, const double *rank, char *msg, size_t cap)
{
    double sum = 0, dangling = 0, diff = 0;
    size_t *deg = calloc(users, sizeof *deg);
    double *next = calloc(users, sizeof *next);
    for (size_t i = 0; i < eu->n; i++) { deg[eu->a[i]]++; deg[ev->a[i]]++; }
    for (size_t v = 0; v < users; v++) { sum += rank[v]; if (deg[v] == 0) dangling += rank[v]; }
    for (size_t i = 0; i < eu->n; i++) { next[ev->a[i]] += rank[eu->a[i]] / (double)deg[eu->a[i]]; next[eu->a[i]] += rank[ev->a[i]] / (double)deg[ev->a[i]]; }
    for (size_t v = 0; v < users; v++) diff += fabs(0.15 / (double)users + 0.85 * (next[v] + dangling / (double)users) - rank[v]);
    snprintf(msg, cap, "ranks sum to %.9f and one more power-iteration step moves them by only %.1e", sum, diff);
    free(deg); free(next);
    return fabs(sum - 1.0) < 1e-6 && diff < 1e-6;
}

static bool verify_louvain(size_t users, const Vec *eu, const Vec *ev, const int *lc, size_t ncl, double q_reported, char *msg, size_t cap)
{
    double *in = calloc(ncl, sizeof *in), *tot = calloc(ncl, sizeof *tot);
    const double m = (double)eu->n;
    for (size_t i = 0; i < eu->n; i++) {
        const int a = lc[eu->a[i]], b = lc[ev->a[i]];
        if (a == b) in[a] += 1;
        tot[a] += 1; tot[b] += 1;
    }
    double q = 0;
    for (size_t c = 0; c < ncl; c++) q += in[c] / m - (tot[c] / (2 * m)) * (tot[c] / (2 * m));
    (void)users;
    snprintf(msg, cap, "modularity recomputed from the edges: %.6f (library reported %.6f)", q, q_reported);
    free(in); free(tot);
    return fabs(q - q_reported) < 1e-6;
}

static bool verify_hops(const Vec *eu, const Vec *ev, const double *dist, const int *par, size_t hub, size_t users, char *msg, size_t cap)
{
    bool ok = dist[hub] == 0;
    for (size_t i = 0; i < eu->n && ok; i++) {              /* friends are at most one hop apart; both reachable or neither */
        const double a = dist[eu->a[i]], b = dist[ev->a[i]];
        if ((a > 1e300) != (b > 1e300) || (a < 1e300 && fabs(a - b) > 1)) ok = false;
    }
    for (size_t v = 0; v < users && ok; v++) if (v != hub && dist[v] < 1e300 && dist[v] != dist[par[v]] + 1) ok = false;
    snprintf(msg, cap, "adjacent users differ by <= 1 hop, and every user is exactly one hop beyond their parent");
    return ok;
}

static bool verify_suggestions(Graph *g, size_t users, char *msg, size_t cap)
{
    unsigned *cnt = calloc(users, sizeof *cnt);
    bool ok = true;
    size_t users_checked = 0;
    for (int trial = 0; trial < 400 && users_checked < 12; trial++) {
        const int u = (int)rnd((unsigned)users);
        const size_t d = graph_degree(g, u);
        if (d < 3 || d > 150) continue;
        int nb[160], nb2[4096];
        graph_neighbors(g, u, nb, NULL, 160);
        Vec touched = {0};
        for (size_t i = 0; i < d; i++) {
            const size_t d2 = graph_degree(g, nb[i]);
            if (d2 > 4096) continue;
            graph_neighbors(g, nb[i], nb2, NULL, 4096);
            for (size_t j = 0; j < d2; j++) {
                if (nb2[j] == u) continue;
                if (cnt[nb2[j]]++ == 0) push(&touched, nb2[j]);
            }
        }
        for (size_t i = 0; i < d; i++) { if (cnt[nb[i]] > 0) { cnt[nb[i]] = 0; } }   /* already friends: not candidates */
        GraphSuggestion sug[5];
        size_t count = 0;
        graph_social_suggest_friends(g, u, 5, sug, &count);
        /* brute-force best: highest count, lowest id */
        unsigned best = 0;
        int best_id = -1;
        for (size_t i = 0; i < touched.n; i++) {
            const int w = touched.a[i];
            if (cnt[w] > best || (cnt[w] == best && best > 0 && w < best_id)) { best = cnt[w]; best_id = w; }
        }
        if (best > 0 && (count == 0 || sug[0].user != best_id || sug[0].mutual != best)) ok = false;
        for (size_t i = 0; i < touched.n; i++) cnt[touched.a[i]] = 0;
        free(touched.a);
        users_checked++;
    }
    snprintf(msg, cap, "top suggestion matches a brute-force count of mutual friends for %zu users", users_checked);
    free(cnt);
    return ok;
}

/* ---- environment ---- */

typedef struct { char cpu[160], kernel[220], compiler[160]; long cores; double ram_gb; bool optimised; double load_start, load_end; } Env;

static double load_now(void)             /* 1-minute load average: how busy the machine is with other work */
{
    double l = 0;
    FILE *f = fopen("/proc/loadavg", "r");
    if (f != NULL) { if (fscanf(f, "%lf", &l) != 1) l = 0; fclose(f); }
    return l;
}

static Env capture_env(void)
{
    Env e = {"unknown", "unknown", "", sysconf(_SC_NPROCESSORS_ONLN), 0, false, 0, 0};
    FILE *f = fopen("/proc/cpuinfo", "r");
    char line[300];
    while (f != NULL && fgets(line, sizeof line, f))
        if (strncmp(line, "model name", 10) == 0) {
            const char *c = strchr(line, ':');
            if (c != NULL) { snprintf(e.cpu, sizeof e.cpu, "%s", c + 2); e.cpu[strcspn(e.cpu, "\n")] = 0; }
            break;
        }
    if (f != NULL) fclose(f);
    f = fopen("/proc/meminfo", "r");
    while (f != NULL && fgets(line, sizeof line, f)) if (strncmp(line, "MemTotal:", 9) == 0) { long kb; sscanf(line + 9, "%ld", &kb); e.ram_gb = (double)kb / 1048576.0; }
    if (f != NULL) fclose(f);
    e.load_start = load_now();
    struct utsname u;
    if (uname(&u) == 0) snprintf(e.kernel, sizeof e.kernel, "%s %s %s", u.sysname, u.release, u.machine);
#ifdef __VERSION__
    snprintf(e.compiler, sizeof e.compiler, "%s", __VERSION__);
#endif
#ifdef __OPTIMIZE__
    e.optimised = true;
#endif
    return e;
}

/* ================================================================== */
/*  The report                                                          */
/* ================================================================== */

typedef struct {
    size_t users, edges, silent, maxdeg, hub, ncomm_generated;
    double build_s, avg_deg;
    struct { char label[56]; double mb, bytes_per_edge; } mem[6];
    int nmem;
    struct { size_t shards, cache, loads; double cache_mb, cut_s, comp_s, pr_s, read_mb, diff; bool ok; } sh;
    struct { char name[56]; size_t n, ok; double p50, p99, max, qps; } q[3];
    int nq;
    struct { char file[40]; char caption[170]; double secs; } pic[24];
    int npic;
    long rss_analysis_mb, rss_final_mb;
    struct { struct { char name[40]; double first_ms, repeat_ms; bool same; } c[3];
             uint64_t hits, misses, stores, invalidations; size_t bytes; bool fresh_ok, ok; } ca;
} Rep;

static void emit_report(FILE *f, const Env *e, const Rep *r, bool has_baseline)
{
    double suite = 0, logsum = 0;
    int scored = 0, passed = 0;
    for (int i = 0; i < g_nstage; i++) {
        passed += g_stage[i].ok;
        if (g_stage[i].scored) { suite += g_stage[i].min; if (g_stage[i].rate > 0) { logsum += log(g_stage[i].rate); scored++; } }
    }
    const double score = scored > 0 ? exp(logsum / scored) : 0;
    fprintf(f, "# Benchmark report: %zu-user social network\n\n", r->users);
    fprintf(f, "Generated by `examples/million_users.c`. Every analysis is one library call on the whole network, timed over %d run(s);\n", g_repeats);
    fprintf(f, "every answer is checked by independent code in the example (not by the library), so a fast wrong answer does not score.\n\n");

    fprintf(f, "## Headline\n\n");
    fprintf(f, "| | |\n|---|---|\n");
    fprintf(f, "| **Suite score** (geometric mean of the analyses' best-run rates) | **%.1f million friendships per second** |\n", score);
    fprintf(f, "| Analysis suite time (sum of the %d scored best runs) | %.2f s |\n", scored, suite);
    fprintf(f, "| Independent checks passed | %d of %d |\n", passed, g_nstage);
    fprintf(f, "| Peak resident memory after the analyses | %ld MB |\n", r->rss_analysis_mb);
    fprintf(f, "| Build | %s |\n\n", e->optimised ? "optimised" : "**NOT OPTIMISED (about 1.3-1.7x slower than -O2): do not compare with optimised runs**");

    fprintf(f, "## Environment\n\n| | |\n|---|---|\n");
    fprintf(f, "| CPU | %s (%ld logical cores; the library is single-threaded) |\n", e->cpu, e->cores);
    fprintf(f, "| Memory | %.1f GB |\n| OS | %s |\n| Compiler | %s |\n| Build flags | `%s` |\n", e->ram_gb, e->kernel, e->compiler, BUILD_CFLAGS);
    fprintf(f, "| Machine load (1-minute average) | %.2f at the start, %.2f at the end (%s) |\n", e->load_start, e->load_end, e->load_start > 0.5 ? "**other programs were using the CPU: timings may be inflated; for comparisons run on an idle machine**" : "machine was idle");
    fprintf(f, "| Run | %d repeats per analysis, seed 88172645463325252, memory budget 4096 MiB, result cache disabled for the scored analyses (every run recomputes; the cache has its own section below) |\n\n", g_repeats);

    fprintf(f, "## The network\n\n| | |\n|---|---|\n");
    fprintf(f, "| Users | %zu |\n| Friendships | %zu |\n| Users with no friends | %zu (%.0f%%) |\n", r->users, r->edges, r->silent, 100.0 * (double)r->silent / (double)r->users);
    fprintf(f, "| Average / largest number of friends | %.1f / %zu |\n| Communities generated | %zu |\n| Time to build it (one run) | %.2f s |\n\n", r->avg_deg, r->maxdeg, r->ncomm_generated, r->build_s);

    fprintf(f, "## Whole-graph analyses\n\n");
    fprintf(f, "Rate = friendships / BEST time (on a shared machine the best run is the least disturbed). Peak memory is the process high-water mark so far (it never goes down).\n\n");
    fprintf(f, "| Analysis | **best (s)** | median (s) | worst (s) | M friendships/s | peak MB | independent check |%s\n", has_baseline ? " vs baseline |" : "");
    fprintf(f, "|---|---:|---:|---:|---:|---:|---|%s\n", has_baseline ? "---:|" : "");
    for (int i = 0; i < g_nstage; i++) {
        const Stage *s = &g_stage[i];
        char rate[24];
        if (s->rate > 0) snprintf(rate, sizeof rate, "%.1f", s->rate); else snprintf(rate, sizeof rate, "-");
        fprintf(f, "| %s%s | **%.3f** | %.3f | %.3f | %s | %ld | %s %s |", s->name, s->scored ? "" : " (not scored)", s->min, s->med, s->max,
                rate, s->rss_mb, s->ok ? "passed:" : "**FAILED:**", s->check);
        if (has_baseline) {
            if (s->baseline > 0) fprintf(f, " %.2fx %s |", s->baseline > s->min ? s->baseline / s->min : s->min / s->baseline, s->baseline > s->min ? "faster" : "slower");
            else fprintf(f, " n/a |");
        }
        fprintf(f, "\n");
    }
    fprintf(f, "\n");

    fprintf(f, "## Per-user queries (latency)\n\n| Query | queries | median | 99th percentile | worst | per second |\n|---|---:|---:|---:|---:|---:|\n");
    for (int i = 0; i < r->nq; i++)
        fprintf(f, "| %s | %zu | %.3f ms | %.3f ms | %.2f ms | %.0f |\n", r->q[i].name, r->q[i].n, r->q[i].p50 * 1e3, r->q[i].p99 * 1e3, r->q[i].max * 1e3, r->q[i].qps);
    fprintf(f, "\n");

    fprintf(f, "## Memory: ways to hold this network\n\n| Representation | memory | bytes per friendship |\n|---|---:|---:|\n");
    for (int i = 0; i < r->nmem; i++) fprintf(f, "| %s | %.0f MB | %.1f |\n", r->mem[i].label, r->mem[i].mb, r->mem[i].bytes_per_edge);
    fprintf(f, "\nThe dense matrix row is the size of its presence bits alone; the library refuses it (matrix limit 175,851 users unweighted).\n\n");

    fprintf(f, "## Low-memory mode (sharding engine, 64 MiB budget)\n\n| | |\n|---|---|\n");
    fprintf(f, "| Friendships cut into | %zu shard files in %.1f s |\n| Shards resident at once | %zu (%.0f MB) |\n", r->sh.shards, r->sh.cut_s, r->sh.cache, r->sh.cache_mb);
    fprintf(f, "| Connected groups | %.2f s |\n| PageRank | %.2f s |\n| Disk traffic | %zu shard loads, %.0f MB read |\n", r->sh.comp_s, r->sh.pr_s, r->sh.loads, r->sh.read_mb);
    fprintf(f, "| Agreement with the in-memory run | largest PageRank difference %.1e (%s) |\n\n", r->sh.diff, r->sh.ok ? "passed" : "FAILED");
    fprintf(f, "The operating system's file cache serves most of those reads here; on a device with a slow SD card the sharded run would be slower.\n\n");

    fprintf(f, "## Result cache (first call against repeat call)\n\nA second copy of the network with the result cache enabled (the scored graph has it off, so the analyses above time the computation). Median of %d repeats; each repeat is compared byte for byte with the first answer.\n\n| Call | first call | repeat (cache hit) | speed-up | repeat identical |\n|---|---:|---:|---:|---|\n", g_repeats);
    for (int k = 0; k < 3; k++)
        fprintf(f, "| %s | %.2f ms | %.4f ms | %.0fx | %s |\n", r->ca.c[k].name, r->ca.c[k].first_ms, r->ca.c[k].repeat_ms,
                r->ca.c[k].first_ms / (r->ca.c[k].repeat_ms > 1e-6 ? r->ca.c[k].repeat_ms : 1e-6), r->ca.c[k].same ? "yes" : "NO");
    fprintf(f, "\nCounters after the run: %llu hits, %llu misses, %llu stores, %llu invalidations, %.1f MB held. After one new friendship the next components call recomputed and found one group fewer (%s), so the cache did not serve a stale answer.\n\nA repeat is not free: the stored answer is copied out, so its cost grows with the result size (4 to 12 MB here).\n\n",
            (unsigned long long)r->ca.hits, (unsigned long long)r->ca.misses, (unsigned long long)r->ca.stores,
            (unsigned long long)r->ca.invalidations, (double)r->ca.bytes / 1048576.0, r->ca.fresh_ok ? "passed" : "FAILED");

    fprintf(f, "## Pictures (docs/images/million/)\n\n");
    for (int i = 0; i < r->npic; i++) fprintf(f, "### %s\n\n![%s](%s)\n\n%s\n\n", r->pic[i].file, r->pic[i].file, r->pic[i].file, r->pic[i].caption);

    fprintf(f, "## How to read and reproduce this\n\n");
    fprintf(f, "- **Reproduce:** `make bench-million` (rebuilds everything with `-O2` in `build-o2/`, then runs this example). Compare with an earlier run: `make bench-million BASELINE=old/results.csv`.\n");
    fprintf(f, "- **Suite score:** the geometric mean, over the scored analyses, of friendships processed per second of wall time in the best of the repeated runs (the run least disturbed by other programs). It is one comparable number; the table above is the real result. PageRank and Louvain do many sweeps per call, so their rates are lower by design.\n");
    fprintf(f, "- **Independent checks:** each analysis is verified by code in this file that uses only the generator's own edge list (union-find for groups, a recomputed modularity for communities, a fresh power-iteration step for PageRank, brute-force counting for clustering and suggestions, ...). A failed check marks the stage FAILED and the run exits non-zero.\n");
    fprintf(f, "- **Caveats:** one machine, one network from one seed; timings include nothing but the library call (no drawing, no I/O); the library is single-threaded so extra cores do not help; an unoptimised (-O0) build measured 1.3 to 1.7 times slower than -O2 on these analyses (the algorithms are memory-bound), so never compare numbers across build flags. On a shared laptop the same analysis varied by up to 2x between whole runs (five complete runs of this benchmark scored 3.7 to 7.1 million friendships per second, because other programs and CPU frequency scaling change the speed), which is why the score uses the best of %d runs and the load average is recorded; expect a spread of tens of percent unless the machine is idle\n", g_repeats);
}

static void write_csv(const char *dir)
{
    char path[600];
    snprintf(path, sizeof path, "%s/results.csv", dir);
    FILE *f = fopen(path, "w");
    if (f == NULL) return;
    fprintf(f, "stage,best_s,median_s,worst_s,repeats,m_friendships_per_s,peak_rss_mb,check_passed\n");
    for (int i = 0; i < g_nstage; i++) {
        const Stage *s = &g_stage[i];
        fprintf(f, "%s,%.6f,%.6f,%.6f,%d,%.3f,%ld,%d\n", s->name, s->min, s->med, s->max, s->reps, s->rate, s->rss_mb, s->ok);
    }
    fclose(f);
}

/* Reads an earlier results.csv and fills Stage.baseline (the median) for stages with the same name. */
static bool load_baseline(const char *path)
{
    FILE *f = fopen(path, "r");
    if (f == NULL) return false;
    char line[400];
    bool any = false;
    while (fgets(line, sizeof line, f)) {
        char *comma = strchr(line, ',');
        if (comma == NULL) continue;
        *comma = 0;
        double best, med;
        if (sscanf(comma + 1, "%lf,%lf", &best, &med) != 2) continue;
        for (int i = 0; i < g_nstage; i++) if (strcmp(g_stage[i].name, line) == 0) { g_stage[i].baseline = med; any = true; }
    }
    fclose(f);
    return any;
}

static double percentile(double *v, size_t n, double p)
{
    qsort(v, n, sizeof *v, cmp_d);
    return v[(size_t)(p * (double)(n - 1))];
}

/* ================================================================== */
/*  The run                                                             */
/* ================================================================== */

static const size_t *g_key;
static int by_key_desc(const void *a, const void *b)
{
    const size_t x = g_key[*(const int *)a], y = g_key[*(const int *)b];
    return (x < y) - (x > y);
}

static int *order_by_key(const size_t *key, size_t n)         /* user IDs sorted by key, largest first */
{
    int *idx = malloc(n * sizeof *idx);
    for (size_t i = 0; i < n; i++) idx[i] = (int)i;
    g_key = key;
    qsort(idx, n, sizeof *idx, by_key_desc);
    return idx;
}

static void must(const char *what, GraphStatus st)
{
    if (st != GRAPH_OK) { fprintf(stderr, "%s failed: %s\n", what, graph_status_string(st)); exit(1); }
}

static Rep g_rep;

static void drawn(Graph *g, const char *dir, const char *file, int size, const VizOptions *base, const char *title, const char *caption)
{
    char path[600];
    snprintf(path, sizeof path, "%s/%s", dir, file);
    VizOptions o = *base;
    o.width = o.height = size;
    o.format = VIZ_FORMAT_PNG;
    o.title  = title;
    const double t = now();
    must(file, graph_visualize(g, path, &o));
    const double secs = now() - t;
    printf("  drew %s  (%zu vertices, %zu edges, %.1f s)\n", file, graph_vertex_count(g), graph_edge_count(g), secs);
    snprintf(g_rep.pic[g_rep.npic].file, sizeof g_rep.pic[0].file, "%s", file);
    snprintf(g_rep.pic[g_rep.npic].caption, sizeof g_rep.pic[0].caption, "%s", caption);
    g_rep.pic[g_rep.npic++].secs = secs;
}

/* a chart saved by the toolkit: remember it for the report */
static void chart_note(const char *file, const char *caption)
{
    snprintf(g_rep.pic[g_rep.npic].file, sizeof g_rep.pic[0].file, "%s", file);
    snprintf(g_rep.pic[g_rep.npic].caption, sizeof g_rep.pic[0].caption, "%s", caption);
    g_rep.pic[g_rep.npic++].secs = 0;
}

int main(int argc, char **argv)
{
    const char *dir = argc > 1 ? argv[1] : "docs/images/million";
    const size_t users = argc > 2 ? (size_t)strtoull(argv[2], NULL, 10) : 1000000;
    const char *baseline = argc > 3 && argv[3][0] ? argv[3] : NULL;
    if (getenv("BENCH_REPEATS") != NULL) g_repeats = atoi(getenv("BENCH_REPEATS"));
    if (g_repeats < 1) g_repeats = 1;
    if (g_repeats > MAXREP) g_repeats = MAXREP;
    const size_t ncomm = users / 500 > 20 ? users / 500 : 20;          /* about 500 people per community on average */
    char title[200], sub[400], msg[200];
    double ts[MAXREP], t;
    Rep *R = &g_rep;
    const Env env = capture_env();

    printf("Benchmark: a social network of %zu users\n", users);
    printf("  cpu: %s | %ld cores | %.1f GB | %s\n  build: %s | %s\n", env.cpu, env.cores, env.ram_gb, env.kernel, BUILD_CFLAGS, env.optimised ? "optimised" : "NOT OPTIMISED");
    if (env.load_start > 0.5) printf("  *** load average %.2f: other programs are using this machine, timings may be inflated ***\n", env.load_start);
    if (!env.optimised) printf("  *** this build is not optimised (about 1.3-1.7x slower than -O2): do not compare with optimised runs; use `make bench-million` ***\n");
    printf("  %d repeat(s) per analysis\n\n", g_repeats);

    /* ---------------------------------------------------------------- the network */
    int *comm = malloc(users * sizeof *comm);
    Vec eu = {0}, ev = {0};
    t = now();
    Graph *g = make_network(users, ncomm, comm, &eu, &ev);
    R->build_s = now() - t;
    const size_t edges = graph_edge_count(g);
    size_t maxdeg = 0, hub = 0, silent = 0;
    for (size_t v = 0; v < users; v++) {
        const size_t d = graph_degree(g, (int)v);
        if (d > maxdeg) { maxdeg = d; hub = v; }
        silent += d == 0;
    }
    R->users = users; R->edges = edges; R->silent = silent; R->maxdeg = maxdeg; R->hub = hub; R->ncomm_generated = ncomm;
    R->avg_deg = 2.0 * (double)edges / (double)users;
    printf("The network: %zu users, %zu friendships, %zu with no friends, built in %.1f s; most friends: %zu (user %zu)\n\n", users, edges, silent, R->build_s, maxdeg, hub);

    printf("Analyses (best / median / worst of %d)\n", g_repeats);

    /* ---------------------------------------------------------------- degrees */
    size_t *hist = calloc(maxdeg + 2, sizeof *hist);
    REPEAT(ts, g_repeats, must("degree histogram", graph_degree_histogram(g, hist, maxdeg + 2)));
    {
        size_t sum = 0, ends_sum = 0;
        for (size_t d = 0; d <= maxdeg; d++) { sum += hist[d]; ends_sum += d * hist[d]; }
        snprintf(msg, sizeof msg, "users sum to %zu, friendship ends sum to 2 x %zu, silent users = %zu", sum, edges, hist[0]);
        stage_done("degree histogram", ts, g_repeats, edges, true, sum == users && ends_sum == 2 * edges && hist[0] == silent, msg);
    }
    size_t cum = 0, median = 0, p99 = 0;
    for (size_t d = 0; d <= maxdeg; d++) {
        cum += hist[d];
        if (!median && cum * 2 >= users) median = d;
        if (!p99 && cum * 100 >= users * 99) p99 = d;
    }
    {
        double *xs = malloc((maxdeg + 1) * sizeof *xs), *ys = malloc((maxdeg + 1) * sizeof *ys);
        size_t k = 0;
        for (size_t d = 1; d <= maxdeg; d++) if (hist[d] > 0) { xs[k] = (double)d; ys[k++] = (double)hist[d]; }
        snprintf(sub, sizeof sub, "%zu users: %zu with no friends (not shown), median %zu, 99th percentile %zu, most popular %zu. Heavy-tailed: a few hubs have many times the typical number of friends (this model is not a clean power law).", users, hist[0], median, p99, maxdeg);
        chart_points(dir, "01_degree_distribution.png", "How many friends does a user have?", sub, "friends (degree)", "number of users", xs, ys, k, true, true, false);
        chart_note("01_degree_distribution.png", "Friends per user on log-log axes. Heavy-tailed: most users have a handful of friends and a few hubs have hundreds. The curve bends down faster than a power law would, because this generated network caps how many links a newcomer makes.");
        free(xs); free(ys);
    }

    /* ---------------------------------------------------------------- connected groups */
    int *comp = malloc(users * sizeof *comp);
    size_t ngroups = 0;
    REPEAT(ts, g_repeats, ngroups = graph_connected_components(g, comp));
    stage_done("connected components", ts, g_repeats, edges, true, verify_components(users, &eu, &ev, comp, ngroups, msg, sizeof msg), msg);
    {
        size_t *gsize = calloc(ngroups, sizeof *gsize);
        for (size_t v = 0; v < users; v++) gsize[comp[v]]++;
        size_t giant = 0, c_iso = 0, c_small = 0, c_mid = 0, c_big = 0, u_iso = 0, u_small = 0, u_mid = 0, u_big = 0;
        for (size_t c = 0; c < ngroups; c++) if (gsize[c] > giant) giant = gsize[c];
        for (size_t c = 0; c < ngroups; c++) {
            const size_t sz = gsize[c];
            if (sz == giant) continue;
            if (sz == 1) { c_iso++; u_iso += sz; } else if (sz < 10) { c_small++; u_small += sz; }
            else if (sz < 100) { c_mid++; u_mid += sz; } else { c_big++; u_big += sz; }
        }
        const char *labels[5] = {"alone (no friends)", "2-9 users", "10-99 users", "100+ (not giant)", "the giant group"};
        const double v[5] = {(double)u_iso, (double)u_small, (double)u_mid, (double)u_big, (double)giant};
        snprintf(sub, sizeof sub, "%zu connected groups. Bars show how many USERS live in each kind of group (log scale); the groups: %zu alone, %zu small, %zu medium, %zu large, 1 giant.", ngroups, c_iso, c_small, c_mid, c_big);
        chart_bars(dir, "02_connected_groups.png", "Who is connected to whom: one giant group, and a long tail", sub, "kind of group", "users (log scale)", labels, v, 5, true, 4);
        chart_note("02_connected_groups.png", "How many users live in each kind of connected group. The giant group holds nearly everyone who has a friend; the silent users are groups of one.");
        free(gsize);
    }

    /* ---------------------------------------------------------------- k-core */
    int *core = malloc(users * sizeof *core);
    REPEAT(ts, g_repeats, must("core numbers", graph_core_numbers(g, core)));
    stage_done("k-core decomposition", ts, g_repeats, edges, true, verify_core(g, users, core, msg, sizeof msg), msg);
    {
        int maxcore = 0;
        for (size_t v = 0; v < users; v++) if (core[v] > maxcore) maxcore = core[v];
        const size_t nb = (size_t)maxcore + 1;
        double *v = calloc(nb, sizeof *v);
        char (*lab)[24] = malloc(nb * sizeof *lab);
        const char **lp = malloc(nb * sizeof *lp);
        for (size_t u = 0; u < users; u++) v[core[u]]++;
        for (size_t i = 0; i < nb; i++) { snprintf(lab[i], 24, "%zu", i); lp[i] = lab[i]; }
        snprintf(sub, sizeof sub, "A user is in the k-core if they keep at least k friends after repeatedly removing everyone with fewer. Deepest core: %d.", maxcore);
        chart_bars(dir, "03_k_core.png", "How deeply embedded are users?", sub, "core number k", "users (log scale)", lp, v, nb, true, maxcore);
        chart_note("03_k_core.png", "Users per core number. Most users sit in the shallow cores; the highlighted bar is the deepest, most tightly knit core.");
        free(v); free(lab); free(lp);
    }

    /* ---------------------------------------------------------------- clustering */
    double *local = malloc(users * sizeof *local), avgc = 0;
    long long triangles = 0;
    REPEAT(ts, g_repeats, { must("clustering", graph_clustering_coefficient(g, local, &avgc)); must("triangles", graph_triangle_count(g, &triangles)); });
    stage_done("triangles + clustering", ts, g_repeats, edges, true, verify_clustering(g, users, local, triangles, msg, sizeof msg), msg);
    {
        enum { NB = 12 };
        double sum[NB] = {0}, cnt[NB] = {0}, xs[NB], ys[NB];
        for (size_t v = 0; v < users; v++) {
            const size_t d = graph_degree(g, (int)v);
            if (d < 2) continue;
            int b = (int)floor(log2((double)d));
            if (b >= NB) b = NB - 1;
            sum[b] += local[v]; cnt[b]++;
        }
        size_t k = 0;
        for (int b = 0; b < NB; b++) if (cnt[b] > 0) { xs[k] = pow(2, b); ys[k++] = sum[b] / cnt[b]; }
        snprintf(sub, sizeof sub, "%lld triangles; average over all users %.3f. Popular users sit in looser neighbourhoods: the more friends, the smaller the share of them who know each other.", triangles, avgc);
        chart_points(dir, "04_clustering_by_degree.png", "Do a user's friends know each other?", sub, "friends (degree, buckets of powers of two)", "average clustering coefficient", xs, ys, k, true, false, true);
        chart_note("04_clustering_by_degree.png", "Average clustering coefficient by number of friends. It falls as degree rises: popular users connect people who do not know each other.");
    }
    free(local);

    /* ---------------------------------------------------------------- PageRank */
    double *rank = malloc(users * sizeof *rank);
    size_t rounds = 0;
    REPEAT(ts, g_repeats, must("pagerank", graph_pagerank(g, 0.85, 1e-8, 100, rank, &rounds)));
    stage_done("PageRank (to 1e-8)", ts, g_repeats, edges, true, verify_pagerank(users, &eu, &ev, rank, msg, sizeof msg), msg);
    size_t *rkey = malloc(users * sizeof *rkey);
    for (size_t v = 0; v < users; v++) rkey[v] = (size_t)(rank[v] * 1e12);
    int *top = order_by_key(rkey, users);
    free(rkey);
    {
        enum { TOPN = 20 };
        char (*lab)[56] = malloc(TOPN * sizeof *lab);
        const char *lp[TOPN];
        double v[TOPN];
        for (int i = 0; i < TOPN; i++) { snprintf(lab[i], 56, "user %d (%zu friends)", top[i], graph_degree(g, top[i])); lp[i] = lab[i]; v[i] = rank[top[i]] * (double)users; }
        snprintf(sub, sizeof sub, "PageRank (%zu rounds) of the 20 most influential users, as a multiple of the average user's rank. Influence follows friend count but not exactly.", rounds);
        chart_hbars(dir, "05_pagerank_top20.png", "The 20 most influential users", sub, "x average", lp, v, TOPN, false, 0);
        chart_note("05_pagerank_top20.png", "The 20 users with the highest PageRank, each shown with how many friends they have.");
        free(lab);
        const size_t sample = users < 40000 ? users : 40000;
        double *xs = malloc(sample * sizeof *xs), *ys = malloc(sample * sizeof *ys);
        size_t k = 0;
        for (size_t i = 0; i < sample; i++) {
            const size_t u = (size_t)((double)i / (double)sample * (double)users);
            const size_t d = graph_degree(g, (int)u);
            if (d == 0) continue;
            xs[k] = (double)d; ys[k++] = rank[u] * (double)users;
        }
        snprintf(sub, sizeof sub, "%zu randomly chosen users with friends. Rank grows with friend count, with a band of spread around the trend.", k);
        chart_points(dir, "06_pagerank_vs_degree.png", "Friends versus influence", sub, "friends (degree)", "PageRank (x average)", xs, ys, k, true, true, false);
        chart_note("06_pagerank_vs_degree.png", "Each dot is a user. Influence rises with friend count, but at any friend count it varies: who your friends are matters.");
        free(xs); free(ys);
    }

    /* ---------------------------------------------------------------- communities */
    int *lc = malloc(users * sizeof *lc);
    double modularity = 0;
    size_t ncl = 0;
    REPEAT(ts, g_repeats, must("louvain", graph_louvain(g, lc, &modularity, &ncl)));
    stage_done("Louvain communities", ts, g_repeats, edges, true, verify_louvain(users, &eu, &ev, lc, ncl, modularity, msg, sizeof msg), msg);
    size_t *csize = calloc(ncl, sizeof *csize);
    for (size_t v = 0; v < users; v++) if (lc[v] >= 0) csize[lc[v]]++;
    int *corder = order_by_key(csize, ncl);
    {
        enum { TOPC = 25 };
        double v[TOPC];
        char (*lab)[24] = malloc(TOPC * sizeof *lab);
        const char *lp[TOPC];
        const size_t show = ncl < TOPC ? ncl : TOPC;
        for (size_t i = 0; i < show; i++) { v[i] = (double)csize[corder[i]]; snprintf(lab[i], 24, "#%zu", i + 1); lp[i] = lab[i]; }
        snprintf(sub, sizeof sub, "%zu communities found; modularity %.3f (above 0.3 means real community structure). The 25 largest are shown.", ncl, modularity);
        chart_bars(dir, "07_largest_communities.png", "The 25 largest communities", sub, "community, by size", "users", lp, v, show, false, -1);
        chart_note("07_largest_communities.png", "Sizes of the 25 largest communities Louvain found.");
        free(lab);
        const size_t maxsz = csize[corder[0]];
        size_t *ofsize = calloc(maxsz + 1, sizeof *ofsize);
        for (size_t c = 0; c < ncl; c++) if (csize[c] > 0) ofsize[csize[c]]++;
        double *xs = malloc((maxsz + 1) * sizeof *xs), *ys = malloc((maxsz + 1) * sizeof *ys);
        size_t k = 0;
        for (size_t s2 = 1; s2 <= maxsz; s2++) if (ofsize[s2]) { xs[k] = (double)s2; ys[k++] = (double)ofsize[s2]; }
        chart_points(dir, "08_community_sizes.png", "How big are communities?", "Most communities are small and a few are large: how many communities have each size.", "users in the community", "number of communities", xs, ys, k, true, true, false);
        chart_note("08_community_sizes.png", "Number of communities of each size, log-log: many small communities, few large ones.");
        free(ofsize); free(xs); free(ys);
    }

    /* ---------------------------------------------------------------- reach */
    double *dist = malloc(users * sizeof *dist);
    int *par = malloc(users * sizeof *par);
    REPEAT(ts, g_repeats, graph_dijkstra(g, (int)hub, dist, par));       /* every friendship costs 1, so distance = hops */
    stage_done("hop distances from the top hub", ts, g_repeats, edges, true, verify_hops(&eu, &ev, dist, par, hub, users, msg, sizeof msg), msg);
    {
        size_t atdist[40] = {0}, reach = 0;
        int lastd = 0;
        for (size_t v = 0; v < users; v++) if (dist[v] < 39) { atdist[(int)dist[v]]++; reach++; }
        for (int d = 0; d < 40; d++) if (atdist[d]) lastd = d;
        char (*lab)[24] = malloc(40 * sizeof *lab);
        const char *lp[40];
        double v[40];
        for (int d = 0; d <= lastd; d++) { snprintf(lab[d], 24, "%d", d); lp[d] = lab[d]; v[d] = (double)atdist[d]; }
        snprintf(sub, sizeof sub, "How many users are exactly d friendship-steps from user %zu, the most-connected user. %zu of %zu users (%.0f%%) can be reached at all.", hub, reach, users, 100.0 * (double)reach / (double)users);
        chart_bars(dir, "09_hops.png", "Six degrees of separation", sub, "hops from the top hub", "users", lp, v, (size_t)lastd + 1, false, -1);
        chart_note("09_hops.png", "Users at each distance from the most-connected user: almost everyone is within a handful of hops.");
        free(lab);
    }

    R->rss_analysis_mb = rss_peak_mb();
    printf("\n");

    /* ---------------------------------------------------------------- per-user queries: latency */
    printf("Per-user queries\n");
    {
        const size_t NQ = 1000, NP = 300;
        double *tq = malloc(NQ * sizeof *tq), *tp = malloc(NP * sizeof *tp);
        GraphSuggestion sug[5];
        size_t count = 0, done = 0, found = 0;
        double total = 0;
        for (size_t i = 0; i < NQ; i++) {                                /* friend suggestions for random active users */
            int u;
            do u = (int)rnd((unsigned)users); while (graph_degree(g, u) == 0);
            t = now();
            graph_social_suggest_friends(g, u, 5, sug, &count);
            tq[done++] = now() - t;
            total += tq[done - 1];
        }
        R->q[R->nq].n = NQ; R->q[R->nq].ok = NQ; R->q[R->nq].qps = (double)NQ / total;
        R->q[R->nq].max = 0;
        for (size_t i = 0; i < NQ; i++) if (tq[i] > R->q[R->nq].max) R->q[R->nq].max = tq[i];
        R->q[R->nq].p50 = percentile(tq, NQ, 0.5); R->q[R->nq].p99 = percentile(tq, NQ, 0.99);
        snprintf(R->q[R->nq].name, sizeof R->q[0].name, "people you may know (random active user)");
        R->nq++;
        total = 0; done = 0;
        for (size_t i = 0; i < NP; i++) {                                /* shortest chain between two random active users */
            int a, b;
            do a = (int)rnd((unsigned)users); while (graph_degree(g, a) == 0);
            do b = (int)rnd((unsigned)users); while (graph_degree(g, b) == 0);
            GraphPath path = {0};
            t = now();
            const GraphStatus ps = graph_social_path(g, a, b, SIZE_MAX, &path);
            tp[done++] = now() - t;
            total += tp[done - 1];
            found += ps == GRAPH_OK;
            graph_path_free(&path);
        }
        R->q[R->nq].n = NP; R->q[R->nq].ok = found; R->q[R->nq].qps = (double)NP / total;
        R->q[R->nq].max = 0;
        for (size_t i = 0; i < NP; i++) if (tp[i] > R->q[R->nq].max) R->q[R->nq].max = tp[i];
        R->q[R->nq].p50 = percentile(tp, NP, 0.5); R->q[R->nq].p99 = percentile(tp, NP, 0.99);
        snprintf(R->q[R->nq].name, sizeof R->q[0].name, "shortest chain of friends (%zu of %zu pairs connected)", found, NP);
        R->nq++;
        for (int i = 0; i < R->nq; i++) printf("  %-62s median %.3f ms, 99th percentile %.3f ms, %.0f per second\n", R->q[i].name, R->q[i].p50 * 1e3, R->q[i].p99 * 1e3, R->q[i].qps);
        ts[0] = R->q[0].p50;
        stage_done("friend suggestion (median query)", ts, 1, 0, false, verify_suggestions(g, users, msg, sizeof msg), msg);
        free(tq); free(tp);
    }
    printf("\n");

    /* ---------------------------------------------------------------- memory */
    printf("Memory\n");
    {
        GraphConfig wc = graph_default_config(false);
        wc.adaptive = false;
        Graph *w = graph_create_with(&wc, users);
        for (size_t i = 0; i < eu.n; i++) graph_add_edge(w, eu.a[i], ev.a[i], 1.0 + (double)(i % 5));
        const double mw = (double)w->be->memory_bytes(w->impl);
        graph_destroy(w);
        const double mu = (double)g->be->memory_bytes(g->impl);
        graph_shrink_to_fit(g);
        const double ms = (double)g->be->memory_bytes(g->impl);
        const double dense = (double)users * (double)users / 8.0;
        const struct { const char *label; double bytes; } rows[4] = {{"weighted list", mw}, {"unweighted list", mu}, {"unweighted list, shrunk", ms}, {"dense matrix (refused)", dense}};
        for (int i = 0; i < 4; i++) {
            snprintf(R->mem[R->nmem].label, sizeof R->mem[0].label, "%s", rows[i].label);
            R->mem[R->nmem].mb = rows[i].bytes / 1048576.0;
            R->mem[R->nmem++].bytes_per_edge = rows[i].bytes / (double)edges;
            printf("  %-46s %10.0f MB  %9.1f bytes per friendship\n", rows[i].label, rows[i].bytes / 1048576.0, rows[i].bytes / (double)edges);
        }
        Graph *probe = graph_copy(g);
        printf("  converting to a matrix: %s\n\n", graph_convert(probe, GRAPH_MATRIX) ? "unexpectedly allowed" : graph_status_string(graph_last_error()));
        graph_destroy(probe);
    }


    /* ---------------------------------------------------------------- the sharding engine */
    printf("Low-memory mode: the sharding engine on a 64 MiB budget\n");
    {
        FILE *f = tmpfile();
        for (size_t i = 0; i < eu.n; i++) fprintf(f, "%d %d\n", eu.a[i], ev.a[i]);
        fprintf(f, "# vertices %zu\n", users);
        fflush(f);
        graph_set_memory_budget((size_t)64 << 20);
        ShardConfig cfg = shard_default_config();
        cfg.force = true;
        Graph *none = NULL;
        ShardedGraph *sg = NULL;
        t = now();
        must("shard", graph_open_edge_list(f, &cfg, &none, &sg));
        R->sh.cut_s = now() - t;
        size_t sgroups = 0, sr = 0;
        double *srank = malloc(users * sizeof *srank);
        int *scomp = malloc(users * sizeof *scomp);
        t = now(); must("shard components", shard_components(sg, scomp, &sgroups)); R->sh.comp_s = now() - t;
        t = now(); must("shard pagerank", shard_pagerank(sg, 0.85, 1e-8, 100, srank, &sr)); R->sh.pr_s = now() - t;
        double worst = 0;
        for (size_t v = 0; v < users; v++) worst = fmax(worst, fabs(srank[v] - rank[v]));
        const ShardStats st = shard_stats(sg);
        R->sh.shards = st.shards; R->sh.cache = st.cache_capacity; R->sh.loads = st.loads;
        R->sh.cache_mb = (double)(st.cache_capacity * st.shard_bytes_max) / 1048576.0;
        R->sh.read_mb = (double)st.bytes_read / 1048576.0; R->sh.diff = worst;
        R->sh.ok = sgroups == ngroups && worst < 1e-6;
        printf("  %zu friendships cut into %zu shard files in %.1f s; cache holds %zu (%.0f MB)\n", edges, st.shards, R->sh.cut_s, st.cache_capacity, R->sh.cache_mb);
        printf("  groups %.2f s (%zu, in memory %zu); PageRank %.2f s (%zu rounds), largest difference from in memory %.1e; %zu loads, %.0f MB read\n\n", R->sh.comp_s, sgroups, ngroups, R->sh.pr_s, sr, worst, st.loads, R->sh.read_mb);
        snprintf(R->mem[R->nmem].label, sizeof R->mem[0].label, "sharded: cache in RAM");
        R->mem[R->nmem].mb = R->sh.cache_mb;
        R->mem[R->nmem++].bytes_per_edge = R->sh.cache_mb * 1048576.0 / (double)edges;
        shard_destroy(sg);
        free(srank); free(scomp);
        fclose(f);
        graph_set_memory_budget(GRAPH_DEFAULT_MEMORY_BUDGET);
    }

    /* ---------------------------------------------------------------- the result cache as a feature */
    printf("Result cache: first call against repeat call\n");
    {
        /* A second copy of the network with the cache ON (the scored graph above has it off, so its timings stay
         * honest).  The default 4 MiB limit is too small for 1,000,000-entry results, so raise it. */
        GraphConfig cc = graph_default_config(false);
        cc.unweighted = true;
        cc.adaptive = false;
        cc.cache_entries = 32;
        cc.cache_bytes = (size_t)256 << 20;
        Graph *cg = graph_create_with(&cc, users);
        if (cg == NULL) { fprintf(stderr, "cannot create the cached graph\n"); return 1; }
        for (size_t i = 0; i < eu.n; i++) graph_add_edge(cg, eu.a[i], ev.a[i], 1.0);
        int *o1 = malloc(users * sizeof *o1), *o2 = malloc(users * sizeof *o2);
        int *p1 = malloc(users * sizeof *p1), *p2 = malloc(users * sizeof *p2);
        double *d1 = malloc(users * sizeof *d1), *d2 = malloc(users * sizeof *d2);
        const char *names[3] = { "connected components", "BFS from the most popular user", "Dijkstra from the most popular user" };
        size_t ngroups = 0;
        R->ca.ok = true;
        for (int k = 0; k < 3; k++) {
            snprintf(R->ca.c[k].name, sizeof R->ca.c[k].name, "%s", names[k]);
            t = now();
            if (k == 0) ngroups = graph_connected_components(cg, o1);
            else if (k == 1) (void)graph_bfs(cg, (int)hub, o1, p1);
            else (void)graph_dijkstra(cg, (int)hub, d1, p1);
            R->ca.c[k].first_ms = (now() - t) * 1e3;
            bool same = true;
            for (int r = 0; r < g_repeats; r++) {
                t = now();
                if (k == 0) (void)graph_connected_components(cg, o2);
                else if (k == 1) (void)graph_bfs(cg, (int)hub, o2, p2);
                else (void)graph_dijkstra(cg, (int)hub, d2, p2);
                ts[r] = (now() - t) * 1e3;
                /* a fast wrong answer must not count: the repeat has to match the first call byte for byte */
                if (k == 0) same = same && memcmp(o1, o2, users * sizeof *o1) == 0;
                else if (k == 1) same = same && memcmp(o1, o2, users * sizeof *o1) == 0 && memcmp(p1, p2, users * sizeof *p1) == 0;
                else same = same && memcmp(d1, d2, users * sizeof *d1) == 0 && memcmp(p1, p2, users * sizeof *p1) == 0;
            }
            qsort(ts, (size_t)g_repeats, sizeof *ts, cmp_d);
            R->ca.c[k].repeat_ms = ts[g_repeats / 2];
            R->ca.c[k].same = same;
            R->ca.ok = R->ca.ok && same;
            printf("  %-38s first %9.2f ms   repeat %9.4f ms   %8.0fx   identical: %s\n", names[k], R->ca.c[k].first_ms,
                   R->ca.c[k].repeat_ms, R->ca.c[k].first_ms / (R->ca.c[k].repeat_ms > 1e-6 ? R->ca.c[k].repeat_ms : 1e-6), same ? "yes" : "NO");
        }
        /* The cache must never serve a stale answer: join two silent users (they had no friends) and the number of
         * groups has to fall by exactly one, which can only happen if the old entry was dropped. */
        CacheStats before, after;
        graph_cache_stats(cg, &before);
        int s1 = -1, s2 = -1;
        for (size_t v = 4; v < users && s2 < 0; v += 5) { if (graph_degree(cg, (int)v) != 0) continue; if (s1 < 0) s1 = (int)v; else s2 = (int)v; }
        R->ca.fresh_ok = false;
        if (s2 >= 0 && graph_add_edge(cg, s1, s2, 1.0)) {
            const size_t again = graph_connected_components(cg, o2);
            graph_cache_stats(cg, &after);
            R->ca.fresh_ok = again + 1 == ngroups && after.misses == before.misses + 1 && after.invalidations > before.invalidations;
        }
        R->ca.ok = R->ca.ok && R->ca.fresh_ok;
        graph_cache_stats(cg, &after);
        R->ca.hits = after.hits; R->ca.misses = after.misses; R->ca.stores = after.stores;
        R->ca.invalidations = after.invalidations; R->ca.bytes = after.bytes;
        /* with g_repeats repeats per call, exactly 3 misses (the first calls) and 3 * g_repeats hits must precede the edit */
        R->ca.ok = R->ca.ok && before.misses == 3 && before.hits == (uint64_t)(3 * g_repeats);
        printf("  after one new friendship the repeat recomputes (fresh: %s); hits %llu, misses %llu, stores %llu, invalidations %llu\n\n",
               R->ca.fresh_ok ? "yes" : "NO", (unsigned long long)after.hits, (unsigned long long)after.misses,
               (unsigned long long)after.stores, (unsigned long long)after.invalidations);
        const char *cl[6]; double cv[6];
        char lab[6][64];
        for (int k = 0; k < 3; k++) {
            snprintf(lab[2 * k], sizeof lab[0], "%s: first call", names[k]);
            snprintf(lab[2 * k + 1], sizeof lab[0], "%s: repeat", names[k]);
            cl[2 * k] = lab[2 * k]; cl[2 * k + 1] = lab[2 * k + 1];
            cv[2 * k] = R->ca.c[k].first_ms; cv[2 * k + 1] = R->ca.c[k].repeat_ms;
        }
        chart_hbars(dir, "18_result_cache.png", "What does the result cache save?", "Milliseconds for one call on 1,000,000 users (log scale). The first call computes; a repeat on the unchanged graph is served from the cache.", "ms", cl, cv, 6, true, -1);
        chart_note("18_result_cache.png", "First call against repeat call for three cached analyses; the repeat copies the stored answer out instead of recomputing.");
        free(o1); free(o2); free(p1); free(p2); free(d1); free(d2);
        graph_destroy(cg);
    }

    /* ---------------------------------------------------------------- memory and time charts */
    {
        const char *ml[6]; double mv[6];
        for (int i = 0; i < R->nmem; i++) { ml[i] = R->mem[i].label; mv[i] = R->mem[i].bytes_per_edge; }
        chart_hbars(dir, "10_memory.png", "What does it cost to hold the network?", "Bytes of memory per friendship (log scale). The dense matrix is the textbook way and needs terabytes of presence bits.", "bytes", ml, mv, (size_t)R->nmem, true, 2);
        chart_note("10_memory.png", "Memory per friendship for each way of holding the network, on a log scale.");
        const char *tl[MAXSTAGE]; double tv[MAXSTAGE];
        int nt = 0;
        for (int i = 0; i < g_nstage; i++) if (g_stage[i].scored) { tl[nt] = g_stage[i].name; tv[nt++] = g_stage[i].med; }
        tl[nt] = "sharded: PageRank (64 MiB budget)"; tv[nt++] = R->sh.pr_s;
        chart_hbars(dir, "11_timings.png", "How long does each analysis take?", "Median seconds for one call on all users (log scale), built with the optimiser on.", "s", tl, tv, (size_t)nt, true, -1);
        chart_note("11_timings.png", "Median time of each whole-graph analysis; the last bar is PageRank in the 64 MiB sharded mode.");
    }

    /* ---------------------------------------------------------------- pictures of the network itself */
    printf("Pictures of the network\n");
    VizOptions o = graph_viz_default_options();
    o.group_name = "Community";
    {   /* the most popular user and their friends */
        const size_t d0 = graph_degree(g, (int)hub);
        int *nb = malloc(d0 * sizeof *nb);
        graph_neighbors(g, (int)hub, nb, NULL, d0);
        size_t *dk = malloc(d0 * sizeof *dk);
        for (size_t i = 0; i < d0; i++) dk[i] = graph_degree(g, nb[i]);
        int *idx = order_by_key(dk, d0);
        const size_t keep = d0 < 250 ? d0 : 250;
        int *vs = malloc((keep + 1) * sizeof *vs);
        vs[0] = (int)hub;
        for (size_t i = 0; i < keep; i++) vs[i + 1] = nb[idx[i]];
        Graph *ego = NULL;
        must("ego", graph_induced_subgraph(g, vs, keep + 1, &ego));
        char (*lab)[24] = malloc((keep + 1) * sizeof *lab);
        const char **lp = malloc((keep + 1) * sizeof *lp);
        for (size_t i = 0; i <= keep; i++) { snprintf(lab[i], 24, "%d", vs[i]); lp[i] = lab[i]; }
        o.labels = lp;
        o.color = VIZ_COLOR_COMMUNITY;
        o.label_top = 12;
        o.iterations = 800;                      /* a hub pulls 250 friends hard; the default settles them against the frame */
        snprintf(title, sizeof title, "User %zu (%zu friends) and their %zu best-connected friends", hub, d0, keep);
        drawn(ego, dir, "12_ego_network.png", 1100, &o, title, "The most popular user and the friends who themselves have the most friends, coloured by community.");
        o.labels = NULL;
        o.iterations = 0;
        graph_destroy(ego);
        free(nb); free(dk); free(idx); free(vs); free(lab); free(lp);
    }
    {   /* the backbone: the best-connected 3000 users */
        size_t *deg = malloc(users * sizeof *deg);
        for (size_t v = 0; v < users; v++) deg[v] = graph_degree(g, (int)v);
        int *by_deg = order_by_key(deg, users);
        size_t keep = users < 3000 ? users : 3000;
        Graph *bb = NULL;
        must("backbone", graph_induced_subgraph(g, by_deg, keep, &bb));
        {   /* Many of the best-connected users are not friends with each other.  The layout would scatter those
             * loose users as isolated dots and squash the real core into a corner, so keep the largest connected piece. */
            int *bc = malloc(keep * sizeof *bc);
            const size_t nbc = graph_connected_components(bb, bc);
            size_t *bsz = calloc(nbc, sizeof *bsz), best = 0, cnt = 0;
            for (size_t v = 0; v < keep; v++) bsz[bc[v]]++;
            for (size_t c = 1; c < nbc; c++) if (bsz[c] > bsz[best]) best = c;
            int *vs2 = malloc(keep * sizeof *vs2);
            for (size_t v = 0; v < keep; v++) if ((size_t)bc[v] == best) vs2[cnt++] = (int)v;
            Graph *core_g = NULL;
            must("backbone core", graph_induced_subgraph(bb, vs2, cnt, &core_g));
            graph_destroy(bb);
            bb = core_g;
            keep = cnt;
            free(bc); free(bsz); free(vs2);
        }
        o.color = VIZ_COLOR_COMMUNITY;
        o.label_top = -1;
        snprintf(title, sizeof title, "The connected core of the best-connected users (%zu)", keep);
        drawn(bb, dir, "13_backbone_communities.png", 1100, &o, title, "The 3000 users with the most friends, keeping the largest group of them that are friends with each other (the others have no friend among the 3000 and would only appear as scattered dots): the network's backbone, coloured by community.");
        o.color = VIZ_COLOR_DEGREE;
        snprintf(title, sizeof title, "The same %zu users: red = most friends", keep);
        drawn(bb, dir, "14_backbone_degree.png", 1100, &o, title, "The same backbone coloured by popularity: the hubs sit at the centre.");
        int *cm = malloc(keep * sizeof *cm), *newid = malloc(keep * sizeof *newid);
        size_t nc2 = 0;
        must("louvain backbone", graph_louvain(bb, cm, NULL, &nc2));
        size_t next = 0;
        for (size_t c = 0; c < nc2; c++) for (size_t v = 0; v < keep; v++) if ((size_t)cm[v] == c) newid[v] = (int)next++;
        GraphConfig cfg2 = graph_default_config(false);
        cfg2.adaptive = false;
        Graph *sorted = graph_create_with(&cfg2, keep);
        for (size_t v = 0; v < keep; v++) {
            const size_t d = graph_degree(bb, (int)v);
            int *nb = malloc((d + 1) * sizeof *nb);
            graph_neighbors(bb, (int)v, nb, NULL, d);
            for (size_t i = 0; i < d; i++) if ((int)v < nb[i]) graph_add_edge(sorted, newid[v], newid[nb[i]], 1.0);
            free(nb);
        }
        VizOptions mo = graph_viz_default_options();
        mo.layout = VIZ_LAYOUT_MATRIX;
        drawn(sorted, dir, "15_backbone_matrix.png", 1000, &mo, "The backbone as a matrix, grouped by community", "The same users as an adjacency matrix, renumbered so each community is a block on the diagonal.");
        graph_destroy(sorted); graph_destroy(bb);
        free(cm); free(newid); free(deg); free(by_deg);
    }
    {   /* the communities as a graph: one dot per community, a line where friendships cross */
        const size_t K = ncl < 250 ? ncl : 250;
        int *rankof = malloc(ncl * sizeof *rankof);
        for (size_t c = 0; c < ncl; c++) rankof[c] = -1;
        for (size_t i = 0; i < K; i++) rankof[corder[i]] = (int)i;
        unsigned *cross = calloc(K * K, sizeof *cross);
        for (size_t i = 0; i < eu.n; i++) {
            const int a = rankof[lc[eu.a[i]]], b = rankof[lc[ev.a[i]]];
            if (a >= 0 && b >= 0 && a != b) { cross[(size_t)a * K + (size_t)b]++; cross[(size_t)b * K + (size_t)a]++; }
        }
        GraphConfig cfg3 = graph_default_config(false);
        cfg3.adaptive = false;
        Graph *q = graph_create_with(&cfg3, K);
        /* Any single threshold fails: with 15% of friendships crossing communities, almost every pair has some.
         * Keep, for each community, only its 3 strongest ties (the pairs with the most friendships between them). */
        for (size_t a = 0; a < K; a++)
            for (int pick = 0; pick < 3; pick++) {
                size_t best = K;
                unsigned bestc = 0;
                for (size_t b = 0; b < K; b++) {
                    if (b == a || cross[a * K + b] <= bestc || graph_has_edge(q, (int)a, (int)b, NULL)) continue;
                    best = b; bestc = cross[a * K + b];
                }
                if (best < K) graph_add_edge(q, (int)a, (int)best, (double)bestc);
            }
        o.color = VIZ_COLOR_DEGREE;
        o.label_top = 0;
        snprintf(title, sizeof title, "Each community and its 3 strongest ties");
        drawn(q, dir, "16_community_graph.png", 1000, &o, title, "One dot per community (the 250 largest); each is joined to the 3 communities it shares the most friendships with. Cross-community friendships favour already popular users, so the biggest communities become hubs of this graph.");
        graph_destroy(q); free(cross); free(rankof);
    }
    {   /* everyone who has at least one friend */
        int *keep = malloc(users * sizeof *keep);
        size_t kc = 0;
        for (size_t v = 0; v < users; v++) if (graph_degree(g, (int)v) > 0) keep[kc++] = (int)v;
        Graph *people = NULL;
        must("people", graph_induced_subgraph(g, keep, kc, &people));
        free(keep);
        o.color = VIZ_COLOR_DEGREE;
        o.label_top = -1;
        snprintf(title, sizeof title, "All %zu users who have friends: red = very popular", kc);
        drawn(people, dir, "17_whole_network.png", 1600, &o, title, "Every user who has at least one friend, drawn with the library's multilevel layout. Red = very popular.");
        graph_destroy(people);
    }
    R->rss_final_mb = rss_peak_mb();

    /* ---------------------------------------------------------------- the report */
    Env env_end = env;
    env_end.load_end = load_now();
    const bool has_baseline = baseline != NULL && load_baseline(baseline);
    char path[600];
    snprintf(path, sizeof path, "%s/REPORT.md", dir);
    FILE *rf = fopen(path, "w");
    if (rf != NULL) { emit_report(rf, &env_end, R, has_baseline); fclose(rf); }
    write_csv(dir);
    printf("\n==================================================================\n\n");
    emit_report(stdout, &env_end, R, has_baseline);
    printf("\nWrote %s/REPORT.md and %s/results.csv\n", dir, dir);
    int failed = 0;
    for (int i = 0; i < g_nstage; i++) failed += !g_stage[i].ok;
    failed += !R->sh.ok;
    failed += !R->ca.ok;
    return failed == 0 ? 0 : 1;
}
