/*
 * graph_viz.c — layout and rendering.  See graph_viz.h.
 *
 * Pipeline:  Graph --snapshot--> dense CSR (indexes 0..n-1, no dead IDs)
 *            --layout--> positions in [0,1]^2 --render--> SVG text or RGB raster --> PNG
 */

#include "graph_viz.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "graph_algorithms.h"
#include "png_writer.h"
#include "viz_font.h"

#define MAX_DIM          8192
#define MAX_RASTER_PIXELS 25000000u     /* 12 bytes per pixel in the float canvas */
#define AUTO_SVG_VERTICES 1500
#define AUTO_SVG_EDGES    8000
#define AUTO_LABELS       150
#define REPEL_CAP         12

/* ------------------------------------------------------------------ */
/*  Options                                                             */
/* ------------------------------------------------------------------ */

VizOptions graph_viz_default_options(void)
{
    VizOptions o;
    memset(&o, 0, sizeof(o));
    o.format      = VIZ_FORMAT_AUTO;
    o.layout      = VIZ_LAYOUT_FORCE;
    o.color       = VIZ_COLOR_AUTO;
    o.width       = 1200;
    o.height      = 1200;
    o.iterations  = 0;
    o.seed        = 1;
    o.show_labels = -1;
    o.show_weights = false;
    o.legend      = true;
    o.label_top   = -1;
    o.group_name  = "Group";
    return o;
}

static bool options_valid(const VizOptions *o)
{
    return o->width >= 64 && o->width <= MAX_DIM && o->height >= 64 && o->height <= MAX_DIM &&
           o->iterations >= 0;
}

/* ------------------------------------------------------------------ */
/*  Snapshot: a dense copy of the structure                             */
/* ------------------------------------------------------------------ */

typedef struct {
    size_t  n;          /* alive vertices                               */
    size_t  bound;      /* id_bound of the graph                        */
    int    *id;         /* dense index -> vertex ID                     */
    int    *idx;        /* vertex ID -> dense index, -1 if dead         */
    size_t *off;        /* n + 1 offsets into adj / w                   */
    int    *adj;        /* targets as dense indexes                     */
    double *w;          /* weights, or NULL when not requested          */
    size_t  arcs;       /* stored adjacency entries                     */
    bool    directed;
} Snap;

static void snap_free(Snap *s)
{
    free(s->id);
    free(s->idx);
    free(s->off);
    free(s->adj);
    free(s->w);
    memset(s, 0, sizeof(*s));
}

static GraphStatus snap_build(Graph *g, bool want_weights, Snap *s)
{
    memset(s, 0, sizeof(*s));
    s->directed = graph_is_directed(g);
    s->bound    = graph_id_bound(g);
    s->n        = graph_vertex_count(g);
    s->id  = malloc((s->n + 1) * sizeof(int));
    s->idx = malloc((s->bound + 1) * sizeof(int));
    s->off = calloc(s->n + 2, sizeof(size_t));
    if (s->id == NULL || s->idx == NULL || s->off == NULL) {
        snap_free(s);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    size_t k = 0, max_deg = 0, total = 0;
    for (size_t v = 0; v < s->bound; ++v) {
        if (graph_is_alive(g, (int)v)) {
            s->idx[v] = (int)k;
            s->id[k++] = (int)v;
            const size_t d = graph_degree(g, (int)v);
            max_deg = d > max_deg ? d : max_deg;
            total += d;
        } else {
            s->idx[v] = -1;
        }
    }
    s->arcs = total;
    s->adj  = malloc((total + 1) * sizeof(int));
    s->w    = want_weights ? malloc((total + 1) * sizeof(double)) : NULL;
    int    *nb = malloc((max_deg + 1) * sizeof(int));
    double *wt = malloc((max_deg + 1) * sizeof(double));
    if (s->adj == NULL || (want_weights && s->w == NULL) || nb == NULL || wt == NULL) {
        free(nb);
        free(wt);
        snap_free(s);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    size_t pos = 0;
    for (size_t i = 0; i < s->n; ++i) {
        s->off[i] = pos;
        const size_t d = graph_neighbors(g, s->id[i], nb, wt, max_deg + 1);
        for (size_t j = 0; j < d && j <= max_deg; ++j) {
            s->adj[pos] = s->idx[nb[j]];
            if (want_weights) {
                s->w[pos] = wt[j];
            }
            ++pos;
        }
    }
    s->off[s->n] = pos;
    free(nb);
    free(wt);
    return GRAPH_OK;
}

/** Edges to draw: undirected pairs once, every arc if directed; self-loops skipped. */
static bool draw_arc(const Snap *s, size_t i, int j)
{
    if ((size_t)j == i) {
        return false;
    }
    return s->directed || (size_t)j > i;
}

static size_t count_drawn(const Snap *s)
{
    size_t c = 0;
    for (size_t i = 0; i < s->n; ++i) {
        for (size_t p = s->off[i]; p < s->off[i + 1]; ++p) {
            c += draw_arc(s, i, s->adj[p]);
        }
    }
    return c;
}

/* ------------------------------------------------------------------ */
/*  Random numbers (xorshift64*)                                        */
/* ------------------------------------------------------------------ */

typedef struct { uint64_t s; } Rng;

static double rng_unit(Rng *r)       /* [0, 1) */
{
    r->s ^= r->s >> 12;
    r->s ^= r->s << 25;
    r->s ^= r->s >> 27;
    return (double)((r->s * 0x2545F4914F6CDD1DULL) >> 11) / 9007199254740992.0;
}

/* ------------------------------------------------------------------ */
/*  Force-directed layout                                               */
/* ------------------------------------------------------------------ */

/** A plain CSR adjacency (indexes 0..n-1); the layout works on one of these per level. */
typedef struct {
    size_t  n;
    size_t *off;          /* n + 1 */
    int    *adj;
} Csr;

static int auto_iterations(size_t n)
{
    return n <= 200 ? 400 : n <= 2000 ? 250 : n <= 20000 ? 120 : 60;
}

/**
 * Start positions: grow each component outward from a random root, putting
 * every vertex a short random step from the vertex that discovered it.
 */
static bool initial_positions(const Csr *c, double k, Rng *rng, double *x, double *y)
{
    int  *queue = malloc((c->n + 1) * sizeof(int));
    bool *seen  = calloc(c->n + 1, sizeof(bool));
    if (queue == NULL || seen == NULL) {
        free(queue);
        free(seen);
        return false;
    }
    for (size_t root = 0; root < c->n; ++root) {
        if (seen[root]) {
            continue;
        }
        x[root] = 0.1 + 0.8 * rng_unit(rng);
        y[root] = 0.1 + 0.8 * rng_unit(rng);
        seen[root] = true;
        size_t head = 0, tail = 0;
        queue[tail++] = (int)root;
        while (head < tail) {
            const int u = queue[head++];
            for (size_t p = c->off[u]; p < c->off[u + 1]; ++p) {
                const int v = c->adj[p];
                if (seen[v]) {
                    continue;
                }
                seen[v] = true;
                const double ang = 6.283185307179586 * rng_unit(rng);
                const double rad = k * (0.5 + rng_unit(rng));
                x[v] = fmin(1.0, fmax(0.0, x[u] + rad * cos(ang)));
                y[v] = fmin(1.0, fmax(0.0, y[u] + rad * sin(ang)));
                queue[tail++] = v;
            }
        }
    }
    free(queue);
    free(seen);
    return true;
}

/**
 * fr_run - Fruchterman-Reingold iterations on positions already in x[], y[]
 * (all kept inside [0,1]).  Repulsion only looks at nearby grid cells.
 * @t0: starting step length, cooled linearly to ~0.
 */
static GraphStatus fr_run(const Csr *c, double attract_scale, int iterations, double t0,
                          Rng *rng, double *x, double *y)
{
    const size_t n = c->n;
    if (n < 2 || iterations <= 0) {
        return GRAPH_OK;
    }
    const double k = 0.75 * sqrt(1.0 / (double)n);      /* ideal edge length */
    const double cutoff = 2.0 * k;                       /* repulsion range */
    int G = (int)(1.0 / cutoff);
    G = G < 1 ? 1 : G > 1024 ? 1024 : G;
    const size_t cells = (size_t)G * (size_t)G;

    double *dx   = malloc(n * sizeof(double));
    double *dy   = malloc(n * sizeof(double));
    int    *head = malloc(cells * sizeof(int));
    int    *next = malloc(n * sizeof(int));
    if (dx == NULL || dy == NULL || head == NULL || next == NULL) {
        free(dx); free(dy); free(head); free(next);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }

    for (int it = 0; it < iterations; ++it) {
        const double t = t0 * (1.0 - (double)it / (double)iterations) + 0.002 * t0;
        memset(dx, 0, n * sizeof(double));
        memset(dy, 0, n * sizeof(double));

        /* Bin vertices into the grid. */
        for (size_t q = 0; q < cells; ++q) { head[q] = -1; }
        const size_t shift = (size_t)it * 2654435761u % n;       /* rotate the insertion order */
        for (size_t q = 0; q < n; ++q) {
            const size_t i = (q + shift) % n;
            int cx = (int)(x[i] * G), cy = (int)(y[i] * G);
            cx = cx >= G ? G - 1 : cx < 0 ? 0 : cx;
            cy = cy >= G ? G - 1 : cy < 0 ? 0 : cy;
            const size_t cell = (size_t)cy * (size_t)G + (size_t)cx;
            next[i] = head[cell];
            head[cell] = (int)i;
        }
        /* Repulsion: only against vertices in the 3x3 surrounding cells.  A crowded cell
         * (early iterations of a dense core) would cost O(n^2), so each vertex looks at no
         * more than REPEL_CAP vertices per cell; the sample differs every iteration because
         * the lists are rebuilt in a different order. */
        for (size_t i = 0; i < n; ++i) {
            int cx = (int)(x[i] * G), cy = (int)(y[i] * G);
            cx = cx >= G ? G - 1 : cx < 0 ? 0 : cx;
            cy = cy >= G ? G - 1 : cy < 0 ? 0 : cy;
            for (int gy = cy - 1; gy <= cy + 1; ++gy) {
                for (int gx = cx - 1; gx <= cx + 1; ++gx) {
                    if (gx < 0 || gy < 0 || gx >= G || gy >= G) { continue; }
                    int seen = 0;
                    for (int j = head[(size_t)gy * (size_t)G + (size_t)gx]; j >= 0 && seen < REPEL_CAP;
                         j = next[j], ++seen) {
                        if ((size_t)j == i) { continue; }
                        double ex = x[i] - x[j], ey = y[i] - y[j];
                        double d2 = ex * ex + ey * ey;
                        if (d2 >= cutoff * cutoff) { continue; }
                        if (d2 < 1e-12) {                          /* identical spot: nudge apart */
                            ex = (rng_unit(rng) - 0.5) * 1e-4;
                            ey = (rng_unit(rng) - 0.5) * 1e-4;
                            d2 = ex * ex + ey * ey;
                        }
                        const double f = k * k / d2;               /* (k^2/d) * (1/d) */
                        dx[i] += ex * f;
                        dy[i] += ey * f;
                    }
                }
            }
        }
        /* Attraction along edges. */
        for (size_t i = 0; i < n; ++i) {
            for (size_t p = c->off[i]; p < c->off[i + 1]; ++p) {
                const size_t j = (size_t)c->adj[p];
                if (j == i) { continue; }
                const double ex = x[i] - x[j], ey = y[i] - y[j];
                const double d = sqrt(ex * ex + ey * ey);
                const double f = attract_scale * d / k;            /* (d^2/k) * (1/d) */
                dx[i] -= ex * f;  dy[i] -= ey * f;
                dx[j] += ex * f;  dy[j] += ey * f;
            }
        }
        /* Weak gravity keeps separate components from drifting apart, then move. */
        for (size_t i = 0; i < n; ++i) {
            dx[i] += (0.5 - x[i]) * 0.6 * k;
            dy[i] += (0.5 - y[i]) * 0.6 * k;
            const double len = sqrt(dx[i] * dx[i] + dy[i] * dy[i]);
            if (len > 0.0) {
                const double step = len < t ? len : t;
                x[i] = fmin(1.0, fmax(0.0, x[i] + dx[i] / len * step));
                y[i] = fmin(1.0, fmax(0.0, y[i] + dy[i] / len * step));
            }
        }
    }
    free(dx); free(dy); free(head); free(next);
    return GRAPH_OK;
}

/* ---- multilevel: coarsen, lay out the small graph, refine on the way back up ---- */

typedef struct {
    Csr  csr;
    int *to_coarse;       /* this level's vertex -> next level's vertex (NULL on the last) */
} Level;

static void csr_free(Csr *c)
{
    free(c->off);
    free(c->adj);
    memset(c, 0, sizeof(*c));
}

/** Merge pairs of adjacent vertices (random greedy matching) into a coarser graph. */
static bool coarsen(const Csr *fine, Rng *rng, Csr *coarse, int **to_coarse)
{
    const size_t n = fine->n;
    int *map = malloc((n + 1) * sizeof(int));
    int *order = malloc((n + 1) * sizeof(int));
    if (map == NULL || order == NULL) {
        free(map); free(order);
        return false;
    }
    for (size_t i = 0; i < n; ++i) { order[i] = (int)i; map[i] = -1; }
    for (size_t i = n; i > 1; --i) {                       /* Fisher-Yates shuffle */
        const size_t j = (size_t)(rng_unit(rng) * (double)i);
        const int t = order[i - 1]; order[i - 1] = order[j]; order[j] = t;
    }
    size_t nc = 0;
    for (size_t q = 0; q < n; ++q) {
        const int u = order[q];
        if (map[u] >= 0) { continue; }
        map[u] = (int)nc;
        for (size_t p = fine->off[u]; p < fine->off[u + 1]; ++p) {
            const int v = fine->adj[p];
            if (v != u && map[v] < 0) {                    /* first free neighbour joins u */
                map[v] = (int)nc;
                break;
            }
        }
        ++nc;
    }
    /* Build the coarse adjacency: every fine arc becomes a coarse arc (self-loops dropped). */
    size_t *off = calloc(nc + 2, sizeof(size_t));
    int    *adj = malloc((fine->off[n] + 1) * sizeof(int));
    if (off == NULL || adj == NULL) {
        free(map); free(order); free(off); free(adj);
        return false;
    }
    for (size_t u = 0; u < n; ++u) {
        for (size_t p = fine->off[u]; p < fine->off[u + 1]; ++p) {
            off[map[u] + 1] += map[fine->adj[p]] != map[u];
        }
    }
    for (size_t c = 0; c < nc; ++c) { off[c + 1] += off[c]; }
    size_t *fill = malloc((nc + 1) * sizeof(size_t));
    if (fill == NULL) {
        free(map); free(order); free(off); free(adj);
        return false;
    }
    memcpy(fill, off, nc * sizeof(size_t));
    for (size_t u = 0; u < n; ++u) {
        for (size_t p = fine->off[u]; p < fine->off[u + 1]; ++p) {
            if (map[fine->adj[p]] != map[u]) {
                adj[fill[map[u]]++] = map[fine->adj[p]];
            }
        }
    }
    free(fill);
    free(order);
    coarse->n   = nc;
    coarse->off = off;
    coarse->adj = adj;
    *to_coarse  = map;
    return true;
}

static int cmp_double(const void *a, const void *b)
{
    const double l = *(const double *)a, r = *(const double *)b;
    return (l > r) - (l < r);
}

#define COARSEST_SIZE 120
#define MAX_LEVELS    24

static GraphStatus layout_force(const Snap *s, int iterations, uint64_t seed, double *x, double *y)
{
    const size_t n = s->n;
    if (n == 0) {
        return GRAPH_OK;
    }
    if (n == 1) {
        x[0] = y[0] = 0.5;
        return GRAPH_OK;
    }
    Rng rng = { seed * 0x9E3779B97F4A7C15ULL + 0x1234567ULL };
    for (int i = 0; i < 4; ++i) { (void)rng_unit(&rng); }
    const double attract_scale = s->directed ? 1.0 : 0.5;     /* undirected arcs are stored twice */

    /* Level 0 borrows the snapshot's arrays; the others are built here. */
    Level lv[MAX_LEVELS];
    memset(lv, 0, sizeof(lv));
    lv[0].csr = (Csr){ n, s->off, s->adj };
    int levels = 1;
    while (levels < MAX_LEVELS && lv[levels - 1].csr.n > COARSEST_SIZE) {
        Csr next;
        int *map;
        if (!coarsen(&lv[levels - 1].csr, &rng, &next, &map)) { break; }
        if ((double)next.n > 0.85 * (double)lv[levels - 1].csr.n) {      /* stalled (e.g. a star) */
            csr_free(&next);
            free(map);
            break;
        }
        lv[levels - 1].to_coarse = map;
        lv[levels].csr = next;
        ++levels;
    }

    GraphStatus st = GRAPH_OK;
    double *cx = NULL, *cy = NULL;                          /* positions of the coarser level */
    for (int l = levels - 1; l >= 0 && st == GRAPH_OK; --l) {
        const Csr *c = &lv[l].csr;
        double *px = l == 0 ? x : malloc((c->n + 1) * sizeof(double));
        double *py = l == 0 ? y : malloc((c->n + 1) * sizeof(double));
        if (px == NULL || py == NULL) {
            if (l != 0) { free(px); free(py); }
            st = GRAPH_ERR_OUT_OF_MEMORY;
            break;
        }
        const double k = 0.75 * sqrt(1.0 / (double)c->n);
        int its;
        double t0;
        if (l == levels - 1) {                              /* coarsest: full layout */
            if (!initial_positions(c, k, &rng, px, py)) { st = GRAPH_ERR_OUT_OF_MEMORY; }
            its = levels == 1 ? iterations : (iterations > 300 ? iterations : 300);
            t0  = fmin(0.1, 4.0 * k);
        } else {                                            /* finer: inherit, jitter, polish */
            for (size_t i = 0; i < c->n; ++i) {
                const int p = lv[l].to_coarse[i];
                const double ang = 6.283185307179586 * rng_unit(&rng);
                px[i] = fmin(1.0, fmax(0.0, cx[p] + 0.5 * k * cos(ang)));
                py[i] = fmin(1.0, fmax(0.0, cy[p] + 0.5 * k * sin(ang)));
            }
            its = l == 0 ? (iterations / 3 > 30 ? iterations / 3 : 30) : 40;
            t0  = fmin(0.1, 1.5 * k);
        }
        if (st == GRAPH_OK) {
            st = fr_run(c, attract_scale, its, t0, &rng, px, py);
        }
        if (cx != NULL) { free(cx); free(cy); }
        cx = px;
        cy = py;
        if (l == 0) { cx = cy = NULL; }                     /* x, y belong to the caller */
    }
    free(cx);
    free(cy);
    for (int l = 0; l < levels; ++l) {
        free(lv[l].to_coarse);
        if (l > 0) { csr_free(&lv[l].csr); }
    }
    if (st != GRAPH_OK) {
        return st;
    }

    /* Fit the picture to [0,1].  For big graphs a few far-flung vertices would shrink the
     * dense part to a dot, so the 1st..99th percentile box is used and outliers are pulled
     * onto its border. */
    double lox = INFINITY, hix = -INFINITY, loy = INFINITY, hiy = -INFINITY;
    if (n >= 500) {
        double *tmp = malloc(n * sizeof(double));
        if (tmp != NULL) {
            memcpy(tmp, x, n * sizeof(double));
            qsort(tmp, n, sizeof(double), cmp_double);
            lox = tmp[n / 100]; hix = tmp[n - 1 - n / 100];
            memcpy(tmp, y, n * sizeof(double));
            qsort(tmp, n, sizeof(double), cmp_double);
            loy = tmp[n / 100]; hiy = tmp[n - 1 - n / 100];
            free(tmp);
        }
    }
    if (!(hix > lox) || !(hiy > loy)) {
        lox = loy = INFINITY;
        hix = hiy = -INFINITY;
        for (size_t i = 0; i < n; ++i) {
            lox = fmin(lox, x[i]); hix = fmax(hix, x[i]);
            loy = fmin(loy, y[i]); hiy = fmax(hiy, y[i]);
        }
    }
    const double span = fmax(fmax(hix - lox, hiy - loy), 1e-12);
    const double ox = 0.5 - (hix - lox) / (2.0 * span), oy = 0.5 - (hiy - loy) / (2.0 * span);
    for (size_t i = 0; i < n; ++i) {
        x[i] = fmin(1.0, fmax(0.0, ox + (x[i] - lox) / span));
        y[i] = fmin(1.0, fmax(0.0, oy + (y[i] - loy) / span));
    }
    return GRAPH_OK;
}

GraphStatus graph_layout(Graph *g, double *x, double *y, const VizOptions *opt)
{
    if (g == NULL || x == NULL || y == NULL) {
        return graph_set_error(GRAPH_ERR_INVALID_ARGUMENT);
    }
    VizOptions o = opt != NULL ? *opt : graph_viz_default_options();
    Snap s;
    GraphStatus st = snap_build(g, false, &s);
    if (st != GRAPH_OK) {
        return graph_set_error(st);
    }
    double *px = malloc((s.n + 1) * sizeof(double));
    double *py = malloc((s.n + 1) * sizeof(double));
    if (px == NULL || py == NULL) {
        free(px); free(py); snap_free(&s);
        return graph_set_error(GRAPH_ERR_OUT_OF_MEMORY);
    }
    st = layout_force(&s, o.iterations > 0 ? o.iterations : auto_iterations(s.n), o.seed, px, py);
    for (size_t v = 0; v < s.bound; ++v) {
        x[v] = y[v] = NAN;
    }
    for (size_t i = 0; i < s.n; ++i) {
        x[s.id[i]] = px[i];
        y[s.id[i]] = py[i];
    }
    free(px); free(py); snap_free(&s);
    return st == GRAPH_OK ? GRAPH_OK : graph_set_error(st);
}

/* ------------------------------------------------------------------ */
/*  Colours and the legend model                                        */
/* ------------------------------------------------------------------ */

typedef struct { unsigned char r, g, b; } Rgb;

static const Rgb PALETTE[12] = {
    { 31, 119, 180 }, { 255, 127, 14 }, { 44, 160, 44 }, { 214, 39, 40 },
    { 148, 103, 189 }, { 140, 86, 75 }, { 227, 119, 194 }, { 127, 127, 127 },
    { 188, 189, 34 }, { 23, 190, 207 }, { 174, 199, 232 }, { 255, 187, 120 },
};
static const Rgb GREY = { 150, 160, 180 };
static const Rgb INK  = { 25, 30, 45 };

/** Temperature scale: 0 = cold (deep blue) ... 1 = hot (dark red). */
static Rgb heat(double t)
{
    static const Rgb stops[8] = { { 49, 54, 149 }, { 69, 117, 180 }, { 116, 173, 209 }, { 254, 224, 144 },
                                  { 253, 174, 97 }, { 244, 109, 67 }, { 215, 48, 39 }, { 165, 0, 38 } };
    t = fmin(1.0, fmax(0.0, t)) * 7.0;
    const int i = t >= 7.0 ? 6 : (int)t;
    const double f = t - i;
    Rgb c;
    c.r = (unsigned char)(stops[i].r + f * (stops[i + 1].r - stops[i].r));
    c.g = (unsigned char)(stops[i].g + f * (stops[i + 1].g - stops[i].g));
    c.b = (unsigned char)(stops[i].b + f * (stops[i + 1].b - stops[i].b));
    return c;
}

/** Popularity 0..1 on a log scale, so the long tail of ordinary people is still visible. */
static double heat_t(size_t degree, size_t max_degree)
{
    return max_degree == 0 ? 0.0 : log1p((double)degree) / log1p((double)max_degree);
}

#define MAX_GROUPS 12

typedef struct {
    Rgb    color;
    size_t size;        /* people in this group                  */
    int    leader;      /* dense index of its most-connected member */
} Group;

typedef enum { LEG_NONE, LEG_GROUPS, LEG_HEAT } LegendKind;

typedef struct {
    LegendKind kind;
    Group      group[MAX_GROUPS];
    size_t     ngroups;
    size_t     other;           /* people in no coloured group (grey) */
    size_t     max_degree;
} Legend;

static size_t degree_of(const Snap *s, size_t i) { return s->off[i + 1] - s->off[i]; }

/** Fills col[i] and slot[i] (group slot or -1) and the legend model. */
static GraphStatus assign_colors(Graph *g, const Snap *s, VizColor scheme, Rgb *col, int *slot, Legend *lg)
{
    memset(lg, 0, sizeof(*lg));
    for (size_t i = 0; i < s->n; ++i) { col[i] = GREY; slot[i] = -1; }
    if (scheme == VIZ_COLOR_NONE || s->n == 0) {
        return GRAPH_OK;
    }
    int *label = malloc((s->bound + 1) * sizeof(int));
    if (label == NULL) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    bool categorical = false;
    if (scheme == VIZ_COLOR_AUTO || scheme == VIZ_COLOR_COMPONENT) {
        const size_t comps = graph_connected_components(g, label);
        categorical = comps > 1 || scheme == VIZ_COLOR_COMPONENT;
        if (!categorical) {
            scheme = VIZ_COLOR_COMMUNITY;                 /* one component: communities say more */
        }
    }
    if (scheme == VIZ_COLOR_COMMUNITY && s->n > 50000) {
        scheme = VIZ_COLOR_DEGREE;      /* Louvain would dominate the running time */
    }
    if (scheme == VIZ_COLOR_COMMUNITY) {
        double q;
        size_t count;
        if (graph_louvain(g, label, &q, &count) == GRAPH_OK) {
            categorical = true;
        } else {
            scheme = VIZ_COLOR_DEGREE;
        }
    }
    if (categorical) {
        /* Up to 12 groups get their own colour; with more, only the 11 largest do and the
         * rest share grey, so a hundred tiny communities do not turn into confetti. */
        size_t groups = 0;
        for (size_t i = 0; i < s->n; ++i) {
            groups = (size_t)label[s->id[i]] + 1 > groups ? (size_t)label[s->id[i]] + 1 : groups;
        }
        size_t *size = calloc(groups + 1, sizeof(size_t));
        int    *slot_of = malloc((groups + 1) * sizeof(int));
        if (size == NULL || slot_of == NULL) {
            free(size); free(slot_of); free(label);
            return GRAPH_ERR_OUT_OF_MEMORY;
        }
        size_t nonempty = 0;
        for (size_t i = 0; i < s->n; ++i) {
            nonempty += size[label[s->id[i]]]++ == 0;
        }
        for (size_t l = 0; l < groups; ++l) { slot_of[l] = -1; }
        const size_t colours = nonempty <= MAX_GROUPS ? nonempty : MAX_GROUPS - 1;
        for (size_t c = 0; c < colours; ++c) {          /* repeatedly take the largest unassigned */
            size_t best = groups;
            for (size_t l = 0; l < groups; ++l) {
                if (slot_of[l] < 0 && size[l] > 0 && (best == groups || size[l] > size[best])) {
                    best = l;
                }
            }
            slot_of[best] = (int)c;
            lg->group[c].color  = PALETTE[c];
            lg->group[c].size   = size[best];
            lg->group[c].leader = -1;
        }
        lg->ngroups = colours;
        for (size_t i = 0; i < s->n; ++i) {
            const int sl = slot_of[label[s->id[i]]];
            slot[i] = sl;
            col[i]  = sl >= 0 ? PALETTE[sl] : GREY;
            if (sl < 0) {
                ++lg->other;
            } else if (lg->group[sl].leader < 0 || degree_of(s, i) > degree_of(s, (size_t)lg->group[sl].leader)) {
                lg->group[sl].leader = (int)i;
            }
        }
        lg->kind = lg->ngroups >= 2 ? LEG_GROUPS : LEG_NONE;
        free(size);
        free(slot_of);
    } else if (scheme == VIZ_COLOR_DEGREE) {
        size_t maxd = 1;
        for (size_t i = 0; i < s->n; ++i) { maxd = degree_of(s, i) > maxd ? degree_of(s, i) : maxd; }
        for (size_t i = 0; i < s->n; ++i) { col[i] = heat(heat_t(degree_of(s, i), maxd)); }
        lg->kind       = LEG_HEAT;
        lg->max_degree = maxd;
    }
    free(label);
    return GRAPH_OK;
}

/* ------------------------------------------------------------------ */
/*  Geometry, drawing order, text tags                                  */
/* ------------------------------------------------------------------ */

static double clampd(double v, double lo, double hi) { return v < lo ? lo : v > hi ? hi : v; }

/** Text/line scale: 1 for pictures up to ~1500 px, 2 up to ~2500 px, ... */
static int text_scale(const VizOptions *o)
{
    const int m = o->width < o->height ? o->width : o->height;
    const int t = (m + 500) / 1000;
    return t < 1 ? 1 : t > 4 ? 4 : t;
}

typedef struct {
    double x, y, side;     /* the square the layout is drawn in        */
    double panel_x, panel_w;   /* legend panel (panel_w == 0: none)     */
    int    ts;
} Frame;

static Frame make_frame(const VizOptions *o, bool legend)
{
    Frame f;
    f.ts = text_scale(o);
    const double W = o->width, H = o->height, m = 0.035 * fmin(W, H);
    f.panel_w = legend ? fmin(0.38 * W, 290.0 * f.ts) : 0.0;
    f.panel_x = W - f.panel_w;
    const double title_h = o->title != NULL ? 24.0 * f.ts : 0.0;
    const double aw = W - f.panel_w - 2 * m, ah = H - 2 * m - title_h;
    f.side = fmax(fmin(aw, ah), 10.0);
    f.x = m + (aw - f.side) / 2;
    f.y = m + title_h + (ah - f.side) / 2;
    return f;
}

/** Vertex radius (pixels): grows with the square root of its degree. */
static double node_radius(const Snap *s, size_t i, double max_r, size_t maxdeg)
{
    return max_r * (0.4 + 0.6 * sqrt((double)degree_of(s, i) / (double)(maxdeg ? maxdeg : 1)));
}

static double max_node_radius(size_t n, double side)
{
    return clampd(side / (2.4 * sqrt((double)(n ? n : 1))), 1.5, 16.0);
}

static double edge_alpha(size_t drawn)
{
    return clampd(0.55 * sqrt(200.0 / (double)(drawn ? drawn : 1)), 0.01, 0.55);
}

static int cmp_u64(const void *a, const void *b)
{
    const uint64_t l = *(const uint64_t *)a, r = *(const uint64_t *)b;
    return (l > r) - (l < r);
}

/** Vertex indexes sorted by degree, lowest first, so the popular ones are painted last (on top). */
static int *order_by_degree(const Snap *s)
{
    uint64_t *key = malloc((s->n + 1) * sizeof(uint64_t));
    int *order = malloc((s->n + 1) * sizeof(int));
    if (key == NULL || order == NULL) {
        free(key); free(order);
        return NULL;
    }
    for (size_t i = 0; i < s->n; ++i) { key[i] = ((uint64_t)degree_of(s, i) << 32) | (uint64_t)i; }
    qsort(key, s->n, sizeof(uint64_t), cmp_u64);
    for (size_t i = 0; i < s->n; ++i) { order[i] = (int)(key[i] & 0xFFFFFFFFu); }
    free(key);
    return order;
}

static void vlabel(const VizOptions *o, int id, char *buf, size_t n)
{
    if (o->labels != NULL && o->labels[id] != NULL) {
        snprintf(buf, n, "%s", o->labels[id]);
    } else {
        snprintf(buf, n, "%d", id);
    }
}

static void fmt_count(char *buf, size_t n, size_t v)
{
    char raw[32];
    snprintf(raw, sizeof raw, "%zu", v);
    const size_t len = strlen(raw);
    size_t o = 0;
    for (size_t i = 0; i < len && o + 2 < n; ++i) {
        if (i > 0 && (len - i) % 3 == 0) { buf[o++] = ','; }
        buf[o++] = raw[i];
    }
    buf[o] = '\0';
}

/** A piece of text to print on the picture: a community name or a hub's name. */
typedef struct {
    double x, y;          /* top-left of the text box, pixels */
    char   text[44];
    Rgb    color;
    bool   group;         /* a community tag (drawn larger)   */
} Tag;

typedef struct { double x0, y0, x1, y1; } Box;

static bool overlaps(const Box *a, const Box *b)
{
    return a->x0 < b->x1 && b->x0 < a->x1 && a->y0 < b->y1 && b->y0 < a->y1;
}

static Rgb darker(Rgb c)
{
    return (Rgb){ (unsigned char)(c.r * 0.55), (unsigned char)(c.g * 0.55), (unsigned char)(c.b * 0.55) };
}

/**
 * Chooses the texts to write on the picture: one tag per community at the centre of its
 * members, then the names of the most-connected people.  Tags that would overlap an
 * earlier one are dropped.  Returns the number of tags written to @out (capacity @cap).
 */
static size_t build_tags(const Snap *s, const double *px, const double *py, const int *slot, const Legend *lg,
                         const VizOptions *o, const Frame *fr, double max_r, size_t maxdeg,
                         const int *order, int hubs, Tag *out, size_t cap)
{
    size_t count = 0;
    Box boxes[64];
    const double cw = (double)(VIZ_FONT_W * fr->ts), ch = (double)(VIZ_FONT_H * fr->ts);
    const double gw = (double)(VIZ_FONT_W * (fr->ts + 1)), gh = (double)(VIZ_FONT_H * (fr->ts + 1));
    char name[48];
    if (lg->kind == LEG_GROUPS && o->legend) {
        for (size_t c = 0; c < lg->ngroups && count < cap && count < 64; ++c) {
            double sx = 0, sy = 0;
            size_t m = 0;
            for (size_t i = 0; i < s->n; ++i) {
                if (slot[i] == (int)c) { sx += px[i]; sy += py[i]; ++m; }
            }
            if (m == 0) { continue; }
            Tag *t = &out[count];
            snprintf(t->text, sizeof t->text, "%s %zu", o->group_name != NULL ? o->group_name : "Group", c + 1);
            const double w = (double)strlen(t->text) * gw;
            t->x = fr->x + sx / (double)m * fr->side - w / 2;
            t->y = fr->y + sy / (double)m * fr->side - gh / 2;
            t->color = darker(lg->group[c].color);
            t->group = true;
            /* Crowded middle: try the spots just above and below until the tag overlaps nothing. */
            static const double nudge[7] = { 0, 1, -1, 2, -2, 3, -3 };
            const double base_y = t->y;
            bool placed = false;
            for (int k = 0; k < 7 && !placed; ++k) {
                t->y = base_y + nudge[k] * gh * 1.05;
                const Box b = { t->x, t->y, t->x + w, t->y + gh };
                bool clash = false;
                for (size_t j = 0; j < count && !clash; ++j) { clash = overlaps(&b, &boxes[j]); }
                if (!clash) { boxes[count] = b; placed = true; }
            }
            if (placed) { ++count; }               /* otherwise the legend still names the group */
        }
    }
    const size_t first_hub = count;
    for (int k = 0; k < hubs && (size_t)k < s->n && count < cap && count < 64; ++k) {
        const int i = order[s->n - 1 - (size_t)k];
        vlabel(o, s->id[i], name, sizeof name);
        Tag *t = &out[count];
        snprintf(t->text, sizeof t->text, "%.20s", name);
        const double w = (double)strlen(t->text) * cw;
        const double cx = fr->x + px[i] * fr->side, cy = fr->y + py[i] * fr->side;
        const double r = node_radius(s, (size_t)i, max_r, maxdeg);
        t->x = cx + r + 3;
        if (t->x + w > fr->x + fr->side + 0.04 * o->width) { t->x = cx - r - 3 - w; }
        t->y = cy - ch / 2;
        t->color = INK;
        t->group = false;
        const Box b = { t->x, t->y, t->x + w, t->y + ch };
        bool clash = false;
        for (size_t j = 0; j < count && !clash; ++j) { clash = overlaps(&b, &boxes[j]); }
        if (clash) { continue; }
        boxes[count] = b;
        ++count;
    }
    (void)first_hub;
    return count;
}

/** One line of the legend: a colour swatch and up to two lines of text. */
typedef struct {
    Rgb  color;
    char line1[64];
    char line2[64];
} LegendRow;

static size_t legend_rows(const Snap *s, const Legend *lg, const VizOptions *o, const int *order,
                          LegendRow *rows, size_t cap)
{
    size_t n = 0;
    char count[32], name[48];
    if (lg->kind == LEG_GROUPS) {
        for (size_t c = 0; c < lg->ngroups && n < cap; ++c) {
            fmt_count(count, sizeof count, lg->group[c].size);
            rows[n].color = lg->group[c].color;
            snprintf(rows[n].line1, sizeof rows[n].line1, "%s %zu: %s people",
                     o->group_name != NULL ? o->group_name : "Group", c + 1, count);
            const int ld = lg->group[c].leader;
            if (ld >= 0) {
                vlabel(o, s->id[ld], name, sizeof name);
                snprintf(rows[n].line2, sizeof rows[n].line2, "hub: %.14s (%zu)", name, degree_of(s, (size_t)ld));
            } else {
                rows[n].line2[0] = '\0';
            }
            ++n;
        }
        if (lg->other > 0 && n < cap) {
            fmt_count(count, sizeof count, lg->other);
            rows[n].color = GREY;
            snprintf(rows[n].line1, sizeof rows[n].line1, "Other: %s people", count);
            snprintf(rows[n].line2, sizeof rows[n].line2, "smaller groups");
            ++n;
        }
    } else if (lg->kind == LEG_HEAT) {
        for (int k = 0; k < 8 && (size_t)k < s->n && n < cap; ++k) {
            const int i = order[s->n - 1 - (size_t)k];
            vlabel(o, s->id[i], name, sizeof name);
            rows[n].color = heat(heat_t(degree_of(s, (size_t)i), lg->max_degree));
            snprintf(rows[n].line1, sizeof rows[n].line1, "%d. %.16s  %zu", k + 1, name, degree_of(s, (size_t)i));
            rows[n].line2[0] = '\0';
            ++n;
        }
    }
    return n;
}

/** The degrees worth a tick on the colour bar: 1, 2, 5, 10, 20, 50, ... up to the maximum. */
static size_t heat_ticks(size_t max_degree, size_t *out, size_t cap)
{
    size_t n = 0;
    for (size_t base = 1; base <= max_degree && n < cap; base *= 10) {
        const size_t m[3] = { 1, 2, 5 };
        for (int k = 0; k < 3; ++k) {
            if (base * m[k] <= max_degree && n < cap) { out[n++] = base * m[k]; }
        }
    }
    return n;
}

/* ------------------------------------------------------------------ */
/*  SVG                                                                 */
/* ------------------------------------------------------------------ */

static void svg_text(FILE *f, const char *t)
{
    for (; *t; ++t) {
        switch (*t) {
        case '&': fputs("&amp;", f); break;
        case '<': fputs("&lt;", f); break;
        case '>': fputs("&gt;", f); break;
        case '"': fputs("&quot;", f); break;
        default:  fputc(*t, f); break;
        }
    }
}

static void svg_label(FILE *f, double x, double y, int size, Rgb c, const char *text, bool bold)
{
    fprintf(f, "<text x=\"%.1f\" y=\"%.1f\" font-size=\"%d\" fill=\"#%02x%02x%02x\" stroke=\"#ffffff\" "
               "stroke-width=\"%d\" paint-order=\"stroke\"%s>", x, y, size, c.r, c.g, c.b, size > 16 ? 4 : 3,
            bold ? " font-weight=\"bold\"" : "");
    svg_text(f, text);
    fputs("</text>\n", f);
}

static const char *plural_name(const VizOptions *o, char *buf, size_t n)
{
    const char *g = o->group_name != NULL ? o->group_name : "Group";
    const size_t len = strlen(g);
    if (len > 0 && g[len - 1] == 'y') { snprintf(buf, n, "%.*sies", (int)len - 1, g); }
    else { snprintf(buf, n, "%ss", g); }
    return buf;
}

static void svg_legend(FILE *f, const Snap *s, const Legend *lg, const VizOptions *o, const Frame *fr, const int *order)
{
    const int ts = fr->ts;
    const double x0 = fr->panel_x + 12.0 * ts;
    double y = fr->y + 4.0 * ts;
    char plural[40], buf[64];
    LegendRow rows[24];
    const size_t nrows = legend_rows(s, lg, o, order, rows, 24);
    fputs("<g font-family=\"monospace\">\n", f);
    if (lg->kind == LEG_GROUPS) {
        fprintf(f, "<text x=\"%.1f\" y=\"%.1f\" font-size=\"%d\" font-weight=\"bold\" fill=\"#222\">", x0, y + 12.0 * ts, 14 * ts);
        svg_text(f, plural_name(o, plural, sizeof plural));
        fputs(" (colour)</text>\n", f);
        y += 26.0 * ts;
    } else {
        fprintf(f, "<defs><linearGradient id=\"heat\" x1=\"0\" y1=\"1\" x2=\"0\" y2=\"0\">");
        for (int k = 0; k <= 7; ++k) {
            const Rgb c = heat(k / 7.0);
            fprintf(f, "<stop offset=\"%.3f\" stop-color=\"#%02x%02x%02x\"/>", k / 7.0, c.r, c.g, c.b);
        }
        fputs("</linearGradient></defs>\n", f);
        fprintf(f, "<text x=\"%.1f\" y=\"%.1f\" font-size=\"%d\" font-weight=\"bold\" fill=\"#222\">Popularity</text>\n",
                x0, y + 12.0 * ts, 14 * ts);
        y += 24.0 * ts;
        const double bh = 230.0 * ts, bw = 22.0 * ts;
        fprintf(f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"url(#heat)\" stroke=\"#444\"/>\n", x0, y, bw, bh);
        size_t ticks[16];
        const size_t nt = heat_ticks(lg->max_degree, ticks, 16);
        for (size_t k = 0; k < nt; ++k) {
            const double ty = y + bh - heat_t(ticks[k], lg->max_degree) * bh;
            fprintf(f, "<path d=\"M%.1f %.1f h%.1f\" stroke=\"#222\"/>\n", x0 + bw, ty, 5.0 * ts);
            snprintf(buf, sizeof buf, "%zu", ticks[k]);
            fprintf(f, "<text x=\"%.1f\" y=\"%.1f\" font-size=\"%d\" fill=\"#222\">%s</text>\n", x0 + bw + 8.0 * ts, ty + 4.0 * ts, 11 * ts, buf);
        }
        fprintf(f, "<text x=\"%.1f\" y=\"%.1f\" font-size=\"%d\" fill=\"#444\">friends (log scale)</text>\n", x0, y + bh + 16.0 * ts, 11 * ts);
        fprintf(f, "<text x=\"%.1f\" y=\"%.1f\" font-size=\"%d\" fill=\"#444\">hot = very popular</text>\n", x0, y + bh + 30.0 * ts, 11 * ts);
        y += bh + 50.0 * ts;
        fprintf(f, "<text x=\"%.1f\" y=\"%.1f\" font-size=\"%d\" font-weight=\"bold\" fill=\"#222\">Most popular</text>\n", x0, y + 4.0 * ts, 13 * ts);
        y += 18.0 * ts;
    }
    for (size_t r = 0; r < nrows; ++r) {
        const double rh = (rows[r].line2[0] ? 36.0 : 20.0) * ts;
        if (y + rh > o->height) { break; }
        fprintf(f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" rx=\"3\" fill=\"#%02x%02x%02x\" stroke=\"#333\"/>\n",
                x0, y + 2.0 * ts, 14.0 * ts, 14.0 * ts, rows[r].color.r, rows[r].color.g, rows[r].color.b);
        fprintf(f, "<text x=\"%.1f\" y=\"%.1f\" font-size=\"%d\" fill=\"#111\">", x0 + 22.0 * ts, y + 14.0 * ts, 12 * ts);
        svg_text(f, rows[r].line1);
        fputs("</text>\n", f);
        if (rows[r].line2[0]) {
            fprintf(f, "<text x=\"%.1f\" y=\"%.1f\" font-size=\"%d\" fill=\"#555\">", x0 + 22.0 * ts, y + 30.0 * ts, 11 * ts);
            svg_text(f, rows[r].line2);
            fputs("</text>\n", f);
        }
        y += rh;
    }
    fputs("</g>\n", f);
}

static GraphStatus render_svg(FILE *f, const Snap *s, const double *px, const double *py, const Rgb *col,
                              const int *slot, const Legend *lg, const int *order, const VizOptions *o, size_t drawn)
{
    const bool legend = o->legend && lg->kind != LEG_NONE;
    const Frame fr = make_frame(o, legend);
    const int ts = fr.ts;
    size_t maxdeg = 1;
    for (size_t i = 0; i < s->n; ++i) { maxdeg = degree_of(s, i) > maxdeg ? degree_of(s, i) : maxdeg; }
    const double max_r = max_node_radius(s->n, fr.side);
    const bool labels  = o->show_labels == 1 || (o->show_labels < 0 && s->n <= AUTO_LABELS);
    const bool weights = o->show_weights && s->w != NULL && drawn <= AUTO_LABELS;
    const double alpha = edge_alpha(drawn);
#define PX(i) (fr.x + px[i] * fr.side)
#define PY(i) (fr.y + py[i] * fr.side)

    fprintf(f, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%d\" height=\"%d\" viewBox=\"0 0 %d %d\" "
               "font-family=\"monospace\">\n", o->width, o->height, o->width, o->height);
    fputs("<rect width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>\n", f);
    if (s->directed) {
        fputs("<defs><marker id=\"arrow\" viewBox=\"0 0 10 10\" refX=\"10\" refY=\"5\" "
              "markerWidth=\"6\" markerHeight=\"6\" orient=\"auto\">"
              "<path d=\"M0,0 L10,5 L0,10 z\" fill=\"#445\"/></marker></defs>\n", f);
    }
    if (o->title != NULL) {
        fprintf(f, "<text x=\"%.1f\" y=\"%.1f\" font-size=\"%d\" font-weight=\"bold\" fill=\"#222\">", 0.035 * fmin(o->width, o->height),
                0.035 * fmin(o->width, o->height) + 16.0 * ts, 17 * ts);
        svg_text(f, o->title);
        fputs("</text>\n", f);
    }
    fprintf(f, "<g stroke=\"#445566\" stroke-opacity=\"%.3f\" stroke-width=\"%.2f\">\n", alpha, drawn > 5000 ? 0.6 : 1.2);
    for (size_t i = 0; i < s->n; ++i) {
        for (size_t p = s->off[i]; p < s->off[i + 1]; ++p) {
            const int j = s->adj[p];
            if (!draw_arc(s, i, j)) { continue; }
            double x1 = PX(i), y1 = PY(i), x2 = PX(j), y2 = PY(j);
            if (s->directed) {          /* stop at the target's rim so the arrowhead is visible */
                const double dx = x2 - x1, dy = y2 - y1, len = sqrt(dx * dx + dy * dy);
                const double r = node_radius(s, (size_t)j, max_r, maxdeg);
                if (len > r + 1.0) { x2 -= dx / len * r; y2 -= dy / len * r; }
            }
            fprintf(f, "<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\"%s/>\n", x1, y1, x2, y2,
                    s->directed ? " marker-end=\"url(#arrow)\"" : "");
        }
    }
    fputs("</g>\n", f);
    if (weights) {
        fputs("<g font-size=\"10\" fill=\"#a33\">\n", f);
        for (size_t i = 0; i < s->n; ++i) {
            for (size_t p = s->off[i]; p < s->off[i + 1]; ++p) {
                if (!draw_arc(s, i, s->adj[p])) { continue; }
                const size_t j = (size_t)s->adj[p];
                fprintf(f, "<text x=\"%.1f\" y=\"%.1f\">%g</text>\n", (PX(i) + PX(j)) * 0.5, (PY(i) + PY(j)) * 0.5, s->w[p]);
            }
        }
        fputs("</g>\n", f);
    }
    /* Nodes, least connected first so the popular ones end up on top.  In temperature mode the
     * popular ones also get a soft halo. */
    fputs("<g stroke=\"#222\" stroke-width=\"0.8\">\n", f);
    for (size_t q = 0; q < s->n; ++q) {
        const size_t i = (size_t)order[q];
        const double r = node_radius(s, i, max_r, maxdeg);
        if (lg->kind == LEG_HEAT && heat_t(degree_of(s, i), maxdeg) > 0.62) {
            fprintf(f, "<circle cx=\"%.1f\" cy=\"%.1f\" r=\"%.1f\" fill=\"#%02x%02x%02x\" fill-opacity=\"0.25\" stroke=\"none\"/>\n",
                    PX(i), PY(i), r * 2.6, col[i].r, col[i].g, col[i].b);
        }
        fprintf(f, "<circle cx=\"%.1f\" cy=\"%.1f\" r=\"%.1f\" fill=\"#%02x%02x%02x\"/>\n", PX(i), PY(i), r, col[i].r, col[i].g, col[i].b);
    }
    fputs("</g>\n", f);
    if (labels) {
        fprintf(f, "<g font-size=\"%.0f\" fill=\"#111\" text-anchor=\"middle\">\n", clampd(max_r * 0.9, 8.0, 14.0));
        char name[48];
        for (size_t i = 0; i < s->n; ++i) {
            fprintf(f, "<text x=\"%.1f\" y=\"%.1f\">", PX(i), PY(i) - node_radius(s, i, max_r, maxdeg) - 3.0);
            vlabel(o, s->id[i], name, sizeof name);
            svg_text(f, name);
            fputs("</text>\n", f);
        }
        fputs("</g>\n", f);
    }
    {
        const int hubs = o->label_top >= 0 ? o->label_top : (labels || o->show_labels == 0 ? 0 : 12);
        Tag tags[64];
        const size_t nt = build_tags(s, px, py, slot, lg, o, &fr, max_r, maxdeg, order, hubs, tags, 64);
        for (size_t t = 0; t < nt; ++t) {
            if (tags[t].group) {
                svg_label(f, tags[t].x, tags[t].y + (double)(VIZ_FONT_H * (ts + 1)) * 0.8, 13 * (ts + 1), tags[t].color, tags[t].text, true);
            } else {
                svg_label(f, tags[t].x, tags[t].y + (double)(VIZ_FONT_H * ts) * 0.8, 13 * ts, tags[t].color, tags[t].text, false);
            }
        }
    }
    if (legend) { svg_legend(f, s, lg, o, &fr, order); }
    fputs("</svg>\n", f);
#undef PX
#undef PY
    return ferror(f) ? GRAPH_ERR_IO : GRAPH_OK;
}

/* ------------------------------------------------------------------ */
/*  Raster                                                              */
/* ------------------------------------------------------------------ */

typedef struct {
    int    w, h;
    float *px;      /* w * h * 3, linear 0..1 */
} Canvas;

static void blend(Canvas *c, int x, int y, Rgb col, float a)
{
    if (x < 0 || y < 0 || x >= c->w || y >= c->h) {
        return;
    }
    float *p = c->px + ((size_t)y * (size_t)c->w + (size_t)x) * 3;
    p[0] += (col.r / 255.0f - p[0]) * a;
    p[1] += (col.g / 255.0f - p[1]) * a;
    p[2] += (col.b / 255.0f - p[2]) * a;
}

static void draw_line(Canvas *c, double fx0, double fy0, double fx1, double fy1, Rgb col, float a)
{
    int x0 = (int)lround(fx0), y0 = (int)lround(fy0), x1 = (int)lround(fx1), y1 = (int)lround(fy1);
    const int dx = abs(x1 - x0), dy = -abs(y1 - y0);
    const int sx = x0 < x1 ? 1 : -1, sy = y0 < y1 ? 1 : -1;
    int err = dx + dy;
    for (;;) {
        blend(c, x0, y0, col, a);
        if (x0 == x1 && y0 == y1) { break; }
        const int e2 = 2 * err;
        if (e2 >= dy) { err += dy; x0 += sx; }
        if (e2 <= dx) { err += dx; y0 += sy; }
    }
}

static void fill_circle(Canvas *c, double cx, double cy, double r, Rgb fill, bool outline)
{
    const int x0 = (int)floor(cx - r - 1), x1 = (int)ceil(cx + r + 1);
    const int y0 = (int)floor(cy - r - 1), y1 = (int)ceil(cy + r + 1);
    const Rgb dark = { 30, 30, 40 };
    for (int y = y0; y <= y1; ++y) {
        for (int x = x0; x <= x1; ++x) {
            const double d = sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy));
            if (d <= r - (outline ? 1.0 : 0.0)) {
                blend(c, x, y, fill, 1.0f);
            } else if (outline && d <= r) {
                blend(c, x, y, dark, 1.0f);
            } else if (d <= r + 1.0) {                      /* soft edge */
                blend(c, x, y, outline ? dark : fill, (float)(r + 1.0 - d) * 0.5f);
            }
        }
    }
}

static unsigned char to_byte(float v)
{
    return (unsigned char)(clampd(v, 0.0, 1.0) * 255.0 + 0.5);
}

static GraphStatus canvas_png(FILE *f, const Canvas *c)
{
    const size_t n = (size_t)c->w * (size_t)c->h * 3;
    unsigned char *bytes = malloc(n);
    if (bytes == NULL) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    for (size_t i = 0; i < n; ++i) {
        bytes[i] = to_byte(c->px[i]);
    }
    const GraphStatus st = png_write_rgb(f, bytes, (unsigned)c->w, (unsigned)c->h);
    free(bytes);
    return st;
}

static bool canvas_init(Canvas *c, int w, int h)
{
    if ((unsigned long long)w * (unsigned long long)h > MAX_RASTER_PIXELS) {
        return false;
    }
    c->w  = w;
    c->h  = h;
    c->px = malloc((size_t)w * (size_t)h * 3 * sizeof(float));
    if (c->px == NULL) {
        return false;
    }
    for (size_t i = 0; i < (size_t)w * (size_t)h * 3; ++i) { c->px[i] = 1.0f; }
    return true;
}

/* ---- text, using the built-in 8x16 bitmap font ---- */

static int text_pixels(const char *t, int scale)
{
    return (int)strlen(t) * VIZ_FONT_W * scale;
}

static void put_glyph(Canvas *c, int x, int y, char ch, Rgb col, int scale, float a)
{
    if (ch < 32 || ch > 126) { ch = '?'; }
    const unsigned char *g = VIZ_FONT[ch - 32];
    for (int row = 0; row < VIZ_FONT_H; ++row) {
        for (int bit = 0; bit < VIZ_FONT_W; ++bit) {
            if (!((g[row] >> (VIZ_FONT_W - 1 - bit)) & 1u)) { continue; }
            for (int dy = 0; dy < scale; ++dy) {
                for (int dx = 0; dx < scale; ++dx) {
                    blend(c, x + bit * scale + dx, y + row * scale + dy, col, a);
                }
            }
        }
    }
}

/** Draws text with its top-left corner at (x, y); @halo outlines it in white so it reads on any background. */
static void draw_text(Canvas *c, int x, int y, const char *t, Rgb col, int scale, bool halo)
{
    const Rgb white = { 255, 255, 255 };
    if (halo) {
        for (int k = 0; t[k]; ++k) {
            for (int oy = -1; oy <= 1; ++oy) {
                for (int ox = -1; ox <= 1; ++ox) {
                    if (ox || oy) { put_glyph(c, x + k * VIZ_FONT_W * scale + ox * scale, y + oy * scale, t[k], white, scale, 0.9f); }
                }
            }
        }
    }
    for (int k = 0; k < (int)strlen(t); ++k) { put_glyph(c, x + k * VIZ_FONT_W * scale, y, t[k], col, scale, 1.0f); }
}

static void fill_rect(Canvas *c, int x0, int y0, int x1, int y1, Rgb col)
{
    for (int y = y0; y < y1; ++y) { for (int x = x0; x < x1; ++x) { blend(c, x, y, col, 1.0f); } }
}

static void frame_rect(Canvas *c, int x0, int y0, int x1, int y1, Rgb col)
{
    for (int x = x0; x < x1; ++x) { blend(c, x, y0, col, 1.0f); blend(c, x, y1 - 1, col, 1.0f); }
    for (int y = y0; y < y1; ++y) { blend(c, x0, y, col, 1.0f); blend(c, x1 - 1, y, col, 1.0f); }
}

/** A soft round glow (used behind the most popular vertices in temperature mode). */
static void glow(Canvas *c, double cx, double cy, double r, Rgb col, float a)
{
    for (int y = (int)floor(cy - r); y <= (int)ceil(cy + r); ++y) {
        for (int x = (int)floor(cx - r); x <= (int)ceil(cx + r); ++x) {
            const double d = sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy)) / r;
            if (d < 1.0) { blend(c, x, y, col, a * (float)((1.0 - d) * (1.0 - d))); }
        }
    }
}

static void png_legend(Canvas *c, const Snap *s, const Legend *lg, const VizOptions *o, const Frame *fr, const int *order)
{
    const int ts = fr->ts;
    const int x0 = (int)(fr->panel_x + 12.0 * ts);
    int y = (int)(fr->y + 4.0 * ts);
    char plural[40], buf[64];
    LegendRow rows[24];
    const size_t nrows = legend_rows(s, lg, o, order, rows, 24);
    const Rgb line = { 70, 70, 80 }, faint = { 90, 95, 110 };
    if (lg->kind == LEG_GROUPS) {
        snprintf(buf, sizeof buf, "%s (colour)", plural_name(o, plural, sizeof plural));
        draw_text(c, x0, y, buf, INK, ts, false);
        y += (VIZ_FONT_H + 10) * ts;
    } else {
        draw_text(c, x0, y, "Popularity", INK, ts, false);
        y += (VIZ_FONT_H + 6) * ts;
        const int bw = 22 * ts, bh = 230 * ts;
        for (int k = 0; k < bh; ++k) {
            fill_rect(c, x0, y + k, x0 + bw, y + k + 1, heat(1.0 - (double)k / (double)(bh - 1)));
        }
        frame_rect(c, x0, y, x0 + bw, y + bh, line);
        size_t ticks[16];
        const size_t nt = heat_ticks(lg->max_degree, ticks, 16);
        int last_y = -1000;
        for (size_t k = 0; k < nt; ++k) {
            const int ty = y + bh - 1 - (int)(heat_t(ticks[k], lg->max_degree) * (double)(bh - 1));
            fill_rect(c, x0 + bw, ty, x0 + bw + 5 * ts, ty + 1, line);
            if (abs(ty - last_y) >= VIZ_FONT_H * ts - 2 * ts || last_y < -500) {
                snprintf(buf, sizeof buf, "%zu", ticks[k]);
                draw_text(c, x0 + bw + 8 * ts, ty - VIZ_FONT_H * ts / 2, buf, INK, ts, false);
                last_y = ty;
            }
        }
        draw_text(c, x0, y + bh + 6 * ts, "friends (log scale)", faint, ts, false);
        draw_text(c, x0, y + bh + (6 + VIZ_FONT_H) * ts, "hot = very popular", faint, ts, false);
        y += bh + (16 + 2 * VIZ_FONT_H) * ts;
        draw_text(c, x0, y, "Most popular", INK, ts, false);
        y += (VIZ_FONT_H + 6) * ts;
    }
    for (size_t r = 0; r < nrows; ++r) {
        const int rh = ((rows[r].line2[0] ? 2 : 1) * VIZ_FONT_H + 6) * ts;
        if (y + rh > o->height) { break; }
        fill_rect(c, x0, y + 2 * ts, x0 + 14 * ts, y + 16 * ts, rows[r].color);
        frame_rect(c, x0, y + 2 * ts, x0 + 14 * ts, y + 16 * ts, line);
        draw_text(c, x0 + 22 * ts, y, rows[r].line1, INK, ts, false);
        if (rows[r].line2[0]) { draw_text(c, x0 + 22 * ts, y + VIZ_FONT_H * ts, rows[r].line2, faint, ts, false); }
        y += rh;
    }
}

static GraphStatus render_png(FILE *f, const Snap *s, const double *px, const double *py, const Rgb *col,
                              const int *slot, const Legend *lg, const int *order, const VizOptions *o, size_t drawn)
{
    Canvas c;
    if (!canvas_init(&c, o->width, o->height)) {
        return (unsigned long long)o->width * (unsigned long long)o->height > MAX_RASTER_PIXELS
                   ? GRAPH_ERR_TOO_LARGE : GRAPH_ERR_OUT_OF_MEMORY;
    }
    const bool legend = o->legend && lg->kind != LEG_NONE;
    const Frame fr = make_frame(o, legend);
    const int ts = fr.ts;
    size_t maxdeg = 1;
    for (size_t i = 0; i < s->n; ++i) { maxdeg = degree_of(s, i) > maxdeg ? degree_of(s, i) : maxdeg; }
    const double max_r = max_node_radius(s->n, fr.side) * (lg->kind == LEG_HEAT ? 1.0 : 1.0);
    const float  alpha = (float)edge_alpha(drawn);
    const Rgb    ink   = { 50, 60, 85 };
    const bool   arrows = s->directed && s->n <= 300;
#define PX(i) (fr.x + px[i] * fr.side)
#define PY(i) (fr.y + py[i] * fr.side)

    for (size_t i = 0; i < s->n; ++i) {
        for (size_t p = s->off[i]; p < s->off[i + 1]; ++p) {
            const int j = s->adj[p];
            if (!draw_arc(s, i, j)) { continue; }
            const double x1 = PX(i), y1 = PY(i), x2 = PX(j), y2 = PY(j);
            draw_line(&c, x1, y1, x2, y2, ink, alpha);
            if (arrows) {
                const double dx = x2 - x1, dy = y2 - y1, len = sqrt(dx * dx + dy * dy);
                if (len > 12.0) {
                    const double r  = node_radius(s, (size_t)j, max_r, maxdeg) + 1.0;
                    const double ux = dx / len, uy = dy / len;
                    const double tx = x2 - ux * r, ty = y2 - uy * r;
                    draw_line(&c, tx, ty, tx - ux * 9 + uy * 4, ty - uy * 9 - ux * 4, ink, 0.9f);
                    draw_line(&c, tx, ty, tx - ux * 9 - uy * 4, ty - uy * 9 + ux * 4, ink, 0.9f);
                }
            }
        }
    }
    /* Least connected first: the popular vertices are painted last, on top, with a glow in temperature mode. */
    for (size_t q = 0; q < s->n; ++q) {
        const size_t i = (size_t)order[q];
        const double r = node_radius(s, i, max_r, maxdeg);
        if (lg->kind == LEG_HEAT && heat_t(degree_of(s, i), maxdeg) > 0.62) {
            glow(&c, PX(i), PY(i), r * 3.2, col[i], 0.55f);
        }
        fill_circle(&c, PX(i), PY(i), r, col[i], r >= 3.0);
    }
    /* Texts: community names at the centre of each group, then the most popular people. */
    {
        const bool all_labels = o->show_labels == 1 && s->n <= 400;
        const int hubs = all_labels ? (int)s->n : (o->label_top >= 0 ? o->label_top : (o->show_labels == 0 ? 0 : 12));
        Tag tags[400];
        const size_t nt = build_tags(s, px, py, slot, lg, o, &fr, max_r, maxdeg, order, hubs, tags, 400);
        for (size_t t = 0; t < nt; ++t) {
            draw_text(&c, (int)tags[t].x, (int)tags[t].y, tags[t].text, tags[t].color, tags[t].group ? ts + 1 : ts, true);
        }
    }
    if (o->title != NULL) {
        draw_text(&c, (int)(0.035 * fmin(o->width, o->height)), (int)(0.035 * fmin(o->width, o->height)), o->title, INK, ts + 1 > 2 ? ts : ts + 1, false);
    }
    if (legend) { png_legend(&c, s, lg, o, &fr, order); }
    (void)text_pixels;
#undef PX
#undef PY
    const GraphStatus st = canvas_png(f, &c);
    free(c.px);
    return st;
}

/** Adjacency-matrix heat-map: pixel (col, row) darkens with the number of edges it covers. */
static GraphStatus render_matrix_png(FILE *f, const Snap *s, const VizOptions *o)
{
    Canvas c;
    if (!canvas_init(&c, o->width, o->height)) {
        return (unsigned long long)o->width * (unsigned long long)o->height > MAX_RASTER_PIXELS
                   ? GRAPH_ERR_TOO_LARGE : GRAPH_ERR_OUT_OF_MEMORY;
    }
    const size_t pixels = (size_t)c.w * (size_t)c.h;
    unsigned *count = calloc(pixels, sizeof(unsigned));
    if (count == NULL) {
        free(c.px);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    unsigned top = 0;
    const double sx = (double)c.w / (double)(s->n ? s->n : 1), sy = (double)c.h / (double)(s->n ? s->n : 1);
    for (size_t i = 0; i < s->n; ++i) {
        const int row = (int)((double)i * sy);
        for (size_t p = s->off[i]; p < s->off[i + 1]; ++p) {
            const int col = (int)((double)s->adj[p] * sx);
            unsigned *cell = &count[(size_t)row * (size_t)c.w + (size_t)col];
            ++*cell;
            top = *cell > top ? *cell : top;
        }
    }
    /* A matrix smaller than the canvas draws blocks, not single pixels. */
    const double scale = top > 0 ? 1.0 / log1p((double)top) : 0.0;
    const Rgb navy = { 20, 40, 130 };
    for (size_t i = 0; i < s->n && sx > 1.0; ++i) {            /* fill blocks for small graphs */
        for (size_t p = s->off[i]; p < s->off[i + 1]; ++p) {
            const int x0 = (int)((double)s->adj[p] * sx), y0 = (int)((double)i * sy);
            for (int yy = y0; yy < y0 + (int)sy && yy < c.h; ++yy) {
                for (int xx = x0; xx < x0 + (int)sx && xx < c.w; ++xx) {
                    blend(&c, xx, yy, navy, 1.0f);
                }
            }
        }
    }
    if (sx <= 1.0) {
        for (size_t k = 0; k < pixels; ++k) {
            if (count[k] > 0) {
                const float a = (float)pow(log1p((double)count[k]) * scale, 0.7);
                blend(&c, (int)(k % (size_t)c.w), (int)(k / (size_t)c.w), navy, 0.15f + 0.85f * a);
            }
        }
    }
    free(count);
    const GraphStatus st = canvas_png(f, &c);
    free(c.px);
    return st;
}

/* ------------------------------------------------------------------ */
/*  Public entry points                                                 */
/* ------------------------------------------------------------------ */

static VizFormat resolve_format(VizFormat requested, const char *path, size_t n, size_t drawn,
                                VizLayout layout)
{
    if (layout == VIZ_LAYOUT_MATRIX) {
        return VIZ_FORMAT_PNG;
    }
    if (requested != VIZ_FORMAT_AUTO) {
        return requested;
    }
    if (path != NULL) {
        const size_t len = strlen(path);
        if (len >= 4 && strcmp(path + len - 4, ".svg") == 0) { return VIZ_FORMAT_SVG; }
        if (len >= 4 && strcmp(path + len - 4, ".png") == 0) { return VIZ_FORMAT_PNG; }
    }
    return n <= AUTO_SVG_VERTICES && drawn <= AUTO_SVG_EDGES ? VIZ_FORMAT_SVG : VIZ_FORMAT_PNG;
}

static GraphStatus visualize(Graph *g, FILE *out, const char *path, const VizOptions *opt)
{
    if (g == NULL || out == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    VizOptions o = opt != NULL ? *opt : graph_viz_default_options();
    if (!options_valid(&o)) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    Snap s;
    GraphStatus st = snap_build(g, o.show_weights, &s);
    if (st != GRAPH_OK) {
        return st;
    }
    const size_t drawn = count_drawn(&s);
    const VizFormat fmt = resolve_format(o.format, path, s.n, drawn, o.layout);

    if (o.layout == VIZ_LAYOUT_MATRIX) {
        st = render_matrix_png(out, &s, &o);
        snap_free(&s);
        return st;
    }
    double *px  = malloc((s.n + 1) * sizeof(double));
    double *py  = malloc((s.n + 1) * sizeof(double));
    Rgb    *col = malloc((s.n + 1) * sizeof(Rgb));
    int    *slot = malloc((s.n + 1) * sizeof(int));
    int    *order = order_by_degree(&s);
    Legend  lg;
    memset(&lg, 0, sizeof lg);
    if (px == NULL || py == NULL || col == NULL || slot == NULL || order == NULL) {
        free(px); free(py); free(col); free(slot); free(order); snap_free(&s);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    st = layout_force(&s, o.iterations > 0 ? o.iterations : auto_iterations(s.n), o.seed, px, py);
    if (st == GRAPH_OK) {
        st = assign_colors(g, &s, o.color, col, slot, &lg);
    }
    if (st == GRAPH_OK) {
        st = fmt == VIZ_FORMAT_SVG ? render_svg(out, &s, px, py, col, slot, &lg, order, &o, drawn)
                                   : render_png(out, &s, px, py, col, slot, &lg, order, &o, drawn);
    }
    free(px); free(py); free(col); free(slot); free(order); snap_free(&s);
    return st;
}

GraphStatus graph_visualize_to(Graph *g, FILE *out, const VizOptions *opt)
{
    const GraphStatus st = visualize(g, out, NULL, opt);
    return st == GRAPH_OK ? GRAPH_OK : graph_set_error(st);
}

GraphStatus graph_visualize(Graph *g, const char *path, const VizOptions *opt)
{
    if (g == NULL || path == NULL) {
        return graph_set_error(GRAPH_ERR_INVALID_ARGUMENT);
    }
    FILE *f = fopen(path, "wb");
    if (f == NULL) {
        return graph_set_error(GRAPH_ERR_IO);
    }
    GraphStatus st = visualize(g, f, path, opt);
    if (fclose(f) != 0 && st == GRAPH_OK) {
        st = GRAPH_ERR_IO;
    }
    if (st != GRAPH_OK) {
        (void)remove(path);                 /* do not leave a half-written image behind */
        return graph_set_error(st);
    }
    return GRAPH_OK;
}
