/*
 * stub_backend.c — dense-matrix test backend.  See stub_backend.h.
 */

#include "stub_backend.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "list_backend.h"

typedef struct {
    size_t    n;          /* ID space size; the matrices are n x n */
    bool      directed;
    bool     *alive;
    bool     *has;        /* has[u * n + v] */
    double   *w;
    int      *free_stack; /* dead IDs; most recently removed is reused first */
    size_t    n_free;
    size_t    free_cap;
    size_t    edges;
    uint64_t  version;
} Stub;

static StubCostMode g_mode  = STUB_COST_MATRIX;
static double       g_scale = 1.0;
static bool         g_fail_create = false;

void stub_set_cost_mode(StubCostMode mode, double scale)
{
    g_mode  = mode;
    g_scale = scale;
}

void stub_set_fail_create(bool fail)
{
    g_fail_create = fail;
}

#define AT(s, u, v) ((size_t)(u) * (s)->n + (size_t)(v))

static bool valid(const Stub *s, int v)
{
    return s != NULL && v >= 0 && (size_t)v < s->n && s->alive[v];
}

/* ---- lifecycle ---- */

static void *s_create(size_t vertices, bool directed)
{
    if (g_fail_create) {
        return NULL;
    }
    Stub *s = calloc(1, sizeof(*s));
    if (s == NULL) {
        return NULL;
    }
    s->n        = vertices;
    s->directed = directed;
    s->alive    = calloc(vertices != 0 ? vertices : 1, sizeof(*s->alive));
    s->has      = calloc(vertices * vertices + 1, sizeof(*s->has));
    s->w        = calloc(vertices * vertices + 1, sizeof(*s->w));
    if (s->alive == NULL || s->has == NULL || s->w == NULL) {
        free(s->alive); free(s->has); free(s->w); free(s);
        return NULL;
    }
    for (size_t i = 0; i < vertices; ++i) {
        s->alive[i] = true;
    }
    return s;
}

static void s_destroy(void *impl)
{
    Stub *s = impl;
    if (s == NULL) {
        return;
    }
    free(s->alive); free(s->has); free(s->w); free(s->free_stack); free(s);
}

/* ---- topology ---- */

static size_t   s_vertex_count(const void *i) { const Stub *s = i; return s->n - s->n_free; }
static size_t   s_id_bound(const void *i)     { return ((const Stub *)i)->n; }
static size_t   s_edge_count(const void *i)   { return ((const Stub *)i)->edges; }
static bool     s_is_directed(const void *i)  { return ((const Stub *)i)->directed; }
static bool     s_is_alive(const void *i, int v) { return valid(i, v); }
static uint64_t s_version(const void *i)      { return ((const Stub *)i)->version; }

static size_t s_degree(const void *i, int v)
{
    const Stub *s = i;
    size_t d = 0;
    if (valid(s, v)) {
        for (size_t j = 0; j < s->n; ++j) {
            d += s->has[AT(s, v, j)];
        }
    }
    return d;
}

static size_t s_memory_bytes(const void *i)
{
    const Stub *s = i;
    return s->n * s->n * (sizeof(bool) + sizeof(double)) + s->n;
}

/* ---- mutation ---- */

static int s_add_vertex(void *impl)
{
    Stub *s = impl;
    if (s->n_free != 0) {
        const int id = s->free_stack[--s->n_free];
        s->alive[id] = true;
        ++s->version;
        return id;
    }
    const size_t m = s->n + 1;
    bool   *alive = calloc(m, sizeof(*alive));
    bool   *has   = calloc(m * m, sizeof(*has));
    double *w     = calloc(m * m, sizeof(*w));
    if (alive == NULL || has == NULL || w == NULL) {
        free(alive); free(has); free(w);
        return -1;
    }
    memcpy(alive, s->alive, s->n * sizeof(*alive));
    for (size_t u = 0; u < s->n; ++u) {
        memcpy(&has[u * m], &s->has[u * s->n], s->n * sizeof(*has));
        memcpy(&w[u * m],   &s->w[u * s->n],   s->n * sizeof(*w));
    }
    free(s->alive); free(s->has); free(s->w);
    s->alive = alive; s->has = has; s->w = w;
    s->n = m;
    s->alive[m - 1] = true;
    ++s->version;
    return (int)(m - 1);
}

static bool s_remove_vertex(void *impl, int v)
{
    Stub *s = impl;
    if (!valid(s, v)) {
        return false;
    }
    if (s->n_free == s->free_cap) {
        const size_t cap = s->free_cap != 0 ? s->free_cap * 2 : 8;
        int *stack = realloc(s->free_stack, cap * sizeof(*stack));
        if (stack == NULL) {
            return false;
        }
        s->free_stack = stack;
        s->free_cap   = cap;
    }
    for (size_t i = 0; i < s->n; ++i) {
        if (s->has[AT(s, v, i)]) {
            --s->edges;
            s->has[AT(s, v, i)] = false;
        }
        if (i != (size_t)v && s->has[AT(s, i, v)]) {
            if (s->directed) {
                --s->edges;   /* undirected mirrors were counted with the row */
            }
            s->has[AT(s, i, v)] = false;
        }
    }
    s->alive[v] = false;
    s->free_stack[s->n_free++] = v;
    ++s->version;
    return true;
}

static bool s_add_edge(void *impl, int u, int v, double weight)
{
    Stub *s = impl;
    if (!valid(s, u) || !valid(s, v) || !isfinite(weight) || s->has[AT(s, u, v)]) {
        return false;
    }
    s->has[AT(s, u, v)] = true;
    s->w[AT(s, u, v)]   = weight;
    if (!s->directed && u != v) {
        s->has[AT(s, v, u)] = true;
        s->w[AT(s, v, u)]   = weight;
    }
    ++s->edges;
    ++s->version;
    return true;
}

static bool s_remove_edge(void *impl, int u, int v)
{
    Stub *s = impl;
    if (!valid(s, u) || !valid(s, v) || !s->has[AT(s, u, v)]) {
        return false;
    }
    s->has[AT(s, u, v)] = false;
    if (!s->directed && u != v) {
        s->has[AT(s, v, u)] = false;
    }
    --s->edges;
    ++s->version;
    return true;
}

static bool s_has_edge(const void *impl, int u, int v, double *weight)
{
    const Stub *s = impl;
    if (!valid(s, u) || !valid(s, v) || !s->has[AT(s, u, v)]) {
        return false;
    }
    if (weight != NULL) {
        *weight = s->w[AT(s, u, v)];
    }
    return true;
}

/* The stub never compacts: dead IDs just wait on the free stack. */
static bool s_needs_compaction(const void *impl) { (void)impl; return false; }
static bool s_compact(void *impl, int *remap)    { (void)impl; (void)remap; return false; }

/* ---- iteration: a column cursor over the vertex's row ---- */

static void s_iter_begin(const void *impl, int v, GraphIter *it)
{
    const Stub *s = impl;
    it->impl   = impl;
    it->vertex = v;
    it->index  = 0;
    it->end    = valid(s, v) ? s->n : 0;
    it->a = it->b = NULL;
}

static bool s_iter_next(GraphIter *it, int *neighbor, double *weight)
{
    const Stub *s = it->impl;
    while (it->index < it->end) {
        const size_t j = it->index++;
        if (s->has[AT(s, it->vertex, j)]) {
            *neighbor = (int)j;
            *weight   = s->w[AT(s, it->vertex, j)];
            return true;
        }
    }
    return false;
}

/* ---- cost model ---- */

static double s_estimate_cost(const GraphShape *shape, const Workload *wl)
{
    if (g_mode == STUB_COST_SCALED) {
        return g_scale * list_backend()->estimate_cost(shape, wl);
    }
    const double n = (double)shape->id_bound;
    double cost = 0.0;
    cost += wl->n[OP_HAS_EDGE]      * 2.0;
    cost += wl->n[OP_ADD_EDGE]      * 3.0;
    cost += wl->n[OP_REMOVE_EDGE]   * 3.0;
    cost += wl->n[OP_ADD_VERTEX]    * (10.0 + 0.5 * n);
    cost += wl->n[OP_REMOVE_VERTEX] * (10.0 + 2.0 * n);
    cost += wl->n[OP_NEIGHBOR_SCAN] * (10.0 + 0.3 * n);
    cost += wl->n[OP_ALGO_RUN]      * (0.5 * n * n + 10.0 * n);
    cost += wl->n[OP_WEIGHTED_RUN]  * (0.5 * n * n + 10.0 * n);
    return cost;
}

static size_t s_estimate_memory(const GraphShape *shape)
{
    if (g_mode == STUB_COST_SCALED) {
        return list_backend()->estimate_memory(shape);
    }
    return shape->id_bound * shape->id_bound * (sizeof(bool) + sizeof(double));
}

static double s_estimate_build_cost(const GraphShape *shape)
{
    if (g_mode == STUB_COST_SCALED) {
        return list_backend()->estimate_build_cost(shape);
    }
    const double n = (double)shape->id_bound;
    return n * n + 5.0 * shape_entries(shape) + 10.0 * n;
}

/* ---- vtable ---- */

static const GraphBackend STUB_BACKEND = {
    .name                = "matrix-stub",
    .create              = s_create,
    .destroy             = s_destroy,
    .vertex_count        = s_vertex_count,
    .id_bound            = s_id_bound,
    .edge_count          = s_edge_count,
    .is_directed         = s_is_directed,
    .is_alive            = s_is_alive,
    .degree              = s_degree,
    .version             = s_version,
    .memory_bytes        = s_memory_bytes,
    .add_vertex          = s_add_vertex,
    .remove_vertex       = s_remove_vertex,
    .add_edge            = s_add_edge,
    .remove_edge         = s_remove_edge,
    .has_edge            = s_has_edge,
    .needs_compaction    = s_needs_compaction,
    .compact             = s_compact,
    .iter_begin          = s_iter_begin,
    .iter_next           = s_iter_next,
    .estimate_cost       = s_estimate_cost,
    .estimate_memory     = s_estimate_memory,
    .estimate_build_cost = s_estimate_build_cost,
};

const GraphBackend *stub_backend(void)
{
    return &STUB_BACKEND;
}
