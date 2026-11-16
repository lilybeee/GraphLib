/*
 * graph_compare.c — benchmark one algorithm on both representations.
 * See graph_compare in graph.h.
 *
 * The graph is cloned into a list and a matrix, so the caller's graph, its
 * statistics and its cache are untouched and nothing here can trigger a
 * migration.  The timed code calls the backend kernels directly, which is
 * what the facade would call, minus the cache and observer bookkeeping.
 */

#define _POSIX_C_SOURCE 200809L

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "adjacency_list.h"
#include "adjacency_matrix.h"
#include "graph_internal.h"
#include "graph_status.h"
#include "list_algorithms.h"
#include "list_backend.h"
#include "matrix_algorithms.h"
#include "matrix_backend.h"
#include "migrate.h"

/* ------------------------------------------------------------------ */
/*  One representation's working set                                    */
/* ------------------------------------------------------------------ */

typedef struct {
    const GraphBackend *be;
    void               *impl;
    size_t              id_bound;
    size_t              vertices;
    int                *order;       /* vertices          */
    int                *parent;      /* id_bound          */
    int                *component;   /* id_bound          */
    double             *distance;    /* id_bound          */
    GraphEdge          *edges;       /* vertices - 1      */
    /* results of the last run */
    bool                ok;
    bool                flag;        /* has_cycle                       */
    size_t              count;       /* visited / components / mst edges */
    double              total;       /* MST weight                       */
    long long           triangles;   /* triangle count                   */
} Side;

static void side_free(Side *s)
{
    if (s->impl != NULL) {
        s->be->destroy(s->impl);
    }
    free(s->order);
    free(s->parent);
    free(s->component);
    free(s->distance);
    free(s->edges);
    memset(s, 0, sizeof(*s));
}

static bool side_init(Side *s, const GraphBackend *be, const GraphBackend *from, const void *src)
{
    memset(s, 0, sizeof(*s));
    s->be   = be;
    s->impl = graph_clone(from, src, be);
    if (s->impl == NULL) {
        return false;
    }
    s->id_bound  = be->id_bound(s->impl);
    s->vertices  = be->vertex_count(s->impl);
    const size_t ib = s->id_bound > 0 ? s->id_bound : 1;
    const size_t vc = s->vertices > 0 ? s->vertices : 1;
    s->order     = calloc(vc, sizeof(int));
    s->parent    = calloc(ib, sizeof(int));
    s->component = calloc(ib, sizeof(int));
    s->distance  = calloc(ib, sizeof(double));
    s->edges     = calloc(vc, sizeof(GraphEdge));
    if (s->order == NULL || s->parent == NULL || s->component == NULL ||
        s->distance == NULL || s->edges == NULL) {
        side_free(s);
        return false;
    }
    return true;
}

static void run_once(Side *s, GraphAlgo algo, int source)
{
    const GraphBackend *be = s->be;
    s->ok = true;
    switch (algo) {
    case GRAPH_ALGO_BFS:
        if (be->bfs == NULL) { s->ok = false; break; }
        s->count = be->bfs(s->impl, source, s->order, s->parent);
        break;
    case GRAPH_ALGO_DFS:
        if (be->dfs == NULL) { s->ok = false; break; }
        s->count = be->dfs(s->impl, source, s->order, s->parent);
        break;
    case GRAPH_ALGO_COMPONENTS:
        if (be->connected_components == NULL) { s->ok = false; break; }
        s->count = be->connected_components(s->impl, s->component);
        break;
    case GRAPH_ALGO_CYCLE:
        if (be->has_cycle == NULL) { s->ok = false; break; }
        s->flag = be->has_cycle(s->impl);
        break;
    case GRAPH_ALGO_TOPO:
        if (be->topological_sort == NULL) { s->ok = false; break; }
        s->ok = be->topological_sort(s->impl, s->order);
        break;
    case GRAPH_ALGO_DIJKSTRA:
        if (be->dijkstra == NULL) { s->ok = false; break; }
        s->ok = be->dijkstra(s->impl, source, s->distance, s->parent);
        break;
    case GRAPH_ALGO_MST:
        if (be->mst_kruskal == NULL) { s->ok = false; break; }
        s->count = 0;
        s->total = be->mst_kruskal(s->impl, s->edges, &s->count);
        s->ok    = isfinite(s->total);
        break;
    case GRAPH_ALGO_PAGERANK: {
        const bool m = be == matrix_backend();
        s->ok = (m ? matrix_pagerank(s->impl, 0.85, 1e-10, 200, s->distance, NULL)
                   : list_pagerank(s->impl, 0.85, 1e-10, 200, s->distance, NULL)) == GRAPH_OK;
        break;
    }
    case GRAPH_ALGO_TRIANGLES: {
        const bool m = be == matrix_backend();
        s->ok = (m ? matrix_triangle_count(s->impl, &s->triangles)
                   : list_triangle_count(s->impl, &s->triangles)) == GRAPH_OK;
        break;
    }
    case GRAPH_ALGO_CORE_NUMBERS: {
        const bool m = be == matrix_backend();
        s->ok = (m ? matrix_core_numbers(s->impl, s->component)
                   : list_core_numbers(s->impl, s->component)) == GRAPH_OK;
        break;
    }
    case GRAPH_ALGO_BELLMAN_FORD: {
        const bool m = be == matrix_backend();
        s->ok = (m ? matrix_bellman_ford(s->impl, source, s->distance, s->parent)
                   : list_bellman_ford(s->impl, source, s->distance, s->parent)) == GRAPH_OK;
        break;
    }
    default:
        s->ok = false;
        break;
    }
}

static double now_ns(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec * 1e9 + (double)ts.tv_nsec;
}

/** Best-of-N wall time of one run, in nanoseconds. */
static double time_side(Side *s, GraphAlgo algo, int source, int reps)
{
    double best = INFINITY;
    for (int i = 0; i < reps; ++i) {
        const double t0 = now_ns();
        run_once(s, algo, source);
        const double dt = now_ns() - t0;
        if (dt < best) {
            best = dt;
        }
    }
    return best;
}

/* ------------------------------------------------------------------ */
/*  Agreement                                                           */
/* ------------------------------------------------------------------ */

static bool close_enough(double a, double b)
{
    if (isinf(a) || isinf(b)) {
        return isinf(a) && isinf(b);
    }
    return fabs(a - b) <= 1e-9 * fmax(1.0, fmax(fabs(a), fabs(b)));
}

/** The vertices listed in a->order and b->order are the same set. */
static bool same_visit_set(const Side *a, const Side *b)
{
    if (a->count != b->count) {
        return false;
    }
    bool *mark = calloc(a->id_bound > 0 ? a->id_bound : 1, sizeof(bool));
    if (mark == NULL) {
        return false;
    }
    for (size_t i = 0; i < a->count; ++i) {
        mark[a->order[i]] = true;
    }
    bool same = true;
    for (size_t i = 0; i < b->count && same; ++i) {
        same = mark[b->order[i]];
    }
    free(mark);
    return same;
}

/** Component labels describe the same partition (labels may differ). */
static bool same_partition(const Side *a, const Side *b)
{
    if (a->count != b->count) {
        return false;
    }
    int *map = malloc((a->count > 0 ? a->count : 1) * sizeof(int));
    if (map == NULL) {
        return false;
    }
    for (size_t i = 0; i < a->count; ++i) {
        map[i] = -1;
    }
    bool same = true;
    for (size_t v = 0; v < a->id_bound && same; ++v) {
        const int ca = a->component[v];
        const int cb = b->component[v];
        if ((ca < 0) != (cb < 0)) {
            same = false;
        } else if (ca >= 0) {
            if (map[ca] < 0) {
                map[ca] = cb;
            } else if (map[ca] != cb) {
                same = false;
            }
        }
    }
    free(map);
    return same;
}

/** s->order is a valid topological order of s's graph. */
static bool valid_topo(const Side *s)
{
    int *pos = malloc((s->id_bound > 0 ? s->id_bound : 1) * sizeof(int));
    if (pos == NULL) {
        return false;
    }
    for (size_t i = 0; i < s->id_bound; ++i) {
        pos[i] = -1;
    }
    for (size_t i = 0; i < s->vertices; ++i) {
        pos[s->order[i]] = (int)i;
    }
    bool valid = true;
    for (size_t u = 0; u < s->id_bound && valid; ++u) {
        if (!s->be->is_alive(s->impl, (int)u)) {
            continue;
        }
        GraphIter it;
        int       v;
        s->be->iter_begin(s->impl, (int)u, &it);
        while (valid && s->be->iter_next(&it, &v, NULL)) {
            valid = pos[u] >= 0 && pos[v] >= 0 && pos[u] < pos[v];
        }
    }
    free(pos);
    return valid;
}

static bool sides_agree(GraphAlgo algo, const Side *a, const Side *b)
{
    if (a->ok != b->ok) {
        return false;
    }
    switch (algo) {
    case GRAPH_ALGO_BFS:
    case GRAPH_ALGO_DFS:
        return same_visit_set(a, b);
    case GRAPH_ALGO_COMPONENTS:
        return same_partition(a, b);
    case GRAPH_ALGO_CYCLE:
        return a->flag == b->flag;
    case GRAPH_ALGO_TOPO:
        return !a->ok || (valid_topo(a) && valid_topo(b));
    case GRAPH_ALGO_DIJKSTRA:
        if (!a->ok) {
            return true;
        }
        for (size_t v = 0; v < a->id_bound; ++v) {
            if (!close_enough(a->distance[v], b->distance[v])) {
                return false;
            }
        }
        return true;
    case GRAPH_ALGO_MST:
        return !a->ok || (a->count == b->count && close_enough(a->total, b->total));
    case GRAPH_ALGO_PAGERANK:
    case GRAPH_ALGO_BELLMAN_FORD:
        if (!a->ok) {
            return true;
        }
        for (size_t v = 0; v < a->id_bound; ++v) {
            if (a->be->is_alive(a->impl, (int)v) && !close_enough(a->distance[v], b->distance[v])) {
                return false;
            }
        }
        return true;
    case GRAPH_ALGO_TRIANGLES:
        return !a->ok || a->triangles == b->triangles;
    case GRAPH_ALGO_CORE_NUMBERS:
        if (!a->ok) {
            return true;
        }
        for (size_t v = 0; v < a->id_bound; ++v) {
            if (a->be->is_alive(a->impl, (int)v) && a->component[v] != b->component[v]) {
                return false;
            }
        }
        return true;
    default:
        return false;
    }
}

/* ------------------------------------------------------------------ */
/*  Public                                                              */
/* ------------------------------------------------------------------ */

const char *graph_algo_name(GraphAlgo algo)
{
    static const char *const names[GRAPH_ALGO_COUNT] = {
        "BFS", "DFS", "Connected components", "Cycle detection",
        "Topological sort", "Dijkstra", "Kruskal MST",
        "PageRank", "Triangle count", "k-core numbers", "Bellman-Ford",
    };
    return (unsigned)algo < GRAPH_ALGO_COUNT ? names[algo] : "?";
}

bool graph_compare(const Graph *g, GraphAlgo algo, int source, int repetitions,
                   GraphCompareResult *out)
{
    if (g == NULL || out == NULL || (unsigned)algo >= GRAPH_ALGO_COUNT ||
        g->be->vertex_count(g->impl) == 0) {
        return false;
    }
    if (repetitions <= 0) {
        repetitions = 5;
    }
    if (!g->be->is_alive(g->impl, source)) {
        source = -1;
        for (size_t v = 0; v < g->be->id_bound(g->impl); ++v) {
            if (g->be->is_alive(g->impl, (int)v)) {
                source = (int)v;
                break;
            }
        }
    }

    Side list, matrix;
    if (!side_init(&list, list_backend(), g->be, g->impl)) {
        return false;
    }
    if (!side_init(&matrix, matrix_backend(), g->be, g->impl)) {
        side_free(&list);
        return false;
    }

    memset(out, 0, sizeof(*out));
    out->algo        = algo;
    out->source      = source;
    out->repetitions = repetitions;
    out->list_bytes   = list.be->memory_bytes(list.impl);
    out->matrix_bytes = matrix.be->memory_bytes(matrix.impl);

    out->list_ns   = time_side(&list, algo, source, repetitions);
    out->matrix_ns = time_side(&matrix, algo, source, repetitions);
    out->list_ok   = list.ok;
    out->matrix_ok = matrix.ok;
    out->agreed    = sides_agree(algo, &list, &matrix);
    out->faster    = out->matrix_ns < out->list_ns ? GRAPH_MATRIX : GRAPH_LIST;

    side_free(&list);
    side_free(&matrix);
    return true;
}

bool graph_compare_all(const Graph *g, int source, int repetitions,
                       GraphCompareResult out[GRAPH_ALGO_COUNT])
{
    for (int a = 0; a < GRAPH_ALGO_COUNT; ++a) {
        if (!graph_compare(g, (GraphAlgo)a, source, repetitions, &out[a])) {
            return false;
        }
    }
    return true;
}

void graph_compare_print(const GraphCompareResult *r, FILE *out)
{
    if (r == NULL || out == NULL) {
        return;
    }
    fprintf(out, "Algorithm : %s", graph_algo_name(r->algo));
    if (r->algo == GRAPH_ALGO_BFS || r->algo == GRAPH_ALGO_DFS ||
        r->algo == GRAPH_ALGO_DIJKSTRA || r->algo == GRAPH_ALGO_BELLMAN_FORD) {
        fprintf(out, " (source %d)", r->source);
    }
    fprintf(out, "\n");
    fprintf(out, "  list  : %12.0f ns   %10zu bytes%s\n", r->list_ns, r->list_bytes,
            r->list_ok ? "" : "   (no result)");
    fprintf(out, "  matrix: %12.0f ns   %10zu bytes%s\n", r->matrix_ns, r->matrix_bytes,
            r->matrix_ok ? "" : "   (no result)");
    fprintf(out, "  faster: %s", graph_rep_name(r->faster));
    if (r->list_ns > 0.0 && r->matrix_ns > 0.0) {
        const double ratio = r->list_ns > r->matrix_ns ? r->list_ns / r->matrix_ns
                                                       : r->matrix_ns / r->list_ns;
        fprintf(out, " (%.2fx)", ratio);
    }
    fprintf(out, "   results: %s\n", r->agreed ? "equivalent" : "DIFFER");
}
