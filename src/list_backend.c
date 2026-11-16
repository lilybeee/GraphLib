/*
 * list_backend.c — adapter that fills a GraphBackend vtable with the
 * adjacency-list implementation.
 *
 * The wrappers only cast the opaque handle; all the work happens in
 * adjacency_list.c.  The one non-trivial piece is the cost model at the
 * bottom, which is the list's half of the observer's decision.
 */

#include <math.h>
#include "list_backend.h"

#include "adjacency_list.h"

/* ------------------------------------------------------------------ */
/*  Lifecycle and topology                                              */
/* ------------------------------------------------------------------ */

static void *b_create_unweighted(size_t vertices, bool directed)
{
    return list_graph_create_unweighted(vertices, directed);
}

static bool b_is_unweighted(const void *impl)
{
    return list_graph_is_unweighted(impl);
}

static void *b_create(size_t vertices, bool directed)
{
    return list_graph_create(vertices, directed);
}

static void b_destroy(void *impl)                 { list_graph_destroy(impl); }
static size_t b_vertex_count(const void *impl)    { return list_graph_vertex_count(impl); }
static size_t b_id_bound(const void *impl)        { return list_graph_id_bound(impl); }
static size_t b_edge_count(const void *impl)      { return list_graph_edge_count(impl); }
static bool b_is_directed(const void *impl)       { return list_graph_is_directed(impl); }
static bool b_is_alive(const void *impl, int v)   { return list_graph_is_alive(impl, v); }
static size_t b_degree(const void *impl, int v)   { return list_graph_degree(impl, v); }
static uint64_t b_version(const void *impl)       { return list_graph_version(impl); }
static size_t b_memory_bytes(const void *impl)    { return list_graph_memory_bytes(impl); }

/* ------------------------------------------------------------------ */
/*  Mutation                                                            */
/* ------------------------------------------------------------------ */

static int b_add_vertex(void *impl)               { return list_graph_add_vertex(impl); }
static bool b_remove_vertex(void *impl, int v)    { return list_graph_remove_vertex(impl, v); }

static bool b_add_edge(void *impl, int u, int v, double w)
{
    return list_graph_add_edge(impl, u, v, w);
}

static bool b_add_new_edge(void *impl, int u, int v, double w)
{
    return list_graph_add_new_edge(impl, u, v, w);
}

static bool b_remove_edge(void *impl, int u, int v)
{
    return list_graph_remove_edge(impl, u, v);
}

static bool b_has_edge(const void *impl, int u, int v, double *w)
{
    return list_graph_has_edge(impl, u, v, w);
}

static bool b_needs_compaction(const void *impl)  { return list_graph_needs_compaction(impl); }
static bool b_compact(void *impl, int *remap)     { return list_graph_compact(impl, remap); }

/* ------------------------------------------------------------------ */
/*  Neighbour iteration                                                 */
/* ------------------------------------------------------------------ */

/* The iterator just walks the vertex's two arrays, so advancing is a pair of
 * loads and an increment. */
static void b_iter_begin(const void *impl, int v, GraphIter *it)
{
    const int    *to;
    const double *w;
    const size_t  deg = list_graph_neighbors(impl, v, &to, &w);

    it->impl   = impl;
    it->vertex = v;
    it->index  = 0;
    it->end    = deg;
    it->a      = to;
    it->b      = w;
}

static bool b_iter_next(GraphIter *it, int *neighbor, double *weight)
{
    if (it->index >= it->end) {
        return false;
    }
    *neighbor = ((const int *)it->a)[it->index];
    if (weight != NULL) {
        *weight = it->b != NULL ? ((const double *)it->b)[it->index] : 1.0;
    }
    ++it->index;
    return true;
}

/* ------------------------------------------------------------------ */
/*  Algorithm kernels (all specialised: they walk the arrays directly)  */
/* ------------------------------------------------------------------ */

static size_t b_bfs(const void *impl, int s, int *order, int *parent)
{
    return list_bfs(impl, s, order, parent);
}

static size_t b_dfs(const void *impl, int s, int *order, int *parent)
{
    return list_dfs(impl, s, order, parent);
}

static size_t b_components(const void *impl, int *component)
{
    return list_connected_components(impl, component);
}

static bool b_has_cycle(const void *impl)         { return list_has_cycle(impl); }

static bool b_topological_sort(const void *impl, int *order)
{
    return list_topological_sort(impl, order);
}

static bool b_dijkstra(const void *impl, int s, double *distance, int *parent)
{
    return list_dijkstra(impl, s, distance, parent);
}

static double b_mst(const void *impl, GraphEdge *edges, size_t *count)
{
    return list_mst_kruskal(impl, edges, count);
}

/* ------------------------------------------------------------------ */
/*  Cost model                                                          */
/* ------------------------------------------------------------------ */

/*
 * Per-operation costs in nanoseconds, fitted to timings of the real list and
 * matrix on random graphs of 300, 1000 and 3000 vertices at average degrees
 * 2 to 250 (single machine, about +-20% run to run, so these are shapes, not
 * promises):
 *
 *   point query   28 ns + ~1 ns per neighbour scanned (one cache miss to reach
 *                 the vertex's block, then a linear pass over ints);
 *   neighbour walk  35 ns + ~4 ns per neighbour (it also reads the weight);
 *   traversal     ~3 ns per stored entry + ~5 ns per vertex (BFS, DFS,
 *                 components; the arrays are walked in order);
 *   weighted run  ~6.5 ns per entry + ~120 ns per vertex (the heap).
 *
 * What matters for the decision is the shape of each formula: a list pays for
 * the degree on every point query and for the whole edge set on every
 * whole-graph algorithm, but it is cheap to iterate, and its weights sit next
 * to the neighbour IDs so weighted algorithms stay fast.
 */
enum {
    LIST_SLOT_BYTES  = 32,  /* sizeof(Slot) in adjacency_list.c            */
    LIST_ENTRY_BYTES = 12   /* one int + one double per stored entry       */
};

static double list_estimate_cost(const GraphShape *s, const Workload *w)
{
    const double entries = shape_entries(s);
    const double deg     = shape_avg_degree(s);
    const double n       = (double)s->live;

    double cost = 0.0;
    cost += w->n[OP_HAS_EDGE]      * (28.0 + 1.0 * deg);
    cost += w->n[OP_ADD_EDGE]      * (48.0 + 1.0 * deg);          /* dup scan + push */
    cost += w->n[OP_REMOVE_EDGE]   * (38.0 + 1.0 * deg);
    cost += w->n[OP_ADD_VERTEX]    * 20.0;
    /* Undirected: touch each neighbour's list.  Directed: scan every list. */
    cost += w->n[OP_REMOVE_VERTEX] *
            (s->directed ? 50.0 + 1.0 * entries + 28.0 * n
                         : 50.0 + deg * (40.0 + deg));
    cost += w->n[OP_NEIGHBOR_SCAN] * (35.0 + 4.2 * deg);
    cost += w->n[OP_ALGO_RUN]      * (3.0 * entries + 5.0 * n);
    cost += w->n[OP_WEIGHTED_RUN]  * (6.5 * entries + 120.0 * n);
    /* The three kinds below were calibrated against bench_list / bench_matrix at
     * -O2; only the list-to-matrix ratio matters.  ITERATIVE is priced as a
     * nominal 25 sweeps. */
    cost += w->n[OP_ITERATIVE_RUN] * 25.0 * (12.0 * entries + 60.0 * n);
    cost += w->n[OP_ALL_PAIRS_RUN] * n * (3.0 * entries + 5.0 * n);
    /* Probing every pair of one vertex's neighbours costs a binary search each. */
    cost += w->n[OP_ADJACENCY_RUN] *
            (n * deg * deg * 0.5 * (4.0 + 1.1 * log2(deg + 2.0)) + 35.0 * n);
    return cost;
}

static size_t list_estimate_memory(const GraphShape *s)
{
    /* Slots for the whole ID space, plus entries with ~50% growth slack. */
    return s->id_bound * LIST_SLOT_BYTES +
           (size_t)(shape_entries(s) * (s->unweighted ? LIST_ENTRY_BYTES - 8 : LIST_ENTRY_BYTES) * 1.5);
}

static double list_estimate_build_cost(const GraphShape *s)
{
    /* Measured by rebuilding a matrix as a list (the only migration that ends
     * here): ~70 ns per entry and ~300 ns per vertex.  Most of that is reading
     * the source matrix (a cache miss per weight); bulk loading itself skips
     * the duplicate scan and costs a few nanoseconds per append. */
    return 70.0 * shape_entries(s) + 300.0 * (double)s->id_bound;
}

/* ------------------------------------------------------------------ */
/*  The vtable                                                          */
/* ------------------------------------------------------------------ */

static const GraphBackend LIST_BACKEND = {
    .name                 = "list",
    .create               = b_create,
    .create_unweighted    = b_create_unweighted,
    .is_unweighted        = b_is_unweighted,
    .destroy              = b_destroy,
    .vertex_count         = b_vertex_count,
    .id_bound             = b_id_bound,
    .edge_count           = b_edge_count,
    .is_directed          = b_is_directed,
    .is_alive             = b_is_alive,
    .degree               = b_degree,
    .version              = b_version,
    .memory_bytes         = b_memory_bytes,
    .add_vertex           = b_add_vertex,
    .remove_vertex        = b_remove_vertex,
    .add_edge             = b_add_edge,
    .remove_edge          = b_remove_edge,
    .has_edge             = b_has_edge,
    .needs_compaction     = b_needs_compaction,
    .compact              = b_compact,
    .iter_begin           = b_iter_begin,
    .iter_next            = b_iter_next,
    .add_new_edge         = b_add_new_edge,
    .bfs                  = b_bfs,
    .dfs                  = b_dfs,
    .connected_components = b_components,
    .has_cycle            = b_has_cycle,
    .topological_sort     = b_topological_sort,
    .dijkstra             = b_dijkstra,
    .mst_kruskal          = b_mst,
    .estimate_cost        = list_estimate_cost,
    .estimate_memory      = list_estimate_memory,
    .estimate_build_cost  = list_estimate_build_cost,
};

const GraphBackend *list_backend(void)
{
    return &LIST_BACKEND;
}
