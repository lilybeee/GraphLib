/*
 * matrix_backend.c — adapter that fills a GraphBackend vtable with the
 * adjacency-matrix implementation.
 *
 * The wrappers only cast the opaque handle; all the work happens in
 * adjacency_matrix.c.  The non-trivial piece is the cost model at the bottom,
 * which is the matrix's half of the observer's decision.
 */

#include "matrix_backend.h"

#include <math.h>

#include "adjacency_matrix.h"

/* ------------------------------------------------------------------ */
/*  Lifecycle and topology                                              */
/* ------------------------------------------------------------------ */

static void *b_create_unweighted(size_t vertices, bool directed) {
    return matrix_graph_create_unweighted(vertices, directed);
}

static bool b_is_unweighted(const void *impl) { return matrix_graph_is_unweighted(impl); }

static void *b_create(size_t vertices, bool directed) {
    return matrix_graph_create(vertices, directed);
}

static void b_destroy(void *impl)               { matrix_graph_destroy(impl); }
static size_t b_vertex_count(const void *impl)  { return matrix_graph_vertex_count(impl); }
static size_t b_id_bound(const void *impl)      { return matrix_graph_id_bound(impl); }
static size_t b_edge_count(const void *impl)    { return matrix_graph_edge_count(impl); }
static bool b_is_directed(const void *impl)     { return matrix_graph_is_directed(impl); }
static bool b_is_alive(const void *impl, int v) { return matrix_graph_is_alive(impl, v); }
static size_t b_degree(const void *impl, int v) { return matrix_graph_degree(impl, v); }
static uint64_t b_version(const void *impl)     { return matrix_graph_version(impl); }
static size_t b_memory_bytes(const void *impl)  { return matrix_graph_memory_bytes(impl); }

/* ------------------------------------------------------------------ */
/*  Mutation                                                            */
/* ------------------------------------------------------------------ */

static int b_add_vertex(void *impl)             { return matrix_graph_add_vertex(impl); }
static bool b_remove_vertex(void *impl, int v)  { return matrix_graph_remove_vertex(impl, v); }

static bool b_add_edge(void *impl, int u, int v, double w) {
    return matrix_graph_add_edge(impl, u, v, w);
}

static bool b_remove_edge(void *impl, int u, int v) {
    return matrix_graph_remove_edge(impl, u, v);
}

static bool b_has_edge(const void *impl, int u, int v, double *w) {
    return matrix_graph_has_edge(impl, u, v, w);
}

static bool b_needs_compaction(const void *impl) { return matrix_graph_needs_compaction(impl); }
static bool b_compact(void *impl, int *remap)    { return matrix_graph_compact(impl, remap); }

/* ------------------------------------------------------------------ */
/*  Neighbour iteration                                                 */
/* ------------------------------------------------------------------ */

/* The iterator keeps a column cursor; each step finds the next set bit. */
static void b_iter_begin(const void *impl, int v, GraphIter *it) {
    it->impl = impl;
    it->vertex = v;
    it->index = 0;
    it->end = 0;
    it->a = NULL;
    it->b = NULL;
}

static bool b_iter_next(GraphIter *it, int *neighbor, double *weight) {
    size_t cursor = it->index;
    const bool found = matrix_graph_next_neighbor(it->impl, it->vertex, &cursor, neighbor, weight);
    it->index = cursor;
    return found;
}

/* ------------------------------------------------------------------ */
/*  Algorithm kernels (all specialised: they scan the bitset rows)      */
/* ------------------------------------------------------------------ */

static size_t b_bfs(const void *impl, int s, int *order, int *parent) {
    return matrix_bfs(impl, s, order, parent);
}

static size_t b_dfs(const void *impl, int s, int *order, int *parent) {
    return matrix_dfs(impl, s, order, parent);
}

static size_t b_components(const void *impl, int *component) {
    return matrix_connected_components(impl, component);
}

static bool b_has_cycle(const void *impl) { return matrix_has_cycle(impl); }

static bool b_topological_sort(const void *impl, int *order) {
    return matrix_topological_sort(impl, order);
}

static bool b_dijkstra(const void *impl, int s, double *distance, int *parent) {
    return matrix_dijkstra(impl, s, distance, parent);
}

static double b_mst(const void *impl, GraphEdge *edges, size_t *count) {
    return matrix_mst_kruskal(impl, edges, count);
}

/* ------------------------------------------------------------------ */
/*  Cost model                                                          */
/* ------------------------------------------------------------------ */

/*
 * Per-operation costs in nanoseconds, fitted to timings of the real list and
 * matrix on random graphs of 300, 1000 and 3000 vertices at average degrees
 * 2 to 250 (single machine, about +-20% run to run, so these are shapes, not
 * promises).  What the measurements say:
 *
 *   point query   ~8-15 ns whatever the density (two loads and a bit test);
 *   traversal     ~(6 + 1.9 * words) ns per vertex, where words = V / 64 is
 *                 the row width, and almost nothing per edge: a matrix wins
 *                 once rows are dense enough that the list's per-edge cost
 *                 outweighs the matrix's per-row cost (degree >~ 10-15% of V);
 *   weighted run  each edge's weight is its own cache miss unless the row is
 *                 dense, so the per-edge cost grows with the row width and
 *                 shrinks with the degree: a matrix loses Dijkstra to a list
 *                 unless the graph is very dense;
 *   first use of  a 32 KB tile allocation (~8 us) is charged to add_edge in
 *   a block       proportion to the chance the edge lands in an empty block.
 */
enum {
    MATRIX_TILE_BYTES = 64 * 64 * 8   /* sizeof(Tile) in adjacency_matrix.c */
};

/** Number of 64-column words per row for an ID space of this size. */
static double row_words(const GraphShape *s) {
    return ceil((double)s->id_bound / 64.0);
}

/**
 * Expected number of allocated tiles.  If `entries` cells land uniformly at
 * random in `total` blocks, the expected number of non-empty blocks is
 * total * (1 - e^(-entries / total)).  Real graphs cluster, so this is a
 * pessimistic (upper) estimate.
 */
static double expected_tiles(const GraphShape *s) {
    const double words = row_words(s);
    const double total = words * words;
    if (total == 0.0) {
        return 0.0;
    }
    return total * (1.0 - exp(-shape_entries(s) / total));
}

static double matrix_estimate_cost(const GraphShape *s, const Workload *w) {
    const double entries = shape_entries(s);
    const double deg = shape_avg_degree(s);
    const double n = (double)s->live;
    const double words = row_words(s);

    const double tiles = words * words;
    const double new_tile = tiles == 0.0 ? 0.0 : exp(-entries / tiles);

    double cost = 0.0;
    cost += w->n[OP_HAS_EDGE]      * (8.0 + 0.15 * words);
    cost += w->n[OP_ADD_EDGE]      * (12.0 + 8000.0 * new_tile);
    cost += w->n[OP_REMOVE_EDGE]   * (8.0 + 0.15 * words);
    cost += w->n[OP_ADD_VERTEX]    * 20.0;
    /* Row scan (V/64 words) plus, for a directed graph, the whole column. */
    cost += w->n[OP_REMOVE_VERTEX] *
            (s->directed ? 60.0 + 1.5 * (double)s->id_bound + 2.0 * words + 20.0 * deg
                         : 60.0 + 2.0 * words + 20.0 * deg);
    /* Walking a row also reads each neighbour's weight (a cache miss each). */
    cost += w->n[OP_NEIGHBOR_SCAN] * (60.0 + 4.0 * words + deg * (10.0 + 1.6 * words));
    cost += w->n[OP_ALGO_RUN]      * (n * (6.0 + 1.9 * words) + 0.3 * entries);
    cost += w->n[OP_WEIGHTED_RUN]  *
            (entries * (6.0 + words * (0.2 + 25.0 / (deg + 5.0))) +
             n * (120.0 + 2.0 * words));
    /* Calibrated against bench_list / bench_matrix at -O2 (see list_estimate_cost). */
    cost += w->n[OP_ITERATIVE_RUN] * 25.0 * (27.0 * entries + n * (100.0 + 12.0 * words));
    cost += w->n[OP_ALL_PAIRS_RUN] * n * (n * (6.0 + 1.9 * words) + 0.3 * entries);
    /* has_edge is a single bit test, so pair probing is cheap on a matrix. */
    cost += w->n[OP_ADJACENCY_RUN] *
            (n * deg * deg * 0.5 * (2.5 + 0.05 * words) +
             n * (60.0 + 4.0 * words + deg * (10.0 + 1.6 * words)));
    return cost;
}

/**
 * Expected number of allocated bit strips.  A strip exists once any of its 64
 * rows has an edge; with `entries` cells spread over `words` strips that is
 * words * (1 - e^(-entries / words)) of them.
 */
static double expected_strips(const GraphShape *s) {
    const double words = row_words(s);
    return words == 0.0 ? 0.0 : words * (1.0 - exp(-shape_entries(s) / words));
}

static size_t matrix_estimate_memory(const GraphShape *s) {
    const double words = row_words(s);
    /* Weight tiles (and the sparse index that finds them) exist only for a weighted graph. */
    const double tiles = s->unweighted ? 0.0 : expected_tiles(s);
    return (size_t)(tiles * (MATRIX_TILE_BYTES + 4 * sizeof(uint64_t)) +
                    expected_strips(s) * 64.0 * words * sizeof(uint64_t) +  /* bits */
                    words * sizeof(void *) +              /* strip pointers   */
                    words * sizeof(uint64_t) +            /* alive bits       */
                    words * 64.0 * sizeof(uint32_t));     /* degrees          */
}

static double matrix_estimate_build_cost(const GraphShape *s) {
    const double words = row_words(s);
    /* Measured by rebuilding a list as a matrix: ~14 us per tile (an
     * allocation plus the first touch of about eight pages) and ~16 ns per
     * entry for the bit set and weight store. */
    return 14000.0 * expected_tiles(s) + 16.0 * shape_entries(s) +
           10.0 * (double)s->id_bound + 2.0 * words * words;
}

/* ------------------------------------------------------------------ */
/*  The vtable                                                          */
/* ------------------------------------------------------------------ */

/* add_new_edge is left NULL: the duplicate check is one bit test, so there is
 * nothing for a bulk loader to skip. */
static const GraphBackend MATRIX_BACKEND = {
    .name                 = "matrix",
    .create               = b_create,
    .create_unweighted    = b_create_unweighted,
    .is_unweighted        = b_is_unweighted,
    .max_vertices         = matrix_max_vertices,
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
    .bfs                  = b_bfs,
    .dfs                  = b_dfs,
    .connected_components = b_components,
    .has_cycle            = b_has_cycle,
    .topological_sort     = b_topological_sort,
    .dijkstra             = b_dijkstra,
    .mst_kruskal          = b_mst,
    .estimate_cost        = matrix_estimate_cost,
    .estimate_memory      = matrix_estimate_memory,
    .estimate_build_cost  = matrix_estimate_build_cost,
};

const GraphBackend *matrix_backend(void) {
    return &MATRIX_BACKEND;
}
