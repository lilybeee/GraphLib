/*
 * graph_build.c — the unified-API functions that create or combine graphs.
 * (The rest of the algorithm wrappers are generated; see graph_algorithms.h.)
 *
 * The library functions return a raw ListGraph / MatrixGraph.  Each wrapper here
 * runs the one matching the graph's current representation, then wraps the raw
 * result in a Graph with graph_adopt(), which takes ownership of it.
 *
 * C note: ListGraph * and MatrixGraph * are different types, so they cannot
 * share one out-parameter.  Each call uses a typed local, and the pointer is
 * stored as void * afterwards (any object pointer converts to void * and back).
 */

#include <stdlib.h>

#include "graph_advanced.h"
#include "graph_algorithms.h"

#include "adjacency_list.h"
#include "adjacency_matrix.h"
#include "graph_internal.h"
#include "list_advanced.h"
#include "list_algorithms.h"
#include "list_backend.h"
#include "matrix_advanced.h"
#include "matrix_algorithms.h"
#include "matrix_backend.h"
#include "migrate.h"

static bool is_matrix(const Graph *g)
{
    return g->be == matrix_backend();
}

/** Validates a (graph, out) pair and clears *out. */
static bool begin_build(const Graph *g, Graph **out)
{
    if (out != NULL) {
        *out = NULL;
    }
    if (g == NULL || out == NULL) {
        graph_set_error(GRAPH_ERR_INVALID_ARGUMENT);
        return false;
    }
    return true;
}

/** Turns the library's raw result into a Graph (or frees it on error). */
static GraphStatus finish_build(const Graph *like, const GraphBackend *be, void *raw,
                                GraphStatus status, Graph **out)
{
    if (status != GRAPH_OK) {
        return graph_set_error(status);
    }
    Graph *ng = graph_adopt(be, raw, like);
    if (ng == NULL) {
        return graph_last_error();
    }
    *out = ng;
    return GRAPH_OK;
}

#define UNARY_BUILD(NAME, OP, EXTRA_PARAMS, EXTRA_ARGS)                                   \
    GraphStatus graph_##NAME(Graph *g, EXTRA_PARAMS Graph **out)                          \
    {                                                                                     \
        if (!begin_build(g, out)) {                                                       \
            return GRAPH_ERR_INVALID_ARGUMENT;                                            \
        }                                                                                 \
        GraphStatus status;                                                               \
        void       *raw;                                                                  \
        if (is_matrix(g)) {                                                               \
            MatrixGraph *r = NULL;                                                        \
            status = matrix_##NAME((const MatrixGraph *)g->impl, EXTRA_ARGS &r);          \
            raw    = r;                                                                   \
        } else {                                                                          \
            ListGraph *r = NULL;                                                          \
            status = list_##NAME((const ListGraph *)g->impl, EXTRA_ARGS &r);              \
            raw    = r;                                                                   \
        }                                                                                 \
        graph_after_op(g, OP);                                                            \
        return finish_build(g, g->be, raw, status, out);                                  \
    }

#define NO_PARAMS
#define NO_ARGS
#define SUBGRAPH_PARAMS const int *vertices, size_t count,
#define SUBGRAPH_ARGS   vertices, count,

UNARY_BUILD(transpose, OP_ALGO_RUN, NO_PARAMS, NO_ARGS)
UNARY_BUILD(transitive_closure, OP_ALL_PAIRS_RUN, NO_PARAMS, NO_ARGS)
UNARY_BUILD(transitive_reduction, OP_ALL_PAIRS_RUN, NO_PARAMS, NO_ARGS)
UNARY_BUILD(induced_subgraph, OP_ALGO_RUN, SUBGRAPH_PARAMS, SUBGRAPH_ARGS)

GraphStatus graph_complement(Graph *g, Graph **out)
{
    if (!begin_build(g, out)) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    GraphStatus status;
    void       *raw;
    if (is_matrix(g)) {
        MatrixGraph *r = NULL;
        status = matrix_graph_complement((const MatrixGraph *)g->impl, &r);
        raw    = r;
    } else {
        ListGraph *r = NULL;
        status = list_graph_complement((const ListGraph *)g->impl, &r);
        raw    = r;
    }
    graph_after_op(g, OP_ALGO_RUN);
    return finish_build(g, g->be, raw, status, out);
}

Graph *graph_copy(Graph *g)
{
    if (g == NULL) {
        graph_set_error(GRAPH_ERR_INVALID_ARGUMENT);
        return NULL;
    }
    void *raw = is_matrix(g) ? (void *)matrix_graph_copy((const MatrixGraph *)g->impl)
                             : (void *)list_graph_copy((const ListGraph *)g->impl);
    Graph *ng = graph_adopt(g->be, raw, g);
    return ng;
}

/* ------------------------------------------------------------------ */
/*  Two-graph operations                                                */
/* ------------------------------------------------------------------ */

/**
 * The handle of @b in a's representation.  If @b already uses it, that is b's
 * own handle and *temp stays NULL; otherwise a converted copy is made, stored
 * in *temp (the caller must destroy it with a's backend) and returned.
 */
static void *in_rep_of(const Graph *a, const Graph *b, void **temp)
{
    *temp = NULL;
    if (a->be == b->be) {
        return b->impl;
    }
    *temp = graph_clone(b->be, b->impl, a->be);
    return *temp;
}

#define BINARY_BUILD(NAME, LIBNAME)                                                              \
    GraphStatus graph_##NAME(Graph *a, Graph *b, Graph **out)                             \
    {                                                                                     \
        if (!begin_build(a, out) || b == NULL) {                                          \
            return graph_set_error(GRAPH_ERR_INVALID_ARGUMENT);                           \
        }                                                                                 \
        void *temp;                                                                       \
        void *other = in_rep_of(a, b, &temp);                                             \
        if (other == NULL) {                                                              \
            return graph_set_error(GRAPH_ERR_OUT_OF_MEMORY);                              \
        }                                                                                 \
        GraphStatus status;                                                               \
        void       *raw;                                                                  \
        if (is_matrix(a)) {                                                               \
            MatrixGraph *r = NULL;                                                        \
            status = matrix_##LIBNAME((const MatrixGraph *)a->impl,                       \
                                   (const MatrixGraph *)other, &r);                       \
            raw    = r;                                                                   \
        } else {                                                                          \
            ListGraph *r = NULL;                                                          \
            status = list_##LIBNAME((const ListGraph *)a->impl, (const ListGraph *)other, &r);\
            raw    = r;                                                                   \
        }                                                                                 \
        if (temp != NULL) {                                                               \
            a->be->destroy(temp);                                                         \
        }                                                                                 \
        graph_after_op(a, OP_ALGO_RUN);                                                   \
        return finish_build(a, a->be, raw, status, out);                                  \
    }

BINARY_BUILD(union, graph_union)
BINARY_BUILD(intersection, graph_intersection)

GraphStatus graph_is_isomorphic(Graph *a, Graph *b, bool *result, int *mapping)
{
    if (a == NULL || b == NULL) {
        return graph_set_error(GRAPH_ERR_INVALID_ARGUMENT);
    }
    void *temp;
    void *other = in_rep_of(a, b, &temp);
    if (other == NULL) {
        return graph_set_error(GRAPH_ERR_OUT_OF_MEMORY);
    }
    const GraphStatus status =
        is_matrix(a) ? matrix_is_isomorphic((const MatrixGraph *)a->impl,
                                            (const MatrixGraph *)other, result, mapping)
                     : list_is_isomorphic((const ListGraph *)a->impl,
                                          (const ListGraph *)other, result, mapping);
    if (temp != NULL) {
        a->be->destroy(temp);
    }
    graph_after_op(a, OP_ADJACENCY_RUN);
    return status == GRAPH_OK ? GRAPH_OK : graph_set_error(status);
}

/* ------------------------------------------------------------------ */
/*  Generators                                                          */
/* ------------------------------------------------------------------ */

static Graph *adopt_generated(GraphRep rep, void *raw)
{
    if (raw == NULL) {
        graph_set_error(GRAPH_ERR_INVALID_ARGUMENT);   /* bad parameters or out of memory */
        return NULL;
    }
    return graph_adopt(graph_backend_for(rep), raw, NULL);
}

#define GENERATOR(NAME, PARAMS, ARGS)                                                     \
    Graph *graph_gen_##NAME(GraphRep rep, PARAMS)                                         \
    {                                                                                     \
        void *raw = rep == GRAPH_MATRIX ? (void *)matrix_gen_##NAME ARGS                  \
                                        : (void *)list_gen_##NAME ARGS;                   \
        return adopt_generated(rep, raw);                                                 \
    }

#define P(...) __VA_ARGS__
GENERATOR(erdos_renyi, P(size_t n, double p, unsigned long long seed, bool directed,
                         double wmin, double wmax), (n, p, seed, directed, wmin, wmax))
GENERATOR(barabasi_albert, P(size_t n, size_t m, unsigned long long seed), (n, m, seed))
GENERATOR(grid, P(size_t rows, size_t cols), (rows, cols))
GENERATOR(path, P(size_t n, bool directed), (n, directed))
GENERATOR(cycle, P(size_t n, bool directed), (n, directed))
GENERATOR(star, P(size_t n), (n))
GENERATOR(complete, P(size_t n, bool directed), (n, directed))
GENERATOR(random_dag, P(size_t n, double p, unsigned long long seed, double wmin, double wmax),
          (n, p, seed, wmin, wmax))
GENERATOR(random_tree, P(size_t n, unsigned long long seed), (n, seed))

/* ------------------------------------------------------------------ */
/*  Readers                                                             */
/* ------------------------------------------------------------------ */

GraphStatus graph_read_edge_list(FILE *in, bool directed, GraphRep rep, Graph **out)
{
    if (out == NULL) {
        return graph_set_error(GRAPH_ERR_INVALID_ARGUMENT);
    }
    *out = NULL;
    GraphStatus status;
    void       *raw;
    if (rep == GRAPH_MATRIX) {
        MatrixGraph *r = NULL;
        status = matrix_read_edge_list(in, directed, &r);
        raw    = r;
    } else {
        ListGraph *r = NULL;
        status = list_read_edge_list(in, directed, &r);
        raw    = r;
    }
    return finish_build(NULL, graph_backend_for(rep), raw, status, out);
}

GraphStatus graph_read_named_edge_list(FILE *in, bool directed, GraphRep rep, Graph **out,
                                       char ***names_out, size_t *count)
{
    if (out == NULL) {
        return graph_set_error(GRAPH_ERR_INVALID_ARGUMENT);
    }
    *out = NULL;
    GraphStatus status;
    void       *raw;
    if (rep == GRAPH_MATRIX) {
        MatrixGraph *r = NULL;
        status = matrix_read_named_edge_list(in, directed, &r, names_out, count);
        raw    = r;
    } else {
        ListGraph *r = NULL;
        status = list_read_named_edge_list(in, directed, &r, names_out, count);
        raw    = r;
    }
    return finish_build(NULL, graph_backend_for(rep), raw, status, out);
}

GraphStatus graph_load_binary(FILE *in, GraphRep rep, Graph **out)
{
    if (out == NULL) {
        return graph_set_error(GRAPH_ERR_INVALID_ARGUMENT);
    }
    *out = NULL;
    GraphStatus status;
    void       *raw;
    if (rep == GRAPH_MATRIX) {
        MatrixGraph *r = NULL;
        status = matrix_load_binary(in, &r);
        raw    = r;
    } else {
        ListGraph *r = NULL;
        status = list_load_binary(in, &r);
        raw    = r;
    }
    return finish_build(NULL, graph_backend_for(rep), raw, status, out);
}

void graph_names_free(char **names, size_t count)
{
    list_names_free(names, count);   /* identical to matrix_names_free */
}

/* ------------------------------------------------------------------ */
/*  Helpers that do not take a graph                                    */
/* ------------------------------------------------------------------ */

GraphStatus graph_hungarian(const double *cost, size_t rows, size_t cols,
                            int *assignment, double *total)
{
    const GraphStatus status = list_hungarian(cost, rows, cols, assignment, total);
    return status == GRAPH_OK ? GRAPH_OK : graph_set_error(status);
}

size_t graph_path_from_parent(const int *parent, size_t n, int src, int dst, int *path, size_t cap)
{
    return list_path_from_parent(parent, n, src, dst, path, cap);
}

size_t graph_path_from_next(const int *next, size_t n, int u, int v, int *path, size_t cap)
{
    return list_path_from_next(next, n, u, v, path, cap);
}

/* ------------------------------------------------------------------ */
/*  Advanced algorithms that cannot be generated (graph_advanced.h)     */
/* ------------------------------------------------------------------ */

GraphStatus graph_adv_subgraph_isomorphism(Graph *pattern, Graph *target, bool induced,
                                           bool match_weights, int *mapping, bool *found)
{
    if (pattern == NULL || target == NULL) {
        return graph_set_error(GRAPH_ERR_INVALID_ARGUMENT);
    }
    void *temp;
    void *other = in_rep_of(pattern, target, &temp);
    if (other == NULL) {
        return graph_set_error(GRAPH_ERR_OUT_OF_MEMORY);
    }
    const GraphStatus status =
        is_matrix(pattern)
            ? matrix_adv_subgraph_isomorphism((const MatrixGraph *)pattern->impl,
                                              (const MatrixGraph *)other, induced,
                                              match_weights, mapping, found)
            : list_adv_subgraph_isomorphism((const ListGraph *)pattern->impl,
                                            (const ListGraph *)other, induced,
                                            match_weights, mapping, found);
    if (temp != NULL) {
        pattern->be->destroy(temp);
    }
    graph_after_op(pattern, OP_ADJACENCY_RUN);
    return status == GRAPH_OK ? GRAPH_OK : graph_set_error(status);
}

GraphStatus graph_adv_graph_isomorphism(Graph *a, Graph *b, bool match_weights,
                                        bool *isomorphic, int *mapping)
{
    if (a == NULL || b == NULL) {
        return graph_set_error(GRAPH_ERR_INVALID_ARGUMENT);
    }
    void *temp;
    void *other = in_rep_of(a, b, &temp);
    if (other == NULL) {
        return graph_set_error(GRAPH_ERR_OUT_OF_MEMORY);
    }
    const GraphStatus status =
        is_matrix(a) ? matrix_adv_graph_isomorphism((const MatrixGraph *)a->impl,
                                                    (const MatrixGraph *)other,
                                                    match_weights, isomorphic, mapping)
                     : list_adv_graph_isomorphism((const ListGraph *)a->impl,
                                                  (const ListGraph *)other,
                                                  match_weights, isomorphic, mapping);
    if (temp != NULL) {
        a->be->destroy(temp);
    }
    graph_after_op(a, OP_ADJACENCY_RUN);
    return status == GRAPH_OK ? GRAPH_OK : graph_set_error(status);
}

GraphStatus graph_adv_apsp_path(const double *dist, const int *pred, size_t n,
                                int source, int target, GraphPath *out)
{
    /* Representation-independent, like graph_adv_path_from_parents. */
    return matrix_adv_apsp_path(dist, pred, n, source, target, out);
}

void graph_adv_cycles_free(GraphAdvCycles *cycles)
{
    if (cycles != NULL) {
        free(cycles->offsets);
        free(cycles->members);
        *cycles = (GraphAdvCycles){0};
    }
}

GraphStatus graph_adv_path_from_parents(const int *parent, size_t n, int source, int target,
                                        double cost, GraphPath *out)
{
    /* Representation-independent: the list and matrix versions are identical. */
    return matrix_adv_path_from_parents(parent, n, source, target, cost, out);
}

/*
 * Dynamic connectivity.  The library trackers borrow the raw ListGraph /
 * MatrixGraph, but the facade may migrate the graph to the other representation
 * (freeing the old one) between two calls.  So this wrapper remembers which
 * handle it was built on and rebuilds the inner tracker if the graph now lives
 * somewhere else.  (Rebuilding is O(V + E), the same cost the tracker pays after
 * any edit.)  Limitation: if the graph is converted away and back so that the
 * allocator hands out the very same address again, the change goes unnoticed
 * until the next edit; migrations are also counted to make that very unlikely.
 */
struct GraphAdvDynConn {
    Graph              *g;
    const GraphBackend *be;
    void               *impl;
    size_t              migrations;
    void               *inner;      /* MatrixAdvDynConn * or ListAdvDynConn * */
};

static void dc_drop_inner(GraphAdvDynConn *dc)
{
    if (dc->inner != NULL) {
        if (dc->be == matrix_backend()) {
            matrix_adv_dynconn_destroy((MatrixAdvDynConn *)dc->inner);
        } else {
            list_adv_dynconn_destroy((ListAdvDynConn *)dc->inner);
        }
        dc->inner = NULL;
    }
}

static GraphStatus dc_sync(GraphAdvDynConn *dc)
{
    Graph *g = dc->g;
    if (dc->inner != NULL && dc->be == g->be && dc->impl == g->impl &&
        dc->migrations == graph_migration_count(g)) {
        return GRAPH_OK;
    }
    dc_drop_inner(dc);
    dc->be         = g->be;
    dc->impl       = g->impl;
    dc->migrations = graph_migration_count(g);
    dc->inner = is_matrix(g) ? (void *)matrix_adv_dynconn_create((const MatrixGraph *)g->impl)
                             : (void *)list_adv_dynconn_create((const ListGraph *)g->impl);
    return dc->inner != NULL ? GRAPH_OK : GRAPH_ERR_OUT_OF_MEMORY;
}

GraphAdvDynConn *graph_adv_dynconn_create(Graph *g)
{
    if (g == NULL) {
        graph_set_error(GRAPH_ERR_INVALID_ARGUMENT);
        return NULL;
    }
    GraphAdvDynConn *dc = calloc(1, sizeof *dc);
    if (dc == NULL) {
        graph_set_error(GRAPH_ERR_OUT_OF_MEMORY);
        return NULL;
    }
    dc->g = g;
    const GraphStatus st = dc_sync(dc);
    if (st != GRAPH_OK) {
        graph_set_error(st);
        free(dc);
        return NULL;
    }
    return dc;
}

void graph_adv_dynconn_destroy(GraphAdvDynConn *dc)
{
    if (dc != NULL) {
        dc_drop_inner(dc);
        free(dc);
    }
}

GraphStatus graph_adv_dynconn_connected(GraphAdvDynConn *dc, int u, int v, bool *connected)
{
    if (dc == NULL) {
        return graph_set_error(GRAPH_ERR_INVALID_ARGUMENT);
    }
    GraphStatus st = dc_sync(dc);
    if (st == GRAPH_OK) {
        st = dc->be == matrix_backend()
                 ? matrix_adv_dynconn_connected((MatrixAdvDynConn *)dc->inner, u, v, connected)
                 : list_adv_dynconn_connected((ListAdvDynConn *)dc->inner, u, v, connected);
    }
    return st == GRAPH_OK ? GRAPH_OK : graph_set_error(st);
}

GraphStatus graph_adv_dynconn_components(GraphAdvDynConn *dc, size_t *count)
{
    if (dc == NULL) {
        return graph_set_error(GRAPH_ERR_INVALID_ARGUMENT);
    }
    GraphStatus st = dc_sync(dc);
    if (st == GRAPH_OK) {
        st = dc->be == matrix_backend()
                 ? matrix_adv_dynconn_components((MatrixAdvDynConn *)dc->inner, count)
                 : list_adv_dynconn_components((ListAdvDynConn *)dc->inner, count);
    }
    return st == GRAPH_OK ? GRAPH_OK : graph_set_error(st);
}
