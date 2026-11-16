/*
 * graph.c — the unified Graph facade.  See graph.h.
 *
 * Every public call does the same three things: forward to the active
 * backend through its vtable, tell the observer what happened, and (in
 * adaptive mode) give the observer a chance to migrate the graph.  Migration
 * only ever happens at the END of a call, so no pointer handed out earlier in
 * the same call can dangle.
 */

#include "graph.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "graph_internal.h"
#include "graph_status.h"
#include "adjacency_list.h"
#include "list_backend.h"
#include "memory_budget.h"
#include "matrix_backend.h"
#include "migrate.h"

/* ------------------------------------------------------------------ */
/*  Helpers                                                             */
/* ------------------------------------------------------------------ */

const GraphBackend *graph_backend_for(GraphRep rep)
{
    return rep == GRAPH_MATRIX ? matrix_backend() : list_backend();
}

/* Record a failure code, then return the given value. */
#define FAIL(code, value) do { graph_set_error(code); return (value); } while (0)

/** May this graph live in representation @b?  A representation with a size limit (the matrix: its
 *  presence bits take V^2/8 bytes) is ruled out once the graph has more vertices than the budget
 *  could ever hold in it, whatever its cost model says. */
static bool representation_allows(const Graph *g, const GraphBackend *b)
{
    if (b->max_vertices == NULL) {
        return true;
    }
    const bool weighted = !(g->be->is_unweighted != NULL && g->be->is_unweighted(g->impl));
    return g->be->id_bound(g->impl) <= b->max_vertices(observer_memory_budget(g->obs), weighted);
}

/** Record an operation and, in adaptive mode, let the observer migrate.
 *  A frozen graph is read-only and shared between threads, so it records nothing. */
void graph_after_op(Graph *g, OpKind op)
{
    if (g->frozen) {
        return;
    }
    observer_note(g->obs, op);
    if (!g->adaptive) {
        return;
    }
    /* The matrix is only a candidate while the graph is small enough for it to be possible. */
    const GraphBackend *candidates[2] = { list_backend(), NULL };
    size_t n_candidates = 1;
    if (representation_allows(g, matrix_backend())) {
        candidates[n_candidates++] = matrix_backend();
    }
    (void)observer_maybe_migrate(g->obs, &g->be, &g->impl, candidates, n_candidates, NULL);
}

/** Alive start vertex; otherwise records INVALID_ARGUMENT (NULL graph) or INVALID_VERTEX. */
static bool alive_start(const Graph *g, int v)
{
    if (g == NULL) {
        graph_set_error(GRAPH_ERR_INVALID_ARGUMENT);
        return false;
    }
    if (!g->be->is_alive(g->impl, v)) {
        graph_set_error(GRAPH_ERR_INVALID_VERTEX);
        return false;
    }
    return true;
}

/* ------------------------------------------------------------------ */
/*  Construction / destruction                                          */
/* ------------------------------------------------------------------ */

GraphConfig graph_default_config(bool directed)
{
    GraphConfig cfg = {
        .directed      = directed,
        .start         = GRAPH_LIST,
        .adaptive      = true,
        .observer      = observer_default_config(),
        .cache_entries = 32,
        .cache_bytes   = 4u * 1024u * 1024u,
        .unweighted    = false,
    };
    cfg.observer.memory_budget = graph_memory_budget();   /* 4 GiB unless changed at start-up */
    return cfg;
}

Graph *graph_create_with(const GraphConfig *config, size_t n)
{
    const GraphConfig cfg = config != NULL ? *config : graph_default_config(false);

    Graph *g = calloc(1, sizeof(*g));
    if (g == NULL) {
        return NULL;
    }
    g->be       = graph_backend_for(cfg.start);
    g->adaptive = cfg.adaptive;
    /* An unweighted graph stays unweighted across representations (graph_clone keeps the flavour). */
    g->impl     = cfg.unweighted && g->be->create_unweighted != NULL
                      ? g->be->create_unweighted(n, cfg.directed)
                      : g->be->create(n, cfg.directed);
    g->obs      = observer_create(&cfg.observer);
    if (cfg.cache_entries > 0) {
        g->cache = cache_create(cfg.cache_entries, cfg.cache_bytes);
    }
    if (g->impl == NULL || g->obs == NULL || (cfg.cache_entries > 0 && g->cache == NULL)) {
        graph_destroy(g);
        return NULL;
    }
    return g;
}

Graph *graph_adopt(const GraphBackend *be, void *impl, const Graph *like)
{
    if (impl == NULL) {
        graph_set_error(GRAPH_ERR_OUT_OF_MEMORY);   /* the producer failed to allocate */
        return NULL;
    }
    GraphConfig cfg = graph_default_config(be->is_directed(impl));
    if (like != NULL) {
        cfg.adaptive = like->adaptive;
    }
    Graph *g = calloc(1, sizeof(*g));
    if (g != NULL) {
        g->be       = be;
        g->impl     = impl;
        g->adaptive = cfg.adaptive;
        g->obs      = observer_create(&cfg.observer);
        g->cache    = cache_create(cfg.cache_entries, cfg.cache_bytes);
    }
    if (g == NULL || g->obs == NULL || g->cache == NULL) {
        if (g == NULL) {
            be->destroy(impl);
        } else {
            graph_destroy(g);   /* destroys impl too */
        }
        graph_set_error(GRAPH_ERR_OUT_OF_MEMORY);
        return NULL;
    }
    return g;
}

Graph *graph_create(bool directed)
{
    const GraphConfig cfg = graph_default_config(directed);
    return graph_create_with(&cfg, 0);
}

Graph *graph_create_sized(size_t n, bool directed)
{
    const GraphConfig cfg = graph_default_config(directed);
    return graph_create_with(&cfg, n);
}

void graph_destroy(Graph *g)
{
    if (g == NULL) {
        return;
    }
    if (g->impl != NULL) {
        g->be->destroy(g->impl);
    }
    observer_destroy(g->obs);
    cache_destroy(g->cache);
    free(g);
}

/* ------------------------------------------------------------------ */
/*  Representation control                                              */
/* ------------------------------------------------------------------ */

void graph_set_adaptive(Graph *g, bool adaptive)
{
    if (g != NULL) {
        g->adaptive = adaptive;
    }
}

bool graph_is_adaptive(const Graph *g)
{
    return g != NULL && g->adaptive;
}

bool graph_convert(Graph *g, GraphRep to)
{
    if (g == NULL) {
        FAIL(GRAPH_ERR_INVALID_ARGUMENT, false);
    }
    if (g->frozen) {
        FAIL(GRAPH_ERR_FROZEN, false);
    }
    const GraphBackend *target = graph_backend_for(to);
    if (target == g->be) {
        return true;
    }
    if (!representation_allows(g, target)) {
        FAIL(GRAPH_ERR_TOO_LARGE, false);          /* more vertices than that representation may hold */
    }
    /* Refuse up front if the other representation would not fit: during the copy both
     * exist at once, so the budget has to cover the old graph plus the new one. */
    const size_t budget = observer_memory_budget(g->obs);
    if (budget != 0 && target->estimate_memory != NULL) {
        const GraphShape shape = graph_shape_of(g->be, g->impl);
        if (g->be->memory_bytes(g->impl) + target->estimate_memory(&shape) > budget) {
            FAIL(GRAPH_ERR_TOO_LARGE, false);
        }
    }
    void *fresh = graph_clone(g->be, g->impl, target);
    if (fresh == NULL) {
        if (graph_last_error() != GRAPH_ERR_TOO_LARGE) {
            graph_set_error(GRAPH_ERR_OUT_OF_MEMORY);
        }
        return false;
    }
    g->be->destroy(g->impl);
    g->be   = target;
    g->impl = fresh;
    return true;
}

size_t graph_matrix_max_vertices(bool weighted)
{
    return matrix_backend()->max_vertices(graph_memory_budget(), weighted);
}

bool graph_shrink_to_fit(Graph *g)
{
    if (g == NULL) {
        FAIL(GRAPH_ERR_INVALID_ARGUMENT, false);
    }
    if (g->frozen) {
        FAIL(GRAPH_ERR_FROZEN, false);
    }
    if (g->be != list_backend()) {
        return true;                              /* the matrix has no per-vertex slack to release */
    }
    return list_graph_shrink_to_fit((ListGraph *)g->impl);
}

void graph_freeze(Graph *g)
{
    if (g != NULL) {
        g->frozen = true;
    }
}

void graph_thaw(Graph *g)
{
    if (g != NULL) {
        g->frozen = false;
    }
}

bool graph_is_frozen(const Graph *g)
{
    return g != NULL && g->frozen;
}

GraphRep graph_representation(const Graph *g)
{
    return (g != NULL && g->be == matrix_backend()) ? GRAPH_MATRIX : GRAPH_LIST;
}

const char *graph_rep_name(GraphRep rep)
{
    return rep == GRAPH_MATRIX ? "matrix" : "list";
}

size_t graph_migration_count(const Graph *g)
{
    return g == NULL ? 0 : observer_migrations(g->obs);
}

/* ------------------------------------------------------------------ */
/*  Queries                                                             */
/* ------------------------------------------------------------------ */

size_t graph_vertex_count(const Graph *g)
{
    return g == NULL ? 0 : g->be->vertex_count(g->impl);
}

size_t graph_id_bound(const Graph *g)
{
    return g == NULL ? 0 : g->be->id_bound(g->impl);
}

bool graph_cache_stats(const Graph *g, CacheStats *out)
{
    if (g == NULL || out == NULL || g->cache == NULL) {
        return false;
    }
    *out = cache_stats(g->cache);
    return true;
}

size_t graph_edge_count(const Graph *g)
{
    return g == NULL ? 0 : g->be->edge_count(g->impl);
}

bool graph_is_directed(const Graph *g)
{
    return g != NULL && g->be->is_directed(g->impl);
}

bool graph_is_alive(const Graph *g, int v)
{
    return g != NULL && g->be->is_alive(g->impl, v);
}

size_t graph_degree(const Graph *g, int v)
{
    return graph_is_alive(g, v) ? g->be->degree(g->impl, v) : 0;
}

bool graph_has_edge(Graph *g, int u, int v, double *weight_out)
{
    if (g == NULL) {
        FAIL(GRAPH_ERR_INVALID_ARGUMENT, false);
    }
    const bool found = g->be->has_edge(g->impl, u, v, weight_out);
    graph_after_op(g, OP_HAS_EDGE);
    return found;
}

size_t graph_neighbors(Graph *g, int v, int *neighbors, double *weights, size_t cap)
{
    if (!alive_start(g, v)) {
        return 0;
    }
    size_t   count = 0;
    GraphIter it;
    int      nb;
    double   w;
    g->be->iter_begin(g->impl, v, &it);
    while (g->be->iter_next(&it, &nb, weights != NULL ? &w : NULL)) {
        if (count < cap) {
            if (neighbors != NULL) {
                neighbors[count] = nb;
            }
            if (weights != NULL) {
                weights[count] = w;
            }
        }
        ++count;
    }
    graph_after_op(g, OP_NEIGHBOR_SCAN);
    return count;
}

/* ------------------------------------------------------------------ */
/*  Mutation                                                            */
/* ------------------------------------------------------------------ */

int graph_add_vertex(Graph *g)
{
    if (g == NULL) {
        FAIL(GRAPH_ERR_INVALID_ARGUMENT, -1);
    }
    if (g->frozen) {
        FAIL(GRAPH_ERR_FROZEN, -1);
    }
    const int id = g->be->add_vertex(g->impl);
    graph_after_op(g, OP_ADD_VERTEX);
    return id;
}

bool graph_remove_vertex(Graph *g, int v)
{
    if (g == NULL) {
        FAIL(GRAPH_ERR_INVALID_ARGUMENT, false);
    }
    if (g->frozen) {
        FAIL(GRAPH_ERR_FROZEN, false);
    }
    const bool ok = g->be->remove_vertex(g->impl, v);
    graph_after_op(g, OP_REMOVE_VERTEX);
    return ok;
}

bool graph_add_edge(Graph *g, int u, int v, double w)
{
    if (g == NULL) {
        FAIL(GRAPH_ERR_INVALID_ARGUMENT, false);
    }
    if (g->frozen) {
        FAIL(GRAPH_ERR_FROZEN, false);
    }
    const bool ok = g->be->add_edge(g->impl, u, v, w);
    graph_after_op(g, OP_ADD_EDGE);
    return ok;
}

bool graph_remove_edge(Graph *g, int u, int v)
{
    if (g == NULL) {
        FAIL(GRAPH_ERR_INVALID_ARGUMENT, false);
    }
    if (g->frozen) {
        FAIL(GRAPH_ERR_FROZEN, false);
    }
    const bool ok = g->be->remove_edge(g->impl, u, v);
    graph_after_op(g, OP_REMOVE_EDGE);
    return ok;
}

bool graph_needs_compaction(const Graph *g)
{
    return g != NULL && g->be->needs_compaction != NULL && g->be->needs_compaction(g->impl);
}

bool graph_compact(Graph *g, int *remap)
{
    if (g == NULL) {
        FAIL(GRAPH_ERR_INVALID_ARGUMENT, false);
    }
    if (g->frozen) {
        FAIL(GRAPH_ERR_FROZEN, false);
    }
    if (g->be->compact == NULL) {
        FAIL(GRAPH_ERR_NO_BACKEND_SUPPORT, false);
    }
    return g->be->compact(g->impl, remap);
}

/* ------------------------------------------------------------------ */
/*  Result cache                                                        */
/*                                                                      */
/*  Every cached result is a Header followed by up to two arrays whose  */
/*  lengths depend only on the graph (id_bound / vertex_count), so the   */
/*  caller always knows how many bytes to expect.  Keys are validated    */
/*  against the backend's version counter inside the cache itself.       */
/* ------------------------------------------------------------------ */

typedef struct {
    size_t ret;     /* integer result: visit count, component count, flag */
    double dret;    /* floating result: MST weight                        */
} Header;

static bool cache_fetch(Graph *g, CacheKey key, Header *hdr,
                        void *p1, size_t n1, void *p2, size_t n2)
{
    if (g->cache == NULL || g->frozen) {
        return false;     /* a lookup reorders the LRU list: not allowed on a shared handle */
    }
    size_t size = 0;
    const unsigned char *blob =
        cache_lookup(g->cache, key, g->be->version(g->impl), &size);
    if (blob == NULL || size != sizeof(Header) + n1 + n2) {
        observer_note(g->obs, OP_CACHE_MISS);
        return false;
    }
    memcpy(hdr, blob, sizeof(Header));
    if (n1 > 0 && p1 != NULL) {
        memcpy(p1, blob + sizeof(Header), n1);
    }
    if (n2 > 0 && p2 != NULL) {
        memcpy(p2, blob + sizeof(Header) + n1, n2);
    }
    observer_note(g->obs, OP_CACHE_HIT);
    return true;
}

static void cache_save(Graph *g, CacheKey key, const Header *hdr,
                       const void *p1, size_t n1, const void *p2, size_t n2)
{
    if (g->cache == NULL || g->frozen) {
        return;
    }
    unsigned char *blob = malloc(sizeof(Header) + n1 + n2);
    if (blob == NULL) {
        return;
    }
    memcpy(blob, hdr, sizeof(Header));
    if (n1 > 0) {
        memcpy(blob + sizeof(Header), p1, n1);
    }
    if (n2 > 0) {
        memcpy(blob + sizeof(Header) + n1, p2, n2);
    }
    (void)cache_store(g->cache, key, g->be->version(g->impl), blob, sizeof(Header) + n1 + n2);
    free(blob);
}

/* ------------------------------------------------------------------ */
/*  Algorithms                                                          */
/* ------------------------------------------------------------------ */

typedef size_t (*TraversalFn)(const void *, int, int *, int *);

static size_t traverse(Graph *g, CacheAlgo algo, TraversalFn fn,
                       int start, int *order, int *parent)
{
    if (g == NULL) {
        FAIL(GRAPH_ERR_INVALID_ARGUMENT, 0);
    }
    if (fn == NULL) {
        FAIL(GRAPH_ERR_NO_BACKEND_SUPPORT, 0);
    }
    if (!alive_start(g, start)) {
        return 0;
    }
    const size_t vc = g->be->vertex_count(g->impl);
    const size_t ib = g->be->id_bound(g->impl);
    const CacheKey key = { (uint32_t)algo, start, 0 };


    int *tmp_order  = calloc(vc > 0 ? vc : 1, sizeof(int));
    int *tmp_parent = calloc(ib > 0 ? ib : 1, sizeof(int));
    if (tmp_order == NULL || tmp_parent == NULL) {
        /* No room for the staging buffers: run directly, uncached. */
        free(tmp_order);
        free(tmp_parent);
        const size_t n = fn(g->impl, start, order, parent);
        graph_after_op(g, OP_ALGO_RUN);
        return n;
    }

    Header hdr;
    bool   ran = false;
    if (!cache_fetch(g, key, &hdr, tmp_order, vc * sizeof(int),
                     tmp_parent, ib * sizeof(int))) {
        hdr = (Header){ fn(g->impl, start, tmp_order, tmp_parent), 0.0 };
        cache_save(g, key, &hdr, tmp_order, vc * sizeof(int),
                   tmp_parent, ib * sizeof(int));
        ran = true;
    }
    if (order != NULL) {
        memcpy(order, tmp_order, vc * sizeof(int));
    }
    if (parent != NULL) {
        memcpy(parent, tmp_parent, ib * sizeof(int));
    }
    free(tmp_order);
    free(tmp_parent);
    if (ran) {
        graph_after_op(g, OP_ALGO_RUN);
    }
    return hdr.ret;
}

size_t graph_bfs(Graph *g, int start, int *order, int *parent)
{
    if (g == NULL) {
        FAIL(GRAPH_ERR_INVALID_ARGUMENT, 0);
    }
    return traverse(g, CACHE_BFS, g->be->bfs, start, order, parent);
}

size_t graph_dfs(Graph *g, int start, int *order, int *parent)
{
    if (g == NULL) {
        FAIL(GRAPH_ERR_INVALID_ARGUMENT, 0);
    }
    return traverse(g, CACHE_DFS, g->be->dfs, start, order, parent);
}

size_t graph_connected_components(Graph *g, int *component)
{
    if (g == NULL) {
        FAIL(GRAPH_ERR_INVALID_ARGUMENT, 0);
    }
    if (g->be->connected_components == NULL) {
        FAIL(GRAPH_ERR_NO_BACKEND_SUPPORT, 0);
    }
    const size_t ib = g->be->id_bound(g->impl);
    const CacheKey key = { CACHE_COMPONENTS, 0, 0 };
    Header hdr;

    int *tmp = calloc(ib > 0 ? ib : 1, sizeof(int));
    if (tmp == NULL) {
        const size_t c = g->be->connected_components(g->impl, component);
        graph_after_op(g, OP_ALGO_RUN);
        return c;
    }
    size_t count;
    if (cache_fetch(g, key, &hdr, tmp, ib * sizeof(int), NULL, 0)) {
        count = hdr.ret;
    } else {
        count = g->be->connected_components(g->impl, tmp);
        hdr = (Header){ count, 0.0 };
        cache_save(g, key, &hdr, tmp, ib * sizeof(int), NULL, 0);
        if (component != NULL) {
            memcpy(component, tmp, ib * sizeof(int));
        }
        free(tmp);
        graph_after_op(g, OP_ALGO_RUN);
        return count;
    }
    if (component != NULL) {
        memcpy(component, tmp, ib * sizeof(int));
    }
    free(tmp);
    return count;
}

bool graph_has_cycle(Graph *g)
{
    if (g == NULL) {
        FAIL(GRAPH_ERR_INVALID_ARGUMENT, false);
    }
    if (g->be->has_cycle == NULL) {
        FAIL(GRAPH_ERR_NO_BACKEND_SUPPORT, false);
    }
    const CacheKey key = { CACHE_CYCLE, 0, 0 };
    Header hdr;
    if (cache_fetch(g, key, &hdr, NULL, 0, NULL, 0)) {
        return hdr.ret != 0;
    }
    const bool cyclic = g->be->has_cycle(g->impl);
    hdr = (Header){ cyclic ? 1u : 0u, 0.0 };
    cache_save(g, key, &hdr, NULL, 0, NULL, 0);
    graph_after_op(g, OP_ALGO_RUN);
    return cyclic;
}

bool graph_topological_sort(Graph *g, int *order)
{
    if (g == NULL) {
        FAIL(GRAPH_ERR_INVALID_ARGUMENT, false);
    }
    if (g->be->topological_sort == NULL) {
        FAIL(GRAPH_ERR_NO_BACKEND_SUPPORT, false);
    }
    const size_t vc = g->be->vertex_count(g->impl);
    const CacheKey key = { CACHE_TOPO, 0, 0 };
    Header hdr;

    int *tmp = calloc(vc > 0 ? vc : 1, sizeof(int));
    if (tmp == NULL) {
        const bool ok = g->be->topological_sort(g->impl, order);
        graph_after_op(g, OP_ALGO_RUN);
        return ok;
    }
    bool ok;
    bool ran = false;
    if (cache_fetch(g, key, &hdr, tmp, vc * sizeof(int), NULL, 0)) {
        ok = hdr.ret != 0;
        if (!ok) {   /* the cached answer was a failure; report why again */
            graph_set_error(g->be->is_directed(g->impl) ? GRAPH_ERR_CYCLE : GRAPH_ERR_UNSUPPORTED);
        }
    } else {
        ok  = g->be->topological_sort(g->impl, tmp);
        ran = true;
        hdr = (Header){ ok ? 1u : 0u, 0.0 };
        cache_save(g, key, &hdr, tmp, vc * sizeof(int), NULL, 0);
    }
    if (ok && order != NULL) {
        memcpy(order, tmp, vc * sizeof(int));
    }
    free(tmp);
    if (ran) {
        graph_after_op(g, OP_ALGO_RUN);
    }
    return ok;
}

bool graph_dijkstra(Graph *g, int source, double *distance, int *parent)
{
    if (g == NULL) {
        FAIL(GRAPH_ERR_INVALID_ARGUMENT, false);
    }
    if (g->be->dijkstra == NULL) {
        FAIL(GRAPH_ERR_NO_BACKEND_SUPPORT, false);
    }
    if (distance == NULL) {
        FAIL(GRAPH_ERR_INVALID_ARGUMENT, false);
    }
    if (!alive_start(g, source)) {
        return false;
    }
    const size_t ib = g->be->id_bound(g->impl);
    const CacheKey key = { CACHE_DIJKSTRA, source, 0 };
    Header hdr;

    double *td = calloc(ib > 0 ? ib : 1, sizeof(double));
    int    *tp = calloc(ib > 0 ? ib : 1, sizeof(int));
    if (td == NULL || tp == NULL) {
        free(td);
        free(tp);
        const bool ok = g->be->dijkstra(g->impl, source, distance, parent);
        graph_after_op(g, OP_WEIGHTED_RUN);
        return ok;
    }
    bool ok;
    bool ran = false;
    if (cache_fetch(g, key, &hdr, td, ib * sizeof(double), tp, ib * sizeof(int))) {
        ok = true;      /* only successful runs are stored */
    } else {
        ok  = g->be->dijkstra(g->impl, source, td, tp);
        ran = true;
        if (ok) {
            hdr = (Header){ 1u, 0.0 };
            cache_save(g, key, &hdr, td, ib * sizeof(double), tp, ib * sizeof(int));
        }
    }
    if (ok) {
        memcpy(distance, td, ib * sizeof(double));
        if (parent != NULL) {
            memcpy(parent, tp, ib * sizeof(int));
        }
    }
    free(td);
    free(tp);
    if (ran) {
        graph_after_op(g, OP_WEIGHTED_RUN);
    }
    return ok;
}

double graph_mst_kruskal(Graph *g, GraphEdge *edges, size_t *edge_count)
{
    if (g == NULL) {
        FAIL(GRAPH_ERR_INVALID_ARGUMENT, INFINITY);
    }
    if (g->be->mst_kruskal == NULL) {
        FAIL(GRAPH_ERR_NO_BACKEND_SUPPORT, INFINITY);
    }
    const size_t vc    = g->be->vertex_count(g->impl);
    const size_t slots = vc > 1 ? vc - 1 : 1;
    const CacheKey key = { CACHE_MST, 0, 0 };
    Header hdr;

    GraphEdge *tmp = calloc(slots, sizeof(GraphEdge));
    size_t     tmp_count = 0;
    if (tmp == NULL) {
        size_t local = 0;
        const double w = g->be->mst_kruskal(g->impl, edges, edge_count != NULL ? edge_count : &local);
        graph_after_op(g, OP_WEIGHTED_RUN);
        return w;
    }
    double total;
    bool   ran = false;
    if (cache_fetch(g, key, &hdr, tmp, slots * sizeof(GraphEdge), NULL, 0)) {
        total     = hdr.dret;
        tmp_count = hdr.ret;
        if (!isfinite(total)) {   /* cached failure: report why again */
            graph_set_error(g->be->is_directed(g->impl) ? GRAPH_ERR_UNSUPPORTED : GRAPH_ERR_NOT_FOUND);
        }
    } else {
        total = g->be->mst_kruskal(g->impl, tmp, &tmp_count);
        ran   = true;
        hdr   = (Header){ tmp_count, total };
        cache_save(g, key, &hdr, tmp, slots * sizeof(GraphEdge), NULL, 0);
    }
    if (edges != NULL) {
        memcpy(edges, tmp, tmp_count * sizeof(GraphEdge));
    }
    if (edge_count != NULL) {
        *edge_count = tmp_count;
    }
    free(tmp);
    if (ran) {
        graph_after_op(g, OP_WEIGHTED_RUN);
    }
    return total;
}
