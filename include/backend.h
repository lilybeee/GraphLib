/*
 * backend.h — the interface every graph representation implements.
 *
 * The facade (graph.c, on main) talks to a representation only through a
 * GraphBackend vtable, which is what lets the observer migrate a live graph
 * from one representation to another.
 *
 * The vtable has three tiers:
 *
 *   Required   topology, mutation, and a neighbour iterator.  Everything the
 *              generic algorithms need, so any backend that fills these in is
 *              already complete and correct.
 *
 *   Optional   specialised algorithm kernels (bfs, dijkstra, ...) and
 *              add_new_edge.  A slot left NULL means "use the generic
 *              version".  A backend fills a kernel slot only when it can beat
 *              the generic version by exploiting its layout, e.g. a bitset-row
 *              BFS on a matrix or a heap Dijkstra on a list.
 *
 *   Cost model estimate_cost / estimate_memory / estimate_build_cost.  Each
 *              backend knows how fast it is, so each one reports its own
 *              numbers and the observer (policy.h) just compares them.  A
 *              backend without all three is never chosen as a migration target.
 *
 * ID conventions, array-size rules and complexity notes match
 * adjacency_list.h: IDs live in [0, id_bound), some of them dead.  Migration
 * preserves IDs, so a handle held by the caller survives a change of
 * representation.
 *
 * Every function receives the backend's own handle as 'impl'.
 */

#ifndef BACKEND_H
#define BACKEND_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "graph_types.h"
#include "stats.h"

/*
 * Neighbour iterator.  The caller owns the struct and treats it as opaque;
 * the backend decides what the fields mean (a list keeps array pointers, a
 * matrix keeps a column cursor).  Advancing it must not allocate.
 */
typedef struct GraphIter {
    const void *impl;
    int         vertex;
    size_t      index;
    size_t      end;
    const void *a;
    const void *b;
} GraphIter;

typedef struct GraphBackend {
    const char *name;               /* "list", "matrix", ... */

    /* ---- lifecycle ---- */
    void  *(*create)(size_t vertices, bool directed);
    void   (*destroy)(void *impl);

    /* ---- optional: weight-free graphs and a size limit (NULL = not supported / none) ---- */
    void  *(*create_unweighted)(size_t vertices, bool directed);
    bool   (*is_unweighted)(const void *impl);
    /* The most vertices this representation may hold under @budget bytes, because beyond that
     * it could not be filled within the budget however few edges it starts with (the matrix:
     * V^2/8 bytes of presence bits).  NULL means no such limit (the list grows with its edges). */
    size_t (*max_vertices)(size_t budget, bool weighted);

    /* ---- topology (all O(1) unless noted) ---- */
    size_t   (*vertex_count)(const void *impl);
    size_t   (*id_bound)(const void *impl);
    size_t   (*edge_count)(const void *impl);
    bool     (*is_directed)(const void *impl);
    bool     (*is_alive)(const void *impl, int v);
    size_t   (*degree)(const void *impl, int v);
    uint64_t (*version)(const void *impl);       /* bumps on every mutation */
    size_t   (*memory_bytes)(const void *impl);  /* O(V) */

    /* ---- mutation ---- */
    int   (*add_vertex)(void *impl);             /* new ID or -1 */
    bool  (*remove_vertex)(void *impl, int v);
    bool  (*add_edge)(void *impl, int u, int v, double weight);
    bool  (*remove_edge)(void *impl, int u, int v);
    bool  (*has_edge)(const void *impl, int u, int v, double *weight_out);

    /* ---- compaction (renumbers IDs; see adjacency_list.h) ---- */
    bool  (*needs_compaction)(const void *impl);
    bool  (*compact)(void *impl, int *remap);

    /* ---- neighbour iteration (required) ---- */
    /* iter_next stores the next neighbour and its weight; 'weight' may be
     * NULL when the caller does not need it, which lets a backend skip the
     * weight read (a cache miss on a matrix). */
    void  (*iter_begin)(const void *impl, int v, GraphIter *it);
    bool  (*iter_next)(GraphIter *it, int *neighbor, double *weight);

    /* ---- optional ---- */
    /* add_edge for an edge the caller knows is absent: skips the duplicate
     * check.  Used by migration to rebuild a graph without O(degree) scans. */
    bool   (*add_new_edge)(void *impl, int u, int v, double weight);

    size_t (*bfs)(const void *impl, int start, int *order, int *parent);
    size_t (*dfs)(const void *impl, int start, int *order, int *parent);
    size_t (*connected_components)(const void *impl, int *component);
    bool   (*has_cycle)(const void *impl);
    bool   (*topological_sort)(const void *impl, int *order);
    bool   (*dijkstra)(const void *impl, int source, double *distance, int *parent);
    double (*mst_kruskal)(const void *impl, GraphEdge *edges, size_t *edge_count);

    /* ---- cost model (used by the observer) ---- */
    /* Predicted cost, in nanoseconds, of running 'workload' on a graph of the
     * given shape stored in this representation. */
    double (*estimate_cost)(const GraphShape *shape, const Workload *workload);
    /* Predicted heap footprint of a graph of this shape. */
    size_t (*estimate_memory)(const GraphShape *shape);
    /* Predicted cost, in nanoseconds, of building this representation from an
     * existing graph of the given shape (a migration). */
    double (*estimate_build_cost)(const GraphShape *shape);
} GraphBackend;

#endif /* BACKEND_H */
