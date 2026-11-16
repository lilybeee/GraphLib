/*
 * graph.h — the unified graph API.
 *
 * This is the only header a user needs.  A Graph stores its edges in ONE
 * representation at a time, an adjacency list or an adjacency matrix, and
 * every function below behaves the same whichever is active.
 *
 * Choosing a representation
 * -------------------------
 *   - You may START on either representation (GraphConfig.start).
 *   - Adaptive mode (on by default) lets the library watch the workload and
 *     migrate the graph to whichever representation is predicted to be
 *     cheaper.  Opt out with graph_set_adaptive(g, false) or by creating the
 *     graph with config.adaptive = false; the graph then stays where it is.
 *   - graph_convert() switches representation by hand, in either mode.  In
 *     adaptive mode the observer is still free to move the graph again later;
 *     turn adaptation off first if you want the choice to stick.
 *   - graph_compare() is an explicit benchmark: it runs one algorithm on a
 *     copy of the graph in each representation and reports time and memory.
 *     It never modifies the graph and never triggers a migration.
 *
 * Memory
 * ------
 *   - Every graph may use at most the library's memory budget (memory_budget.h, 4 GiB by
 *     default).  Going over it makes the call fail with GRAPH_ERR_TOO_LARGE and leaves the graph
 *     as it was, so a request that is too big is an error you can handle, not a crash.
 *   - The matrix needs V^2/8 bytes of presence bits, so it has a size limit
 *     (graph_matrix_max_vertices); a larger graph can only be a list.  Adaptive mode never
 *     offers the matrix beyond that limit, and graph_convert refuses it.
 *   - GraphConfig.unweighted stores no weights (every edge is 1.0): about half the list's memory
 *     and no weight tiles in the matrix.  graph_shrink_to_fit() releases the growth slack of a
 *     finished graph.
 *
 * Vertex IDs
 * ----------
 * IDs live in [0, graph_id_bound()).  An ID in that range is alive or dead
 * (removed, waiting to be reused by graph_add_vertex).  IDs never change on
 * their own, not even across a migration.  The one operation that renumbers
 * them is graph_compact(), which is explicit and reports old -> new IDs.
 *
 * Array sizes for the algorithms
 * ------------------------------
 *   - arrays indexed by vertex ID (parent, component, distance) need
 *     graph_id_bound() entries; dead IDs receive -1 (INFINITY for distances);
 *   - arrays listing vertices (order) need graph_vertex_count() entries, and
 *     the MST edge array needs vertex_count - 1 entries.
 *
 * Thread safety: none.  Even "read" calls update workload statistics and may
 * migrate the graph, so callers must serialise access to a Graph.  graph_freeze() makes a
 * read-only handle that several threads may share.
 *
 * Where to go next: graph_algorithms.h (about 90 algorithms), graph_advanced.h (graph_adv_*),
 * graph_social.h (per-user queries), graph_stream.h (edge-list files too big to load) and
 * graph_viz.h (pictures).  The README has runnable examples of each.
 *
 * Link with -lm.
 */

#ifndef GRAPH_H
#define GRAPH_H

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>

#include "cache.h"          /* CacheStats */
#include "graph_results.h"  /* GraphGroups, GraphPath, ... */
#include "graph_status.h"   /* GraphStatus, graph_last_error() */
#include "graph_types.h"   /* GraphEdge */
#include "memory_budget.h"
#include "policy.h"        /* ObserverConfig */

typedef struct Graph Graph;

typedef enum {
    GRAPH_LIST,
    GRAPH_MATRIX
} GraphRep;

/* ------------------------------------------------------------------ */
/*  Construction / destruction                                          */
/* ------------------------------------------------------------------ */

typedef struct {
    bool           directed;
    GraphRep       start;         /* representation the graph begins in      */
    bool           adaptive;      /* allow automatic migration               */
    ObserverConfig observer;      /* thresholds for adaptive mode            */
    size_t         cache_entries; /* cached algorithm results; 0 = no cache  */
    size_t         cache_bytes;   /* payload byte limit for the cache        */
    bool           unweighted;    /* store no weights (every edge is 1.0): about half the memory
                                     in the list, and no 32 KB weight tiles in the matrix; any
                                     other edge weight is refused                         */
} GraphConfig;

/** directed as given; start = GRAPH_LIST; adaptive on; 32 cached results, 4 MiB. */
GraphConfig graph_default_config(bool directed);

/** An empty graph with default settings.  Returns NULL on allocation failure. */
Graph *graph_create(bool directed);

/** A graph with default settings and vertices 0 .. n-1 already alive. */
Graph *graph_create_sized(size_t n, bool directed);

/**
 * graph_create_with - full control.  @config may be NULL for the defaults
 * (undirected).  @n initial vertices.  Returns NULL on failure.
 */
Graph *graph_create_with(const GraphConfig *config, size_t n);

/** Frees everything.  Safe with NULL. */
void graph_destroy(Graph *g);

/* ------------------------------------------------------------------ */
/*  Representation control                                              */
/* ------------------------------------------------------------------ */

/** Opt in (true) or out (false) of automatic migration. */
void graph_set_adaptive(Graph *g, bool adaptive);
bool graph_is_adaptive(const Graph *g);

/**
 * graph_convert - move the graph to representation @to now.  A no-op that
 * returns true if it is already there.  On allocation failure the graph is
 * unchanged and false is returned.  If the other representation would not fit
 * in the memory budget (see memory_budget.h; the old and the new graph exist
 * together during the copy), nothing is copied and graph_last_error() is
 * GRAPH_ERR_TOO_LARGE.  Set the budget to 0 to convert regardless.
 * Complexity: O(V + E).
 */
bool graph_convert(Graph *g, GraphRep to);

/**
 * graph_shrink_to_fit - call once a graph is finished being built to give back the
 * growth slack its adjacency lists keep (typically 20-30% of the graph's memory).
 * Contents are unchanged; later edits still work but regrow the lists.  A no-op on
 * the matrix representation.  Returns false (graph intact) if frozen or out of memory.
 */
bool graph_shrink_to_fit(Graph *g);

/**
 * graph_matrix_max_vertices - the most vertices a graph may have and still be stored as a matrix
 * under the current memory budget (see memory_budget.h and matrix_max_vertices).  Beyond this,
 * graph_convert(g, GRAPH_MATRIX) fails with GRAPH_ERR_TOO_LARGE and adaptive mode never chooses
 * the matrix.  Weighted graphs get a smaller limit than unweighted ones.
 */
size_t graph_matrix_max_vertices(bool weighted);

/** The representation currently in use (for inspection; never needed to operate). */
GraphRep    graph_representation(const Graph *g);
const char *graph_rep_name(GraphRep rep);

/**
 * graph_freeze - make the handle read-only and safe to share between threads.
 *
 * A normal Graph is NOT thread-safe, even for "read" calls: every call updates
 * workload counters, may refresh the result cache and may migrate the graph.
 * A frozen graph does none of that, so any number of threads may call the
 * query and algorithm functions on it at the same time (each call works on
 * its own scratch memory; graph_last_error() is per thread).  While frozen,
 * every function that would modify the graph (add/remove vertex or edge,
 * compact, convert) fails with GRAPH_ERR_FROZEN, and results are not cached.
 *
 * You must make sure no other thread is inside a call on @g when you freeze
 * or thaw it (join the threads, or use your own lock around the switch).
 */
void graph_freeze(Graph *g);
void graph_thaw(Graph *g);
bool graph_is_frozen(const Graph *g);

/** How many automatic migrations have happened. */
size_t graph_migration_count(const Graph *g);

/* ------------------------------------------------------------------ */
/*  Queries                                                             */
/* ------------------------------------------------------------------ */

size_t graph_vertex_count(const Graph *g);   /* alive vertices                  */
size_t graph_id_bound(const Graph *g);       /* one past the largest ID in use  */
size_t graph_edge_count(const Graph *g);     /* undirected edges counted once   */

/**
 * graph_cache_stats - copy the result cache's counters (hits, misses, stores,
 * evictions, invalidations, entries, bytes) into @out.  Returns false, leaving
 * @out untouched, if @g or @out is NULL or the graph has no cache.  @out is a
 * plain copy: it stays valid after later calls and owns nothing.
 */
bool graph_cache_stats(const Graph *g, CacheStats *out);
bool   graph_is_directed(const Graph *g);
bool   graph_is_alive(const Graph *g, int v);
size_t graph_degree(const Graph *g, int v);  /* out-degree; 0 if v is not alive */

/** True iff edge (u -> v) exists; fills *weight_out if non-NULL. */
bool graph_has_edge(Graph *g, int u, int v, double *weight_out);

/**
 * graph_neighbors - copy v's out-neighbours (and weights, if @weights is
 * non-NULL) into caller arrays of capacity @cap.  Returns the number of
 * neighbours v has, which can exceed @cap; only min(result, cap) were written.
 * Size the arrays with graph_degree().  Order is unspecified.
 */
size_t graph_neighbors(Graph *g, int v, int *neighbors, double *weights, size_t cap);

/* ------------------------------------------------------------------ */
/*  Mutation                                                            */
/* ------------------------------------------------------------------ */

/** Add a vertex (reusing a dead ID if there is one).  Returns its ID, or -1. */
int  graph_add_vertex(Graph *g);

/** Remove v and every edge touching it.  False if v is not alive. */
bool graph_remove_vertex(Graph *g, int v);

/**
 * graph_add_edge - insert (u -> v) with weight w (and v -> u if undirected).
 * False if u or v is not alive, w is not finite, or the edge exists.
 */
bool graph_add_edge(Graph *g, int u, int v, double w);

/** Remove (u -> v).  False if it does not exist. */
bool graph_remove_edge(Graph *g, int u, int v);

/** True once more than a quarter of the ID space is dead.  Advisory only. */
bool graph_needs_compaction(const Graph *g);

/**
 * graph_compact - renumber alive vertices to 0 .. vertex_count-1, keeping
 * their relative order.  @remap (optional) needs graph_id_bound() entries
 * measured BEFORE the call; remap[old] = new ID, or -1 if old was dead.
 * Every vertex ID the caller holds is invalid afterwards unless translated.
 * False (graph unchanged) on failure.
 */
bool graph_compact(Graph *g, int *remap);

/* ------------------------------------------------------------------ */
/*  Algorithms                                                          */
/* ------------------------------------------------------------------ */

/**
 * graph_bfs / graph_dfs - traverse from @start.  @order receives visit order,
 * @parent (optional) the discovery tree (-1 for the root and unreached IDs).
 * Returns the number of vertices visited; 0 if @start is not alive.
 */
size_t graph_bfs(Graph *g, int start, int *order, int *parent);
size_t graph_dfs(Graph *g, int start, int *order, int *parent);

/**
 * graph_connected_components - label each alive vertex with a 0-based
 * component number (weak components if directed).  Returns the count.
 */
size_t graph_connected_components(Graph *g, int *component);

/** Cycle test (directed: a true cycle; undirected: any non-tree edge). */
bool graph_has_cycle(Graph *g);

/** Topological order of a DAG.  False if undirected or cyclic. */
bool graph_topological_sort(Graph *g, int *order);

/**
 * graph_dijkstra - single-source shortest paths.  False if @source is not
 * alive or any edge weight is negative.  @parent is optional.
 */
bool graph_dijkstra(Graph *g, int source, double *distance, int *parent);

/**
 * graph_mst_kruskal - minimum spanning tree of an undirected graph.  Returns
 * the total weight, or INFINITY if the graph is directed or disconnected.
 */
double graph_mst_kruskal(Graph *g, GraphEdge *edges, size_t *edge_count);

/* ------------------------------------------------------------------ */
/*  Benchmark / comparison                                              */
/* ------------------------------------------------------------------ */

typedef enum {
    GRAPH_ALGO_BFS,
    GRAPH_ALGO_DFS,
    GRAPH_ALGO_COMPONENTS,
    GRAPH_ALGO_CYCLE,
    GRAPH_ALGO_TOPO,
    GRAPH_ALGO_DIJKSTRA,
    GRAPH_ALGO_MST,
    GRAPH_ALGO_PAGERANK,       /* the library algorithms: see graph_algorithms.h  */
    GRAPH_ALGO_TRIANGLES,      /* (undirected only: "no result" on directed graphs) */
    GRAPH_ALGO_CORE_NUMBERS,
    GRAPH_ALGO_BELLMAN_FORD,
    GRAPH_ALGO_COUNT
} GraphAlgo;

const char *graph_algo_name(GraphAlgo algo);

typedef struct {
    GraphAlgo algo;
    int       source;        /* start/source vertex used (traversals, Dijkstra) */
    int       repetitions;
    double    list_ns;       /* best-of-N wall time on the list                 */
    double    matrix_ns;     /* best-of-N wall time on the matrix               */
    size_t    list_bytes;    /* heap footprint in each representation           */
    size_t    matrix_bytes;
    bool      list_ok;       /* the algorithm succeeded (e.g. a DAG for topo)   */
    bool      matrix_ok;
    bool      agreed;        /* both representations gave an equivalent answer  */
    GraphRep  faster;        /* the one with the lower time                     */
} GraphCompareResult;

/**
 * graph_compare - benchmark @algo on a private copy of @g in each
 * representation.  @g is not modified, its workload statistics are untouched
 * and it will not migrate because of this call.
 *
 * @source:      start vertex for BFS/DFS/Dijkstra; ignored otherwise.  If it is
 *               not alive the lowest alive ID is used.
 * @repetitions: runs per representation (<= 0 means 5); the best time is kept.
 *
 * "Agreed" means equivalent, not byte-identical: traversals may visit
 * neighbours in a different order, so BFS/DFS are compared by the set of
 * vertices reached, components by partition, Dijkstra by distances, MST by
 * total weight and edge count, topological sort by validity.
 *
 * Returns false (and leaves @out unspecified) on allocation failure or an
 * empty graph.  Complexity: O(V + E) per copy plus the timed runs.
 */
bool graph_compare(const Graph *g, GraphAlgo algo, int source, int repetitions,
                   GraphCompareResult *out);

/** graph_compare for every algorithm; @out needs GRAPH_ALGO_COUNT entries. */
bool graph_compare_all(const Graph *g, int source, int repetitions,
                       GraphCompareResult out[GRAPH_ALGO_COUNT]);

/**
 * graph_crossover_report - find where the matrix starts to beat the list.
 *
 * For every size in @sizes and density in @densities (ascending is easiest to
 * read) it builds a random graph (Erdos-Renyi; a random DAG for topological
 * sort), runs graph_compare for every GraphAlgo, and prints one line per
 * (algorithm, size, density) plus a summary of the density at which the faster
 * representation flips.  @csv selects machine-readable output without the
 * summary.  @repetitions as in graph_compare (<= 0 means 5).
 *
 * Takes seconds to minutes for large sizes: the matrix needs O(n^2) bits.
 * Returns false on bad arguments or allocation failure (see graph_last_error()).
 * Use -O2 for meaningful timings.
 */
bool graph_crossover_report(FILE *out, const size_t *sizes, size_t n_sizes,
                            const double *densities, size_t n_dens,
                            unsigned long long seed, int repetitions, bool csv);

/** Print one result as a short human-readable block. */
void graph_compare_print(const GraphCompareResult *result, FILE *out);

#endif /* GRAPH_H */
