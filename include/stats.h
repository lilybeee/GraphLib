/*
 * stats.h — workload counters and graph-shape summaries.
 *
 * The facade bumps a counter for every operation it performs; the observer
 * (policy.h) turns those counters plus the current graph shape into a
 * predicted cost per representation.  Counters are doubles so they can be
 * decayed: old behaviour fades out instead of being cut off at a window edge.
 */

#ifndef STATS_H
#define STATS_H

#include <stdbool.h>
#include <stddef.h>

typedef enum {
    OP_HAS_EDGE,
    OP_ADD_EDGE,
    OP_REMOVE_EDGE,
    OP_ADD_VERTEX,
    OP_REMOVE_VERTEX,
    OP_NEIGHBOR_SCAN,   /* walking one vertex's neighbours (friends-of-v, ...) */
    OP_ALGO_RUN,        /* a traversal-style algorithm that executed: BFS, DFS,
                           components, cycle check, topological sort (no weights) */
    OP_WEIGHTED_RUN,    /* an algorithm that reads every edge weight: Dijkstra, MST */
    OP_ITERATIVE_RUN,   /* repeated full sweeps over every edge: PageRank, HITS, Louvain,
                           k-core, SCC, Bellman-Ford, label propagation          */
    OP_ALL_PAIRS_RUN,   /* one traversal per vertex: closeness, betweenness,
                           eccentricity, diameter, transitive closure            */
    OP_ADJACENCY_RUN,   /* many "is (a, b) an edge?" probes between a vertex's
                           neighbours: triangles, cliques, flow, isomorphism     */
    OP_CACHE_HIT,       /* a result served from the cache (costs the graph 0)  */
    OP_CACHE_MISS,
    OP_COUNT
} OpKind;

typedef struct {
    double n[OP_COUNT];
} Workload;

/* What the graph looks like right now; the inputs every cost model needs. */
typedef struct {
    size_t live;        /* alive vertices                                   */
    size_t id_bound;    /* size of the ID space (matrix side length)        */
    size_t edges;       /* undirected edges counted once                    */
    bool   directed;
    bool   unweighted;  /* stores no weights (every edge is 1.0)            */
} GraphShape;

/** Number of stored adjacency entries (an undirected edge is stored twice). */
static inline double shape_entries(const GraphShape *s)
{
    return s->directed ? (double)s->edges : 2.0 * (double)s->edges;
}

/** Average out-degree over alive vertices. */
static inline double shape_avg_degree(const GraphShape *s)
{
    return s->live == 0 ? 0.0 : shape_entries(s) / (double)s->live;
}

/** Fraction of possible directed arcs that exist, in [0, 1]. */
static inline double shape_density(const GraphShape *s)
{
    const double n = (double)s->live;
    return s->live < 2 ? 0.0 : shape_entries(s) / (n * (n - 1.0));
}

void   workload_note(Workload *w, OpKind op);
void   workload_decay(Workload *w, double factor);
void   workload_reset(Workload *w);

/** Sum of every counter except cache hits/misses (those cost the graph nothing). */
double workload_graph_ops(const Workload *w);

/** Cache hit rate in [0, 1]; 0 when nothing was looked up. */
double workload_hit_rate(const Workload *w);

#endif /* STATS_H */
