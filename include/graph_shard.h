/*
 * graph_shard.h — the sharding engine: run on a graph bigger than the memory you allow.
 *
 * WHEN TO USE IT.  Normally the library keeps the whole graph in RAM (about 25-50 bytes per
 * friendship).  If the graph does not fit the memory budget -- a small device such as a
 * Raspberry Pi, or just a graph with hundreds of millions of edges -- turn the engine on:
 *
 *     graph_set_memory_budget(64u << 20);              // the limit you want, here 64 MB
 *     ShardConfig cfg = shard_default_config();        // budget defaults to the one above
 *     cfg.force = true;                                // optional: shard even if it would fit
 *     Graph *g; ShardedGraph *sg;
 *     graph_open_edge_list(file, &cfg, &g, &sg);       // exactly one of g / sg comes back
 *     if (sg) { shard_pagerank(sg, ...); shard_destroy(sg); }
 *
 * HOW IT WORKS (and why).  This is the "semi-external" model:
 *   - The EDGES -- by far the biggest part -- live in files on disk, cut into SHARDS: each shard
 *     holds the adjacency lists of a contiguous range of user IDs.  They are written once, when
 *     the graph is opened, and never change (the graph is fixed while it is analysed).
 *   - The per-user STATE of an algorithm (PageRank values, component labels, BFS levels) stays
 *     in RAM, 4-24 bytes per user.  This is the key decision: because every user's state is
 *     always in memory, a shard can name any neighbour, in any other shard, without that
 *     shard being loaded -- so there are no "ghost vertices" or cross-shard messages.
 *   - Because shards are read-only, evicting one is just freeing its memory: nothing is ever
 *     written back, and the only "lock" needed is a pin that keeps a shard from being evicted
 *     while it is in use.  At most `cache_capacity` shards are resident at once (an LRU cache
 *     sized from the budget).
 *   - An algorithm is a SWEEP: for each shard in turn, pin it, process its rows, unpin it.
 *     Each sweep reads every shard once, so cost is (sweeps x graph size) of sequential disk
 *     reads.  That is why this mode is slower than in-memory processing -- on purpose: it is
 *     the price of fitting in a small budget.
 *
 * WHAT IT CAN AND CANNOT DO.  It supports whole-graph questions that are sweeps over the edges
 * (degrees, components, PageRank, hop-limited BFS).  Algorithms that jump unpredictably
 * between users (Dijkstra, A*, max flow, betweenness) would reload shards constantly and are
 * not provided; for those the graph has to fit in memory.  The engine ignores edge weights
 * (every edge counts as 1) and merges repeated edges exactly as graph_add_edge does for an
 * in-memory graph, so both modes see the same graph; an undirected file's edges are stored as two
 * arcs.
 *
 * MEMORY.  The budget covers the shard cache plus the per-user arrays; opening fails with
 * GRAPH_ERR_TOO_LARGE if even that cannot fit (about 24 bytes per user are needed, so a 64 MB
 * budget allows roughly 1.4 million users).
 */

#ifndef GRAPH_SHARD_H
#define GRAPH_SHARD_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include "graph.h"
#include "graph_status.h"

typedef struct ShardedGraph ShardedGraph;

typedef struct {
    size_t      memory_budget;  /* bytes the engine may use in RAM; 0 = graph_memory_budget() */
    const char *dir;            /* existing directory for the shard files; NULL = a private
                                   temporary directory that shard_destroy() removes           */
    bool        force;          /* graph_open_edge_list: shard even if the graph would fit    */
} ShardConfig;

/** budget 0 (the library budget), temporary directory, no forcing. */
ShardConfig shard_default_config(void);

/**
 * Cut the edge list in @in into shard files under the budget and return a handle.  Makes three
 * passes over the file (so it must be seekable) and holds only O(users) integers and a few write
 * buffers in memory while doing it.  Errors: GRAPH_ERR_IO (malformed or unwritable),
 * GRAPH_ERR_TOO_LARGE (the users alone exceed the budget), GRAPH_ERR_OUT_OF_MEMORY.
 */
GraphStatus shard_build(FILE *in, const ShardConfig *config, ShardedGraph **out);

/**
 * The mode switch.  Reads the size of the edge list; if the graph fits the budget as an
 * in-memory list and config->force is false, loads it into an ordinary Graph (*graph set,
 * *sharded NULL); otherwise shards it (*sharded set, *graph NULL).  The caller owns whichever
 * came back (graph_destroy / shard_destroy).
 */
GraphStatus graph_open_edge_list(FILE *in, const ShardConfig *config, Graph **graph,
                                 ShardedGraph **sharded);

/** Deletes the shard files (and the directory, if the engine made it) and frees the handle. */
void shard_destroy(ShardedGraph *sg);

/* ---- what the handle knows ---- */

size_t shard_vertex_count(const ShardedGraph *sg);   /* users: largest ID + 1 */
size_t shard_edge_count(const ShardedGraph *sg);     /* distinct edges (repeats merged) */
bool   shard_is_directed(const ShardedGraph *sg);
size_t shard_degree(const ShardedGraph *sg, int user);   /* out-degree (both ends if undirected) */

typedef struct {
    size_t   shards;          /* how many shard files                              */
    size_t   cache_capacity;  /* shards that may be resident at once               */
    size_t   shard_bytes_max; /* size of the biggest shard in memory               */
    size_t   loads;           /* shard reads from disk so far                      */
    uint64_t bytes_read;      /* bytes read from disk so far                       */
} ShardStats;

ShardStats shard_stats(const ShardedGraph *sg);

/* ---- algorithms: each is a sweep over the shards ---- */

/** out_degree[v] for every user (n entries). */
GraphStatus shard_degrees(const ShardedGraph *sg, size_t *out_degree);

/** Connected components, directions ignored; labels 0.. in order of first member (n entries). */
GraphStatus shard_components(ShardedGraph *sg, int *component, size_t *count);

/**
 * PageRank by power iteration, the same definition as graph_pagerank on an unweighted graph
 * (dangling users spread their rank evenly).  rank has n entries and sums to 1.
 * GRAPH_ERR_NO_CONVERGENCE if the L1 change stays above @tolerance after @max_iterations
 * (rank still holds the latest estimate).  *iterations (optional) = rounds used.
 */
GraphStatus shard_pagerank(ShardedGraph *sg, double damping, double tolerance, size_t max_iterations,
                           double *rank, size_t *iterations);

/**
 * Hop distances from @source, giving up beyond @max_hops (SIZE_MAX = no limit): level[v] is the
 * distance, or -1 if not reached (n entries).  *reached (optional) counts users reached,
 * including the source.  Only shards that contain users of the current level are read.
 */
GraphStatus shard_bfs(ShardedGraph *sg, int source, size_t max_hops, int *level, size_t *reached);

#endif /* GRAPH_SHARD_H */
