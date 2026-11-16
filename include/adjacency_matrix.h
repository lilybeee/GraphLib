/*
 * adjacency_matrix.h — Adjacency-matrix graph backend
 *
 * Vertices are integer IDs.  IDs live in the range [0, id_bound); an ID in
 * that range is either alive or dead (a removed vertex whose ID is waiting to
 * be reused by matrix_graph_add_vertex).  Every operation that takes a vertex
 * rejects dead and out-of-range IDs.  The API mirrors adjacency_list.h, so the
 * two backends can be driven by the same code and the same tests.
 *
 * Storage
 * -------
 *   The matrix is stored as two cooperating structures:
 *     - presence bits, in strips of 64 rows.  A row's bits are contiguous
 *       (V / 64 words), so a traversal tests 64 neighbours per word operation
 *       and streams through memory.  A strip is allocated the first time an
 *       edge leaves one of its 64 rows.
 *     - weights, in 64 x 64 tiles allocated the first time an edge lands in
 *       that block.  Tiles never move, so growing the graph copies only the
 *       small bit strips (1 bit per cell), never the 8-byte weights.
 *   Memory therefore tracks the edges: the weights of an empty region cost one
 *   pointer, and the bits cost V^2 / 8 bytes at most.  Undirected edges are
 *   stored in both directions (a self-loop once), but edge_count counts them
 *   once.
 *
 * Removal is lazy, exactly as in the list backend: matrix_graph_remove_vertex
 * clears the vertex's row and column and puts its ID on a free list;
 * matrix_graph_needs_compaction() turns true once more than 25% of the IDs are
 * dead; matrix_graph_compact() then renumbers the survivors densely and reports
 * old -> new IDs.  Compaction is the only operation that changes IDs.
 *
 * Array-size convention for the algorithms below:
 *   - arrays indexed by vertex ID (parent, component, distance) need
 *     matrix_graph_id_bound() entries; dead IDs receive -1 / INFINITY;
 *   - arrays listing vertices (order, edges) need matrix_graph_vertex_count()
 *     entries (vertex_count - 1 for the MST edges).
 * With no removals id_bound == vertex_count, as in a plain fixed-size graph.
 *
 * Neighbours are always visited in ascending ID order.
 *
 * Compile with -lm (uses <math.h> for INFINITY).
 */

#ifndef ADJACENCY_MATRIX_H
#define ADJACENCY_MATRIX_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "graph_types.h"

typedef struct MatrixGraph MatrixGraph;

typedef GraphEdge MatrixEdge;

/* ------------------------------------------------------------------ */
/*  Construction / destruction                                          */
/* ------------------------------------------------------------------ */

/**
 * matrix_graph_create - allocate a new graph with vertices 0 .. vertices-1.
 * @vertices: initial number of vertices (may be 0; more can be added later)
 * @directed: true for a directed graph, false for undirected
 * Returns NULL on allocation failure, or if the tile directory for that many
 * vertices cannot be represented.
 * Complexity: O((V / 64)^2) for the tile directory
 */
MatrixGraph *matrix_graph_create(size_t vertices, bool directed);

/** Frees all memory owned by the graph.  Safe to call with NULL. */
void matrix_graph_destroy(MatrixGraph *graph);

/**
 * matrix_graph_create_unweighted - like matrix_graph_create, but the graph keeps only the
 * presence bits: no 32 KB weight tiles, so every edge reads back as weight 1.0 and adding an
 * edge with any other weight fails with GRAPH_ERR_INVALID_WEIGHT.  An edge then costs one bit
 * instead of (potentially) a whole tile.
 */
MatrixGraph *matrix_graph_create_unweighted(size_t vertices, bool directed);

/** True for a graph made by matrix_graph_create_unweighted. */
bool matrix_graph_is_unweighted(const MatrixGraph *graph);

/**
 * matrix_max_vertices - the most vertices a matrix may represent under @budget bytes
 * (0 = unlimited, capped at INT_MAX).  Once every row strip exists the presence bits take
 * V^2/8 bytes, so this is the largest V whose full matrix fits: sqrt(8 * share * budget),
 * with share 0.5 for a weighted graph (half the budget is kept back for weight tiles) and 0.9
 * for an unweighted one.  Roughly 131 000 vertices weighted and 176 000 unweighted under the
 * default 4 GiB.  A matrix larger than this is refused (GRAPH_ERR_TOO_LARGE) however few edges
 * it starts with, because it could not be filled within the budget.
 */
size_t matrix_max_vertices(size_t budget, bool weighted);

/* ------------------------------------------------------------------ */
/*  Queries                                                             */
/* ------------------------------------------------------------------ */

/** Returns the number of alive vertices (0 if graph is NULL). */
size_t matrix_graph_vertex_count(const MatrixGraph *graph);

/** Returns one past the largest ID ever handed out since the last compaction. */
size_t matrix_graph_id_bound(const MatrixGraph *graph);

/** Returns true iff v is an alive vertex. */
bool matrix_graph_is_alive(const MatrixGraph *graph, int v);

/** Returns the number of edges; undirected edges are counted once. */
size_t matrix_graph_edge_count(const MatrixGraph *graph);

/** Returns true iff the graph is directed. */
bool matrix_graph_is_directed(const MatrixGraph *graph);

/** Returns the out-degree of v (0 if v is not alive).  Complexity: O(1) */
size_t matrix_graph_degree(const MatrixGraph *graph, int v);

/**
 * matrix_graph_next_neighbor - step through v's out-neighbours in ascending ID
 * order.  Start with *cursor = 0.  On success fills *to (and *weight, if
 * non-NULL), advances *cursor past it, and returns true; returns false when
 * there are no more neighbours or v is not alive.
 * Complexity: O(1) per neighbour plus O(V / 64) over the whole walk.
 */
bool matrix_graph_next_neighbor(const MatrixGraph *graph, int v, size_t *cursor,
                                int *to, double *weight);

/**
 * Returns a counter that increases on every successful modification.  Equal
 * values mean the graph is unchanged (what a result cache needs).
 */
uint64_t matrix_graph_version(const MatrixGraph *graph);

/** Approximate heap footprint of the graph in bytes.  Complexity: O(1) */
size_t matrix_graph_memory_bytes(const MatrixGraph *graph);

/* ------------------------------------------------------------------ */
/*  Vertex modification                                                 */
/* ------------------------------------------------------------------ */

/**
 * matrix_graph_add_vertex - add an isolated vertex.
 * Reuses a dead ID if one exists, otherwise grows the ID space by one,
 * enlarging the tile directory when the current capacity is full.
 * Returns the new vertex ID, or -1 on allocation failure / ID exhaustion.
 * Complexity: amortized O(1) (directory growth is O((V / 64)^2), doubling)
 */
int matrix_graph_add_vertex(MatrixGraph *graph);

/**
 * matrix_graph_remove_vertex - remove v and every edge incident to it.
 * Returns false if v is not an alive vertex.  v's ID becomes reusable.
 * Complexity: O(V)
 */
bool matrix_graph_remove_vertex(MatrixGraph *graph, int v);

/**
 * matrix_graph_needs_compaction - true once more than 25% of the ID space is
 * dead.  Purely advisory; the caller decides when to compact.
 */
bool matrix_graph_needs_compaction(const MatrixGraph *graph);

/**
 * matrix_graph_compact - renumber alive vertices to 0 .. vertex_count-1,
 * preserving their relative order, and release the dead rows, columns and any
 * tiles left empty.
 * @remap: optional caller-provided array of size >= the *old* id_bound;
 *         remap[old] receives the new ID, or -1 if old was dead.
 * Returns false (graph unchanged) on allocation failure.  Builds the compacted
 * matrix beside the old one, so it transiently needs both.
 * Complexity: O(V + E + tiles)
 */
bool matrix_graph_compact(MatrixGraph *graph, int *remap);

/* ------------------------------------------------------------------ */
/*  Edge modification                                                   */
/* ------------------------------------------------------------------ */

/**
 * matrix_graph_add_edge - insert an edge (source -> destination).
 * Returns false if:
 *   - an endpoint is not an alive vertex,
 *   - weight is not finite (INFINITY or NaN), or is not 1.0 in an unweighted graph,
 *   - the edge already exists,
 *   - a tile that the edge needs cannot be allocated.
 * For undirected graphs the reverse edge is also inserted.
 * Complexity: O(1), plus one tile allocation the first time a block is used.
 */
bool matrix_graph_add_edge(MatrixGraph *graph, int source, int destination, double weight);

/**
 * matrix_graph_remove_edge - delete the edge (source -> destination).
 * Returns false if the edge does not exist or an endpoint is not alive.
 * Complexity: O(1)
 */
bool matrix_graph_remove_edge(MatrixGraph *graph, int source, int destination);

/**
 * matrix_graph_has_edge - test whether (source -> destination) exists.
 * If weight_out is non-NULL and the edge exists, *weight_out is filled.
 * Complexity: O(1)
 */
bool matrix_graph_has_edge(const MatrixGraph *graph, int source, int destination,
                           double *weight_out);

/**
 * matrix_graph_print - print the matrix of alive vertices to stdout, with "."
 * for a missing edge.  Intended for small graphs.  Complexity: O(V^2)
 */
void matrix_graph_print(const MatrixGraph *graph);

/* ------------------------------------------------------------------ */
/*  Traversals                                                          */
/* ------------------------------------------------------------------ */

/**
 * matrix_bfs - breadth-first search from 'start' (bitset frontier expansion).
 * @order:  caller-provided array of size >= vertex_count; receives vertices
 *          in visit order.
 * @parent: caller-provided array of size >= id_bound, or NULL.
 *          parent[start] = -1; parent[v] = the vertex that discovered v;
 *          unreached and dead IDs are -1.
 * Returns the number of vertices visited.
 * Complexity: O(V * V / 64 + E)
 */
size_t matrix_bfs(const MatrixGraph *graph, int start, int *order, int *parent);

/**
 * matrix_dfs - depth-first search from 'start' (iterative, true pre-order: a
 * vertex is marked when it is first visited, so parent[] is the DFS tree).
 * @order:  same convention as matrix_bfs.
 * @parent: same convention as matrix_bfs.
 * Returns the number of vertices visited.
 * Complexity: O(V * V / 64 + E)
 */
size_t matrix_dfs(const MatrixGraph *graph, int start, int *order, int *parent);

/* ------------------------------------------------------------------ */
/*  Analysis                                                            */
/* ------------------------------------------------------------------ */

/**
 * matrix_connected_components - label every alive vertex with its component.
 * Works on both directed (weakly) and undirected graphs.
 * @component: caller-provided array of size >= id_bound.
 *             component[v] = 0-based component index, numbered in order of
 *             each component's smallest vertex ID; dead IDs get -1.
 * Returns the total number of components (0 if graph is NULL).
 * Complexity: O(V * V / 64 + E) undirected, O(tiles * 64 + E) directed.
 */
size_t matrix_connected_components(const MatrixGraph *graph, int *component);

/**
 * matrix_has_cycle - test whether the graph contains a cycle (iterative DFS).
 * For directed graphs: back-edge detection (true cycle).
 * For undirected graphs: any non-tree edge, including a self-loop.
 * Returns false if graph is NULL or empty.
 * Complexity: O(V * V / 64 + E)
 */
bool matrix_has_cycle(const MatrixGraph *graph);

/**
 * matrix_topological_sort - topological ordering of a DAG (Kahn's algorithm).
 * @order: caller-provided array of size >= vertex_count.
 * Returns true on success.
 * Returns false if:
 *   - the graph is undirected,
 *   - a cycle exists,
 *   - graph is NULL,
 *   - vertex_count > 0 but order is NULL.
 * Complexity: O(tiles * 64 + E)
 */
bool matrix_topological_sort(const MatrixGraph *graph, int *order);

/**
 * matrix_dijkstra - single-source shortest paths (non-negative weights).
 * @source:   source vertex.
 * @distance: caller-provided array of size >= id_bound.
 *            distance[v] = shortest distance from source to v, or INFINITY.
 * @parent:   caller-provided array of size >= id_bound, or NULL.
 *            parent[source] = -1; parent[v] = predecessor of v on shortest path.
 * Returns false if source is not alive, distance is NULL, or any edge weight
 * in the graph is negative (the graph tracks this, so the check is O(1)).
 * Complexity: O(V * V / 64 + (V + E) log V) with an indexed binary heap.
 */
bool matrix_dijkstra(const MatrixGraph *graph, int source, double *distance, int *parent);

/**
 * matrix_mst_kruskal - minimum spanning tree (undirected graphs only).
 * @edges:      caller-provided array of size >= vertex_count - 1.
 * @edge_count: receives the number of MST edges written.
 * Returns the total MST weight, or INFINITY if:
 *   - graph is directed,
 *   - graph is disconnected (so no spanning tree exists),
 *   - graph or edges is NULL (when vertex_count > 1).
 * Complexity: O(tiles * 64 + E log E)
 */
double matrix_mst_kruskal(const MatrixGraph *graph, MatrixEdge *edges, size_t *edge_count);

#endif /* ADJACENCY_MATRIX_H */
