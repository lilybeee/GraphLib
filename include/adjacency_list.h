/*
 * adjacency_list.h — Adjacency-list graph backend
 *
 * Vertices are integer IDs.  IDs live in the range [0, id_bound); an ID in
 * that range is either alive or dead (a removed vertex whose ID is waiting to
 * be reused by list_graph_add_vertex).  Every operation that takes a vertex
 * rejects dead and out-of-range IDs.
 *
 * Storage
 * -------
 *   - Vertex records live in fixed-size segments (LIST_SEG_SHIFT, default
 *     1024 vertices per segment).  Growing the graph allocates one new segment
 *     and never copies existing vertices; compaction frees empty segments.
 *   - Each vertex owns one contiguous block holding its neighbour IDs and
 *     edge weights, so scanning a neighbourhood is a linear walk over memory.
 *   - Every undirected edge is stored in both directions (a self-loop once),
 *     but edge_count counts it once.
 *
 * Removal is lazy: list_graph_remove_vertex frees the vertex's edges and puts
 * its ID on a free list.  Once enough IDs are dead,
 * list_graph_needs_compaction() turns true; list_graph_compact() then
 * renumbers the surviving vertices densely and reports old -> new IDs.
 * Compaction is explicit because it changes IDs; nothing else ever does.
 *
 * Array-size convention for the algorithms below:
 *   - arrays indexed by vertex ID (parent, component, distance) need
 *     list_graph_id_bound() entries; dead IDs receive -1 / INFINITY;
 *   - arrays listing vertices (order, edges) need list_graph_vertex_count()
 *     entries (vertex_count - 1 for the MST edges).
 * With no removals id_bound == vertex_count, as in a plain fixed-size graph.
 *
 * Compile with -lm (uses <math.h> for INFINITY).
 */

#ifndef ADJACENCY_LIST_H
#define ADJACENCY_LIST_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "graph_types.h"

/* ---- Opaque graph handle ---- */
typedef struct ListGraph ListGraph;

/* ---- Edge type (used by Kruskal's MST) ---- */
typedef GraphEdge ListEdge;

/* ------------------------------------------------------------------ */
/*  Construction / destruction                                          */
/* ------------------------------------------------------------------ */

/**
 * list_graph_create - allocate a new graph with vertices 0 .. vertices-1.
 * @vertices: initial number of vertices (may be 0; more can be added later)
 * @directed: true for a directed graph, false for undirected
 * Returns NULL on allocation failure or if vertices exceeds INT_MAX.
 * Complexity: O(V)
 */
ListGraph *list_graph_create(size_t vertices, bool directed);

/**
 * list_graph_destroy - free all memory owned by the graph.
 * Safe to call with NULL.
 * Complexity: O(V + E)
 */
void list_graph_destroy(ListGraph *graph);

/**
 * list_graph_create_unweighted - like list_graph_create, but the graph stores no
 * weights at all: every edge has weight 1.0, which is what list_graph_has_edge
 * and the iterators report, and adding an edge with any other weight fails with
 * GRAPH_ERR_INVALID_WEIGHT.  Each stored neighbour then costs 4 bytes instead of
 * 12 (about 25 bytes per friendship instead of 47), which is what a large
 * social network wants.  list_graph_neighbors returns NULL for the weights.
 */
ListGraph *list_graph_create_unweighted(size_t vertices, bool directed);

/**
 * list_graph_shrink_to_fit - release the spare room every adjacency block keeps
 * for growth (blocks double when full, so on average a third of the entry memory
 * is slack).  Contents are unchanged; pointers from list_graph_neighbors are
 * invalidated; adding edges afterwards works but regrows blocks.  Intended for
 * a graph that is finished being built.  Returns false (graph still valid) if
 * some block could not be reallocated.  Complexity: O(V + E).
 */
bool list_graph_shrink_to_fit(ListGraph *graph);

/** True for a graph made by list_graph_create_unweighted. */
bool list_graph_is_unweighted(const ListGraph *graph);

/* ------------------------------------------------------------------ */
/*  Queries                                                             */
/* ------------------------------------------------------------------ */

/** Returns the number of alive vertices (0 if graph is NULL). */
size_t list_graph_vertex_count(const ListGraph *graph);

/** Returns one past the largest ID ever handed out since the last compaction. */
size_t list_graph_id_bound(const ListGraph *graph);

/** Returns true iff v is an alive vertex. */
bool list_graph_is_alive(const ListGraph *graph, int v);

/**
 * Returns the number of edges.
 * Undirected edges are counted once.
 */
size_t list_graph_edge_count(const ListGraph *graph);

/** Returns true iff the graph is directed. */
bool list_graph_is_directed(const ListGraph *graph);

/** Returns the out-degree of v (0 if v is not alive). */
size_t list_graph_degree(const ListGraph *graph, int v);

/**
 * list_graph_neighbors - borrow v's adjacency arrays.
 * On success *to and *weight point at degree(v) contiguous entries (*weight is
 * NULL for an unweighted graph: every weight is 1.0) and the
 * degree is returned; otherwise both are set to NULL and 0 is returned.
 * The pointers are invalidated by the next modification of the graph.
 * Either output may be NULL if not needed.
 * Complexity: O(1)
 */
size_t list_graph_neighbors(const ListGraph *graph, int v,
                            const int **to, const double **weight);

/**
 * list_graph_next_neighbor - step through v's out-neighbours in ascending
 * ID order (rows are kept sorted).  Start with *cursor = 0.
 * On success fills *to (and *weight, if non-NULL), advances *cursor and returns
 * true; returns false when there are no more neighbours or v is not alive.
 * Do not modify the graph mid-walk.  Complexity: O(1) per neighbour.
 */
bool list_graph_next_neighbor(const ListGraph *graph, int v, size_t *cursor,
                              int *to, double *weight);

/**
 * Returns a counter that increases on every successful modification
 * (vertex/edge added or removed, compaction).  Equal values mean the graph is
 * unchanged, which is what a result cache needs to stay valid.
 */
uint64_t list_graph_version(const ListGraph *graph);

/** Approximate heap footprint of the graph in bytes. Complexity: O(V) */
size_t list_graph_memory_bytes(const ListGraph *graph);

/* ------------------------------------------------------------------ */
/*  Vertex modification                                                 */
/* ------------------------------------------------------------------ */

/**
 * list_graph_add_vertex - add an isolated vertex.
 * Reuses a dead ID if one exists, otherwise grows the ID space by one
 * (allocating a new segment when the current one is full).
 * Returns the new vertex ID, or -1 on allocation failure / ID exhaustion.
 * Complexity: amortized O(1)
 */
int list_graph_add_vertex(ListGraph *graph);

/**
 * list_graph_remove_vertex - remove v and every edge incident to it.
 * Returns false if v is not an alive vertex.  v's ID becomes reusable.
 * Complexity: O(degree(v) * degree of neighbours) for undirected graphs;
 *             O(V + E) for directed graphs (in-edges have no reverse index).
 */
bool list_graph_remove_vertex(ListGraph *graph, int v);

/**
 * list_graph_needs_compaction - true once more than 25% of the ID space is
 * dead.  Purely advisory; the caller decides when to compact.
 */
bool list_graph_needs_compaction(const ListGraph *graph);

/**
 * list_graph_compact - renumber alive vertices to 0 .. vertex_count-1,
 * preserving their relative order, and free the dead tail.
 * @remap: optional caller-provided array of size >= the *old* id_bound;
 *         remap[old] receives the new ID, or -1 if old was dead.
 * Returns false (graph unchanged) on allocation failure.
 * After a successful call id_bound == vertex_count and every previously held
 * vertex ID is invalid unless translated through remap.
 * Complexity: O(V + E)
 */
bool list_graph_compact(ListGraph *graph, int *remap);

/* ------------------------------------------------------------------ */
/*  Edge modification                                                   */
/* ------------------------------------------------------------------ */

/**
 * list_graph_add_edge - insert an edge (u -> v) with the given weight.
 * Returns false if:
 *   - u or v is not an alive vertex,
 *   - weight is not finite (INFINITY or NaN), or is not 1.0 in an unweighted graph,
 *   - the edge already exists.
 * For undirected graphs the reverse edge (v -> u) is also inserted.
 * Rows are kept sorted by neighbour ID.
 * Complexity: O(log degree(u)) to locate, O(degree(u)) worst case to shift (O(1) when
 * edges arrive in ascending order).
 */
bool list_graph_add_edge(ListGraph *graph, int u, int v, double weight);

/**
 * list_graph_add_new_edge - add_edge for bulk loading and migration.
 * The caller guarantees the edge (u -> v) does not exist, which lets this skip
 * the duplicate check.  Calling it for an existing edge stores a duplicate
 * and corrupts the edge count.  Same failure cases otherwise.
 * Complexity: as list_graph_add_edge
 */
bool list_graph_add_new_edge(ListGraph *graph, int u, int v, double weight);

/**
 * list_graph_remove_edge - delete the edge (u -> v).
 * Returns false if the edge does not exist or a vertex is not alive.
 * Rows stay sorted.
 * Complexity: O(log degree(u)) to locate, O(degree(u)) to shift
 */
bool list_graph_remove_edge(ListGraph *graph, int u, int v);

/**
 * list_graph_has_edge - test whether (u -> v) exists.
 * If weight_out is non-NULL and the edge exists, *weight_out is filled.
 * Complexity: O(log degree(u)) (binary search; rows are sorted)
 */
bool list_graph_has_edge(const ListGraph *graph, int u, int v, double *weight_out);

/**
 * list_graph_print - print the adjacency list to stdout for debugging.
 * Complexity: O(V + E)
 */
void list_graph_print(const ListGraph *graph);

/* ------------------------------------------------------------------ */
/*  Traversals                                                          */
/* ------------------------------------------------------------------ */

/**
 * list_bfs - breadth-first search from 'start'.
 * @order:  caller-provided array of size >= vertex_count; receives vertices
 *          in visit order.
 * @parent: caller-provided array of size >= id_bound, or NULL.
 *          parent[start] = -1; parent[v] = the vertex that discovered v;
 *          unreached and dead IDs are -1.
 * Returns the number of vertices visited.
 * Complexity: O(V + E)
 */
size_t list_bfs(const ListGraph *graph, int start, int *order, int *parent);

/**
 * list_dfs - depth-first search from 'start' (iterative: no recursion-depth
 * limit and no per-vertex allocation).  Neighbours are explored in storage
 * order, so the result matches a recursive pre-order DFS.
 * @order:  same convention as list_bfs.
 * @parent: same convention as list_bfs.
 * Returns the number of vertices visited.
 * Complexity: O(V + E)
 */
size_t list_dfs(const ListGraph *graph, int start, int *order, int *parent);

/* ------------------------------------------------------------------ */
/*  Analysis                                                            */
/* ------------------------------------------------------------------ */

/**
 * list_connected_components - label every alive vertex with its component.
 * Works on both directed (weakly) and undirected graphs.
 * @component: caller-provided array of size >= id_bound.
 *             component[v] = 0-based component index, numbered in order of
 *             each component's smallest vertex ID; dead IDs get -1.
 * Returns the total number of components (0 if graph is NULL).
 * Complexity: O(V + E) for undirected, O((V + E) alpha(V)) for directed.
 */
size_t list_connected_components(const ListGraph *graph, int *component);

/**
 * list_has_cycle - test whether the graph contains a cycle (iterative DFS).
 * For directed graphs: back-edge detection (true cycle).
 * For undirected graphs: any non-tree edge, including a self-loop.
 * Returns false if graph is NULL or empty.
 * Complexity: O(V + E)
 */
bool list_has_cycle(const ListGraph *graph);

/**
 * list_topological_sort - topological ordering of a DAG (Kahn's algorithm).
 * @order: caller-provided array of size >= vertex_count.
 *         Receives vertices in topological order.
 * Returns true on success.
 * Returns false if:
 *   - the graph is undirected,
 *   - a cycle exists,
 *   - graph is NULL,
 *   - vertex_count > 0 but order is NULL.
 * Complexity: O(V + E)
 */
bool list_topological_sort(const ListGraph *graph, int *order);

/**
 * list_dijkstra - single-source shortest paths (non-negative weights).
 * @source:   source vertex.
 * @distance: caller-provided array of size >= id_bound.
 *            distance[v] = shortest distance from source to v, or INFINITY.
 * @parent:   caller-provided array of size >= id_bound, or NULL.
 *            parent[source] = -1; parent[v] = predecessor of v on shortest path.
 * Returns false if source is not alive, distance is NULL, or any edge weight
 * in the graph is negative (the graph tracks this, so the check is O(1)).
 * Complexity: O((V + E) log V) with an indexed binary heap, O(V) extra space.
 */
bool list_dijkstra(const ListGraph *graph, int source, double *distance, int *parent);

/**
 * list_mst_kruskal - minimum spanning tree (undirected graphs only).
 * @edges:      caller-provided array of size >= vertex_count - 1.
 * @edge_count: receives the number of MST edges written.
 * Returns the total MST weight, or INFINITY if:
 *   - graph is directed,
 *   - graph is disconnected (so no spanning tree exists),
 *   - graph or edges is NULL (when vertex_count > 1).
 * Complexity: O(E log E) for the sort + near-O(E) union-find.
 */
double list_mst_kruskal(const ListGraph *graph, ListEdge *edges, size_t *edge_count);

#endif /* ADJACENCY_LIST_H */
