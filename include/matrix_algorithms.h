/*
 * matrix_algorithms.h — algorithms on the adjacency-matrix backend.
 *
 * Everything here takes a `const MatrixGraph *` and never modifies it, so any
 * number of threads may run these functions on one graph at the same time as
 * long as nobody mutates it.
 *
 * Conventions
 * -----------
 *   - Vertices are IDs in [0, matrix_graph_id_bound()); dead IDs are skipped.
 *   - Arrays indexed by vertex ID need matrix_graph_id_bound() entries.  Dead
 *     IDs receive -1 (INFINITY for distances, 0 for scores) unless noted.
 *   - "Undirected only" functions return GRAPH_ERR_UNSUPPORTED (or the
 *     documented sentinel) for directed graphs.
 *   - A function that allocates its result says so and names the matching
 *     free function.  Nothing else allocates on the caller's behalf.
 *   - Complexities use V = alive vertices, E = edges, N = id_bound.  Walking
 *     one row of the matrix costs O(N / 64 + degree).
 *
 * Link with -lm.
 */

#ifndef MATRIX_ALGORITHMS_H
#define MATRIX_ALGORITHMS_H

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>

#include "adjacency_matrix.h"
#include "graph_results.h"
#include "graph_status.h"

/* ------------------------------------------------------------------ */
/*  Shared result types                                                 */
/* ------------------------------------------------------------------ */

/** The status type of the algorithm library; see graph_status.h for the codes. */
typedef GraphStatus MatrixStatus;


/** A list of vertex sets (components, cliques, ...).  Free with matrix_groups_free. */
typedef GraphGroups MatrixGroups;

void matrix_groups_free(MatrixGroups *groups);

/** A path with its total weight.  Free with matrix_path_free. */
typedef GraphPath MatrixPath;

void matrix_path_free(MatrixPath *path);

/* ================================================================== */
/*  1. Traversal with callbacks                                         */
/* ================================================================== */

/** Called once per vertex as it is visited; return false to stop early. */
typedef GraphVisitFn MatrixVisitFn;

/**
 * matrix_bfs_visit / matrix_dfs_visit - traverse from @start, calling @fn per
 * vertex.  DFS is iterative pre-order.  No output arrays are needed, so memory
 * is O(N) bits plus the frontier.  Returns the number of vertices visited
 * (0 if @start is not alive or @fn is NULL).  Complexity: O(V + E + N * V / 64).
 */
size_t matrix_bfs_visit(const MatrixGraph *g, int start, MatrixVisitFn fn, void *ctx);
size_t matrix_dfs_visit(const MatrixGraph *g, int start, MatrixVisitFn fn, void *ctx);

/* ================================================================== */
/*  2. Structure and analysis                                           */
/* ================================================================== */

/** Out-degree / (V - 1).  out: N entries. */
MatrixStatus matrix_degree_centrality(const MatrixGraph *g, double *out);

/**
 * Closeness by hop count, Wasserman-Faust normalised so disconnected graphs
 * work: (r / (V - 1)) * (r / sum_of_distances), r = vertices reachable from v.
 * 0 for an isolated vertex.  Complexity: O(V * (V + E)).
 */
MatrixStatus matrix_closeness_centrality(const MatrixGraph *g, double *out);

/**
 * Betweenness (Brandes), unweighted shortest paths.  Raw pair counts; for an
 * undirected graph each pair is counted once.  Complexity: O(V * (V + E)).
 */
MatrixStatus matrix_betweenness_centrality(const MatrixGraph *g, double *out);

/** Number of triangles.  Undirected only.  Complexity: O(sum of degree^2). */
MatrixStatus matrix_triangle_count(const MatrixGraph *g, long long *total);

/**
 * Local clustering coefficient per vertex (0 when degree < 2) and/or the
 * average over alive vertices (zeros included).  Either output may be NULL.
 * Undirected only.
 */
MatrixStatus matrix_clustering_coefficient(const MatrixGraph *g, double *local, double *average);

/**
 * 2-colourability of the underlying undirected graph.  @color (optional, N
 * entries) receives 0/1, or -1 for dead IDs; it is only meaningful if *result.
 */
MatrixStatus matrix_is_bipartite(const MatrixGraph *g, bool *result, int *color);

/** Largest hop distance from each vertex to a vertex it can reach.  N entries. */
MatrixStatus matrix_eccentricity(const MatrixGraph *g, int *ecc);

/**
 * Diameter and radius in hops.  GRAPH_ERR_UNSUPPORTED if some vertex cannot
 * reach another (not connected / not strongly connected); INVALID if empty.
 */
MatrixStatus matrix_diameter_radius(const MatrixGraph *g, int *diameter, int *radius);

/** All edges reversed; same IDs.  *out is a new graph (matrix_graph_destroy). */
MatrixStatus matrix_transpose(const MatrixGraph *g, MatrixGraph **out);

/** Common out-neighbours of u and v, ascending.  out: N entries. */
MatrixStatus matrix_common_neighbors(const MatrixGraph *g, int u, int v, int *out, size_t *count);

/** |N(u) ∩ N(v)| / |N(u) ∪ N(v)| over out-neighbours; 0 when both are empty. */
MatrixStatus matrix_jaccard(const MatrixGraph *g, int u, int v, double *score);

typedef GraphNetworkStats MatrixNetworkStats;

MatrixStatus matrix_network_stats(const MatrixGraph *g, MatrixNetworkStats *out);

/** hist[d] = number of vertices with out-degree d; degrees >= len-1 land in hist[len-1]. */
MatrixStatus matrix_degree_histogram(const MatrixGraph *g, size_t *hist, size_t len);

/* ---- connectivity ---- */

bool matrix_is_connected(const MatrixGraph *g);            /* weakly, if directed; empty = true */
bool matrix_is_strongly_connected(const MatrixGraph *g);
bool matrix_is_dag(const MatrixGraph *g);                  /* directed and acyclic              */
bool matrix_is_tree(const MatrixGraph *g);                 /* undirected, connected, E = V - 1  */
bool matrix_is_forest(const MatrixGraph *g);               /* undirected and acyclic            */

/**
 * Strongly connected components (Tarjan, iterative).  component: N entries,
 * numbered in reverse topological order of the condensation (component 0 has
 * no outgoing edges to other components).  Undirected graphs get their
 * ordinary components.  *count (optional) receives the number of components.
 */
MatrixStatus matrix_strongly_connected_components(const MatrixGraph *g, int *component, size_t *count);

/** Undirected only.  is_cut: N entries; *count (optional) = number of cut vertices. */
MatrixStatus matrix_articulation_points(const MatrixGraph *g, bool *is_cut, size_t *count);

/** Undirected only.  out needs matrix_graph_edge_count() entries; *count = bridges found. */
MatrixStatus matrix_bridges(const MatrixGraph *g, MatrixEdge *out, size_t *count);

/**
 * Biconnected components, undirected only.  Each group is the vertex set of
 * one block (an articulation point appears in every block it joins).
 * Isolated vertices belong to no block.  Free with matrix_groups_free.
 */
MatrixStatus matrix_biconnected_components(const MatrixGraph *g, MatrixGroups *out);

/**
 * Find a cycle.  Writes its vertices in order into @cycle (the last connects
 * back to the first) and returns its length; 0 if the graph is acyclic.  A
 * self-loop is a cycle of length 1.  @cycle needs matrix_graph_vertex_count()
 * entries.  Complexity: O(V + E).
 */
size_t matrix_find_cycle(const MatrixGraph *g, int *cycle);

/**
 * Edge (u, v), u != v, with weight 1.0 iff v is reachable from u.  New graph,
 * same IDs and directedness.  Complexity: O(V * (V + E)).
 */
MatrixStatus matrix_transitive_closure(const MatrixGraph *g, MatrixGraph **out);

/** Smallest graph with the same reachability.  DAGs only.  Edge weights are kept. */
MatrixStatus matrix_transitive_reduction(const MatrixGraph *g, MatrixGraph **out);

/* ================================================================== */
/*  3. Shortest paths                                                   */
/* ================================================================== */

/**
 * Rebuild the vertex list src..dst from a parent array (as filled by the BFS,
 * Dijkstra and Bellman-Ford functions).  @path needs @cap entries (N is
 * always enough).  Returns the length, or 0 if dst is unreachable or cap is
 * too small.
 */
size_t matrix_path_from_parent(const int *parent, size_t n, int src, int dst, int *path, size_t cap);

/** Same for the `next` array of matrix_floyd_warshall / matrix_johnson. */
size_t matrix_path_from_next(const int *next, size_t n, int u, int v, int *path, size_t cap);

/**
 * Bellman-Ford; handles negative weights.  distance/parent: N entries
 * (INFINITY / -1 when unreachable).  GRAPH_ERR_NEGATIVE_CYCLE if one is
 * reachable from @src (any negative undirected edge is a negative cycle).
 * Complexity: O(V * (V + E)).
 */
MatrixStatus matrix_bellman_ford(const MatrixGraph *g, int src, double *distance, int *parent);

/**
 * All-pairs shortest paths.  distance and next (optional) are N*N row-major
 * arrays: distance[u*N+v], next[u*N+v] = first hop after u on a shortest path
 * to v (-1 if none).  distance[u*N+u] = 0.  GRAPH_ERR_TOO_LARGE if N > 8192.
 * Complexity: O(V^3).
 */
MatrixStatus matrix_floyd_warshall(const MatrixGraph *g, double *distance, int *next);

/**
 * All-pairs shortest paths for sparse graphs with negative weights (Johnson):
 * Bellman-Ford once, then Dijkstra from every vertex.  Same output layout as
 * matrix_floyd_warshall.  Complexity: O(V * (V + E) log V).
 */
MatrixStatus matrix_johnson(const MatrixGraph *g, double *distance, int *next);

/** Estimate of the remaining cost from @vertex to @target; must never overestimate. */
typedef GraphHeuristic MatrixHeuristic;

/**
 * A* from @src to @dst.  @h may be NULL (then this is Dijkstra).  Weights must
 * be non-negative.  On success @out (allocated; matrix_path_free) holds the
 * path.  GRAPH_ERR_NOT_FOUND if @dst is unreachable.  An admissible heuristic
 * gives an optimal path even if it is not consistent.
 */
MatrixStatus matrix_a_star(const MatrixGraph *g, int src, int dst, MatrixHeuristic h, void *ctx,
                           MatrixPath *out);

/**
 * Dijkstra from several sources at once: distance to the nearest source.
 * @origin (optional, N entries) = which source that is.  Non-negative weights.
 */
MatrixStatus matrix_dijkstra_multi(const MatrixGraph *g, const int *sources, size_t count,
                                   double *distance, int *parent, int *origin);

/** Shortest paths when every weight is 0 or 1 (else INVALID_ARGUMENT).  O(V + E). */
MatrixStatus matrix_zero_one_bfs(const MatrixGraph *g, int src, double *distance, int *parent);

/**
 * Up to @k shortest loopless src -> dst paths in non-decreasing cost (Yen).
 * @paths: k entries, the first *found of which are filled (each freed with
 * matrix_path_free).  Non-negative weights.  GRAPH_ERR_NOT_FOUND if no path.
 */
MatrixStatus matrix_k_shortest_paths(const MatrixGraph *g, int src, int dst, size_t k,
                                     MatrixPath *paths, size_t *found);

/* ================================================================== */
/*  4. Spanning trees, flow, matching                                   */
/* ================================================================== */

/**
 * Prim's MST on the dense matrix, O(N^2).  Same contract as matrix_mst_kruskal:
 * @edges needs vertex_count - 1 entries; returns the total weight, INFINITY if
 * the graph is directed or disconnected.
 */
double matrix_mst_prim(const MatrixGraph *g, MatrixEdge *edges, size_t *edge_count);

typedef GraphFlowAlgo MatrixFlowAlgo;
#define MATRIX_FLOW_DINIC GRAPH_FLOW_DINIC
#define MATRIX_FLOW_EDMONDS_KARP GRAPH_FLOW_EDMONDS_KARP

/**
 * Maximum flow from @source to @sink; edge weights are capacities (must be
 * >= 0).  An undirected edge carries flow either way.  *value receives the
 * flow; @source_side (optional, N entries) marks the source side of a minimum
 * cut.  Dinic: O(V^2 E); Edmonds-Karp: O(V E^2).
 */
MatrixStatus matrix_max_flow(const MatrixGraph *g, int source, int sink, MatrixFlowAlgo algo,
                             double *value, bool *source_side);

/**
 * Minimum source-sink cut.  *value = cut capacity; @cut (needs edge_count
 * entries) lists the edges leaving the source side; *count = how many.
 */
MatrixStatus matrix_min_cut(const MatrixGraph *g, int source, int sink, double *value,
                            MatrixEdge *cut, size_t *count);

/**
 * Maximum cardinality matching in a bipartite graph (Hopcroft-Karp), O(E sqrt(V)).
 * match: N entries, match[v] = partner or -1.  *size = number of matched pairs.
 * GRAPH_ERR_UNSUPPORTED if the graph is not bipartite.
 */
MatrixStatus matrix_bipartite_matching(const MatrixGraph *g, int *match, size_t *size);

/**
 * Minimum-cost assignment (Hungarian), O(rows^2 * cols).  @cost is row-major
 * rows x cols with rows <= cols; assignment[r] = the column given to row r.
 */
MatrixStatus matrix_hungarian(const double *cost, size_t rows, size_t cols, int *assignment,
                              double *total);

/* ================================================================== */
/*  5. Combinatorial                                                    */
/* ================================================================== */

/**
 * Euler trail that uses every edge once.  If @closed, a circuit is required.
 * @path needs edge_count + 1 entries; *length = edge_count + 1 (0 for a graph
 * with no edges).  GRAPH_ERR_NOT_FOUND if none exists.  Complexity: O(V + E).
 */
MatrixStatus matrix_euler_trail(const MatrixGraph *g, bool closed, int *path, size_t *length);

/**
 * Hamiltonian path (or cycle) by bitmask DP.  At most 24 alive vertices
 * (else TOO_LARGE).  @path needs vertex_count entries.  NOT_FOUND if none.
 * A cycle needs at least 3 vertices (2 if directed).  Complexity: O(2^V V^2).
 */
MatrixStatus matrix_hamiltonian(const MatrixGraph *g, bool cycle, int *path);

/**
 * Proper vertex colouring of the underlying undirected graph; colours start at
 * 0.  color: N entries (-1 for dead).  *ncolors (optional) = colours used.
 * Greedy visits vertices in ascending ID; Welsh-Powell by descending degree;
 * DSatur by saturation (usually fewest colours).
 */
MatrixStatus matrix_color_greedy(const MatrixGraph *g, int *color, size_t *ncolors);
MatrixStatus matrix_color_welsh_powell(const MatrixGraph *g, int *color, size_t *ncolors);
MatrixStatus matrix_color_dsatur(const MatrixGraph *g, int *color, size_t *ncolors);

/**
 * All maximal cliques (Bron-Kerbosch with pivoting), undirected only.  Output
 * size can be exponential.  Free with matrix_groups_free.
 */
MatrixStatus matrix_maximal_cliques(const MatrixGraph *g, MatrixGroups *out);

/** A largest clique.  @clique needs vertex_count entries; *size = its size. */
MatrixStatus matrix_max_clique(const MatrixGraph *g, int *clique, size_t *size);

/** 2-approximate vertex cover (maximal matching).  in_cover: N entries. */
MatrixStatus matrix_vertex_cover_approx(const MatrixGraph *g, bool *in_cover, size_t *size);

/** Greedy maximal independent set (lowest degree first).  in_set: N entries. */
MatrixStatus matrix_independent_set_greedy(const MatrixGraph *g, bool *in_set, size_t *size);

/**
 * Longest weighted path in a DAG.  distance[v] = longest path ending at v
 * (>= 0, a single vertex is a path; dead IDs get -INFINITY); parent: N entries;
 * *end = the vertex with the largest distance (-1 if the graph is empty).
 * GRAPH_ERR_UNSUPPORTED if undirected or cyclic.
 */
MatrixStatus matrix_dag_longest_path(const MatrixGraph *g, double *distance, int *parent, int *end);

/** The critical path itself.  @path needs vertex_count entries. */
MatrixStatus matrix_critical_path(const MatrixGraph *g, int *path, size_t *length, double *cost);

/* ================================================================== */
/*  6. Network analysis                                                 */
/* ================================================================== */

/**
 * PageRank by power iteration.  Edge weights (>= 0) are transition weights;
 * dangling vertices spread rank evenly.  out sums to 1 over alive vertices.
 * *iterations (optional) = rounds used.  GRAPH_ERR_NO_CONVERGENCE if the L1
 * change stays above @tolerance after @max_iterations (out still holds the
 * latest estimate).
 */
MatrixStatus matrix_pagerank(const MatrixGraph *g, double damping, double tolerance,
                             size_t max_iterations, double *out, size_t *iterations);

/** HITS hub and authority scores, each L2-normalised. */
MatrixStatus matrix_hits(const MatrixGraph *g, double tolerance, size_t max_iterations,
                         double *hub, double *authority);

/** Eigenvector centrality (in-edges for directed graphs), L2-normalised. */
MatrixStatus matrix_eigenvector_centrality(const MatrixGraph *g, double tolerance,
                                           size_t max_iterations, double *out);

/** Core number of each vertex (k-core decomposition) of the underlying undirected graph. */
MatrixStatus matrix_core_numbers(const MatrixGraph *g, int *core);

/**
 * Community detection by label propagation, deterministic (ascending vertex
 * order, ties to the smaller label).  community: N entries numbered 0.. by
 * first appearance; *count (optional) = number of communities.
 */
MatrixStatus matrix_label_propagation(const MatrixGraph *g, int *community, size_t *count);

/**
 * Louvain modularity optimisation on the underlying undirected graph (edge
 * weights must be >= 0; self-loops are ignored).  community: N entries (-1 for
 * dead IDs); *modularity (optional) and *count (optional) describe the result.
 * Deterministic: vertices are visited in ascending ID order.
 */
MatrixStatus matrix_louvain(const MatrixGraph *g, int *community, double *modularity, size_t *count);

/** Modularity of a given partition (community has N entries, -1 for dead IDs). */
MatrixStatus matrix_modularity(const MatrixGraph *g, const int *community, double *q);

/** Newman's degree assortativity.  UNSUPPORTED if undefined (e.g. a regular graph). */
MatrixStatus matrix_degree_assortativity(const MatrixGraph *g, double *r);

/* ================================================================== */
/*  7. Graph operations                                                 */
/* ================================================================== */

/** Deep copy with identical IDs.  NULL on failure. */
MatrixGraph *matrix_graph_copy(const MatrixGraph *g);

/**
 * Subgraph induced by vertices[0..count-1] (all alive, no repeats).  Vertex
 * vertices[i] becomes ID i in the result.
 */
MatrixStatus matrix_induced_subgraph(const MatrixGraph *g, const int *vertices, size_t count,
                                     MatrixGraph **out);

/**
 * Union / intersection of two graphs with the same directedness over a shared
 * ID space.  Where both have an edge, @a's weight wins.
 */
MatrixStatus matrix_graph_union(const MatrixGraph *a, const MatrixGraph *b, MatrixGraph **out);
MatrixStatus matrix_graph_intersection(const MatrixGraph *a, const MatrixGraph *b, MatrixGraph **out);

/** Edges that are absent become edges of weight 1.0; no self-loops; same alive vertices. */
MatrixStatus matrix_graph_complement(const MatrixGraph *g, MatrixGraph **out);

/**
 * Structural isomorphism (weights ignored).  Backtracking with degree pruning:
 * fast on ordinary graphs, exponential in the worst case, so meant for small
 * graphs.  mapping (optional, N entries of a) receives a vertex of b per vertex of a.
 */
MatrixStatus matrix_is_isomorphic(const MatrixGraph *a, const MatrixGraph *b, bool *result,
                                  int *mapping);

/* ================================================================== */
/*  8. Generators (deterministic for a given seed)                      */
/* ================================================================== */

/** G(n, p): each possible edge present with probability p, weight uniform in [wmin, wmax]. */
MatrixGraph *matrix_gen_erdos_renyi(size_t n, double p, unsigned long long seed, bool directed,
                                    double wmin, double wmax);
/**
 * Preferential attachment: a complete graph on m + 1 vertices, then each new
 * vertex attaches to m distinct existing ones with probability proportional to
 * their degree.  Undirected.  n <= m gives a complete graph.
 */
MatrixGraph *matrix_gen_barabasi_albert(size_t n, size_t m, unsigned long long seed);
MatrixGraph *matrix_gen_grid(size_t rows, size_t cols);                 /* 4-neighbour lattice */
MatrixGraph *matrix_gen_path(size_t n, bool directed);
/** Undirected rings need n >= 3 (smaller n is just a path); a directed ring of 1 is a self-loop. */
MatrixGraph *matrix_gen_cycle(size_t n, bool directed);
MatrixGraph *matrix_gen_star(size_t n);                                 /* vertex 0 is the hub  */
MatrixGraph *matrix_gen_complete(size_t n, bool directed);
/** Edges only from lower to higher ID, so it is always a DAG. */
MatrixGraph *matrix_gen_random_dag(size_t n, double p, unsigned long long seed,
                                   double wmin, double wmax);
/** Uniformly attached random tree (vertex i joins a random earlier vertex). */
MatrixGraph *matrix_gen_random_tree(size_t n, unsigned long long seed);

/* ================================================================== */
/*  9. Input / output                                                   */
/* ================================================================== */

/**
 * Graphviz DOT.  @labels (optional, N entries, may contain NULLs) names the
 * vertices; @show_weights prints edge weights.  Dead IDs are omitted.
 */
MatrixStatus matrix_write_dot(const MatrixGraph *g, FILE *out, const char *name,
                              const char *const *labels, bool show_weights);

/**
 * Edge list, one `u v [weight]` per line; blank lines and lines starting with
 * `#` are ignored; a missing weight is 1.0.  Vertex IDs are non-negative
 * integers and the graph grows to hold the largest one.  The first line may be
 * `# directed` or `# undirected` to override @directed.
 */
MatrixStatus matrix_read_edge_list(FILE *in, bool directed, MatrixGraph **out);
MatrixStatus matrix_write_edge_list(const MatrixGraph *g, FILE *out);

/**
 * Edge list with arbitrary text names (`alice bob 2.5`).  Names get IDs 0.. in
 * order of first appearance.  *names_out (N strings, each freed individually
 * then the array, or all at once with matrix_names_free) maps ID -> name.
 */
MatrixStatus matrix_read_named_edge_list(FILE *in, bool directed, MatrixGraph **out,
                                         char ***names_out, size_t *count);
void matrix_names_free(char **names, size_t count);

/** Compact binary snapshot: keeps IDs, dead IDs, direction and weights exactly. */
MatrixStatus matrix_save_binary(const MatrixGraph *g, FILE *out);
MatrixStatus matrix_load_binary(FILE *in, MatrixGraph **out);

/** Export-only formats. */
MatrixStatus matrix_write_graphml(const MatrixGraph *g, FILE *out);
MatrixStatus matrix_write_gml(const MatrixGraph *g, FILE *out);
MatrixStatus matrix_write_json(const MatrixGraph *g, FILE *out);

#endif /* MATRIX_ALGORITHMS_H */
