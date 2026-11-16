/*
 * list_advanced.h — declarations for list_advanced.c
 *
 * Twenty-five additional graph algorithms built on the adjacency-list
 * backend (adjacency_list.h) and the helpers in list_util.h.  Every
 * public name is prefixed `list_adv_` / `ListAdv` so it cannot collide
 * with anything in list_algorithms.h.
 *
 * Conventions (same as list_analysis.c):
 *   - functions return a ListStatus (GRAPH_OK on success);
 *   - result arrays are supplied by the caller and sized by
 *     list_graph_id_bound(); dead IDs get a documented sentinel;
 *   - the few results that own heap memory have an explicit *_free function.
 *
 * Compile: cc -std=c11 -Wall -Wextra -Werror -pedantic -c list_advanced.c
 */

#ifndef LIST_ADVANCED_H
#define LIST_ADVANCED_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "adjacency_list.h"
#include "graph_results.h"
#include "graph_status.h"
#include "list_algorithms.h"

/* ------------------------------------------------------------------ */
/*  Result types                                                        */
/* ------------------------------------------------------------------ */

/** A vertex sequence with its total weight.  Release with list_adv_path_free. */
typedef GraphPath ListAdvPath;

typedef GraphAdvCycles ListAdvCycles;

typedef GraphAdvTreeOrder ListAdvTreeOrder;
#define LIST_ADV_TREE_PREORDER GRAPH_ADV_TREE_PREORDER
#define LIST_ADV_TREE_INORDER GRAPH_ADV_TREE_INORDER
#define LIST_ADV_TREE_POSTORDER GRAPH_ADV_TREE_POSTORDER
#define LIST_ADV_TREE_LEVELORDER GRAPH_ADV_TREE_LEVELORDER

typedef GraphAdvCycleKind ListAdvCycleKind;
#define LIST_ADV_CYCLE_NONE GRAPH_ADV_CYCLE_NONE
#define LIST_ADV_CYCLE_ZERO GRAPH_ADV_CYCLE_ZERO
#define LIST_ADV_CYCLE_NEGATIVE GRAPH_ADV_CYCLE_NEGATIVE

typedef GraphAdvCentralityKind ListAdvCentralityKind;
#define LIST_ADV_CENTRALITY_DEGREE GRAPH_ADV_CENTRALITY_DEGREE
#define LIST_ADV_CENTRALITY_CLOSENESS GRAPH_ADV_CENTRALITY_CLOSENESS
#define LIST_ADV_CENTRALITY_BETWEENNESS GRAPH_ADV_CENTRALITY_BETWEENNESS
#define LIST_ADV_CENTRALITY_HARMONIC GRAPH_ADV_CENTRALITY_HARMONIC
#define LIST_ADV_CENTRALITY_EIGENVECTOR GRAPH_ADV_CENTRALITY_EIGENVECTOR

typedef GraphAdvActivation ListAdvActivation;
#define LIST_ADV_ACT_NONE GRAPH_ADV_ACT_NONE
#define LIST_ADV_ACT_RELU GRAPH_ADV_ACT_RELU
#define LIST_ADV_ACT_TANH GRAPH_ADV_ACT_TANH

/** One graph-convolution layer: out = act(A_norm * in * W + bias). */
typedef GraphAdvGnnLayer ListAdvGnnLayer;

/** Heuristic estimate of the remaining cost from v to goal; must be >= 0. */
typedef GraphHeuristic ListAdvHeuristic;

/* ------------------------------------------------------------------ */
/*  Helpers                                                             */
/* ------------------------------------------------------------------ */

void         list_adv_path_free(ListAdvPath *path);
void         list_adv_cycles_free(ListAdvCycles *cycles);
ListStatus list_adv_path_from_parents(const int *parent, size_t n, int source, int target,
                                          double cost, ListAdvPath *out);
ListStatus list_adv_apsp_path(const double *dist, const int *pred, size_t n,
                                  int source, int target, ListAdvPath *out);

/* ------------------------------------------------------------------ */
/*  The algorithms                                                      */
/* ------------------------------------------------------------------ */

/*  1 */ ListStatus list_adv_longest_path(const ListGraph *g, int source,
                                              double *dist, int *parent);
/*  2 */ ListStatus list_adv_pagerank(const ListGraph *g, double damping, double tol,
                                          size_t max_iter, bool weighted,
                                          double *rank, size_t *iterations);
/*  3 */ ListStatus list_adv_reverse_delete(const ListGraph *g, ListEdge *edges,
                                                size_t *edge_count, double *total_weight,
                                                size_t *components);
/*  4 */ ListStatus list_adv_rocha_thatte(const ListGraph *g, size_t max_len,
                                              size_t max_cycles, ListAdvCycles *out);
/*  5 */ ListStatus list_adv_tree_traversal(const ListGraph *g, int root,
                                                ListAdvTreeOrder order, int *out, size_t *count);
/*  6 */ ListStatus list_adv_zero_weight_cycle(const ListGraph *g, double eps,
                                                   ListAdvCycleKind *kind,
                                                   int *cycle, size_t *cycle_len);
/*  7 */ ListStatus list_adv_subgraph_isomorphism(const ListGraph *pattern,
                                                      const ListGraph *target, bool induced,
                                                      bool match_weights, int *mapping,
                                                      bool *found);
/*  8 */ ListStatus list_adv_coffman_graham(const ListGraph *g, size_t processors,
                                                int *label, int *slot, size_t *makespan);
/*  9 */ ListStatus list_adv_crossing_number_heuristic(const ListGraph *g,
                                                           size_t max_passes,
                                                           long long *crossings, int *order);
/* 10 */ ListStatus list_adv_color_refinement(const ListGraph *g, int *color,
                                                  size_t *num_colors, size_t *rounds);
/* 11 */ ListStatus list_adv_color_coding_path(const ListGraph *g, size_t k,
                                                   size_t trials, uint64_t seed,
                                                   bool *found, int *path);
/* 12 */ ListStatus list_adv_centrality(const ListGraph *g, ListAdvCentralityKind kind,
                                            double *out);
         ListStatus list_adv_harmonic_centrality(const ListGraph *g, double *out);
         ListStatus list_adv_eigenvector_centrality(const ListGraph *g, size_t max_iter,
                                                        double tol, double *out,
                                                        size_t *iterations);
/* 13 */ ListStatus list_adv_astar(const ListGraph *g, int source, int target,
                                       ListAdvHeuristic h, void *ctx,
                                       ListAdvPath *path, size_t *expanded);
/* 14 */ ListStatus list_adv_flood(const ListGraph *g, int source, size_t ttl,
                                       int *round, int *parent, size_t *messages,
                                       size_t *reached);
/* 15 */ ListStatus list_adv_force_layout(const ListGraph *g, size_t iterations,
                                              double *x, double *y);
/* 16 */ ListStatus list_adv_floyd_warshall(const ListGraph *g, double *dist, int *pred);
/* 17 */ ListStatus list_adv_gnn_forward(const ListGraph *g, const double *features,
                                             size_t in_dim, const ListAdvGnnLayer *layers,
                                             size_t num_layers, double *output);
/* 18 */ ListStatus list_adv_graph_isomorphism(const ListGraph *a, const ListGraph *b,
                                                   bool match_weights, bool *isomorphic,
                                                   int *mapping);
/* 19 */ ListStatus list_adv_ida_star(const ListGraph *g, int source, int target,
                                          ListAdvHeuristic h, void *ctx, size_t node_limit,
                                          ListAdvPath *path);
/* 20 */ ListStatus list_adv_iddfs(const ListGraph *g, int source, int target,
                                       size_t max_depth, ListAdvPath *path);
/* 21 */ ListStatus list_adv_jump_point_search(const ListGraph *g, size_t width,
                                                   int source, int target, ListAdvPath *path);
/* 22 */ ListStatus list_adv_johnson(const ListGraph *g, double *dist, int *pred);
/* 23 */ ListStatus list_adv_k_shortest_paths(const ListGraph *g, int source, int target,
                                                  size_t k, ListAdvPath *paths, size_t *found);
/* 24 */ ListStatus list_adv_widest_path(const ListGraph *g, int source,
                                             double *width, int *parent);
/* 25 */ ListStatus list_adv_lex_bfs(const ListGraph *g, int start, int *order);

/* ------------------------------------------------------------------ */
/*  Additions: drawing metrics and dynamic connectivity                 */
/* ------------------------------------------------------------------ */

/** Incremental/decremental connectivity tracker.  Create with _create, free with _destroy. */
typedef struct ListAdvDynConn ListAdvDynConn;

/* 26 */ ListStatus list_adv_crossing_count(const ListGraph *g, const int *order,
                                                size_t order_len, long long *crossings);
/* 27 */ ListStatus list_adv_circle_layout(const ListGraph *g, double radius,
                                               double *x, double *y);
/* 28 */ ListAdvDynConn *list_adv_dynconn_create(const ListGraph *g);
         void              list_adv_dynconn_destroy(ListAdvDynConn *dc);
         ListStatus      list_adv_dynconn_connected(ListAdvDynConn *dc, int u, int v,
                                                        bool *connected);
         ListStatus      list_adv_dynconn_components(ListAdvDynConn *dc, size_t *count);

#endif /* LIST_ADVANCED_H */
