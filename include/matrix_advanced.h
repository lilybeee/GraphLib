/*
 * matrix_advanced.h — declarations for matrix_advanced.c
 *
 * Twenty-five additional graph algorithms built on the adjacency-matrix
 * backend (adjacency_matrix.h) and the helpers in matrix_util.h.  Every
 * public name is prefixed `matrix_adv_` / `MatrixAdv` so it cannot collide
 * with anything in matrix_algorithms.h.
 *
 * Conventions (same as matrix_analysis.c):
 *   - functions return a MatrixStatus (GRAPH_OK on success);
 *   - result arrays are supplied by the caller and sized by
 *     matrix_graph_id_bound(); dead IDs get a documented sentinel;
 *   - the few results that own heap memory have an explicit *_free function.
 *
 * Compile: cc -std=c11 -Wall -Wextra -Werror -pedantic -c matrix_advanced.c
 */

#ifndef MATRIX_ADVANCED_H
#define MATRIX_ADVANCED_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "adjacency_matrix.h"
#include "graph_results.h"
#include "graph_status.h"
#include "matrix_algorithms.h"

/* ------------------------------------------------------------------ */
/*  Result types                                                        */
/* ------------------------------------------------------------------ */

/** A vertex sequence with its total weight.  Release with matrix_adv_path_free. */
typedef GraphPath MatrixAdvPath;

typedef GraphAdvCycles MatrixAdvCycles;

typedef GraphAdvTreeOrder MatrixAdvTreeOrder;
#define MATRIX_ADV_TREE_PREORDER GRAPH_ADV_TREE_PREORDER
#define MATRIX_ADV_TREE_INORDER GRAPH_ADV_TREE_INORDER
#define MATRIX_ADV_TREE_POSTORDER GRAPH_ADV_TREE_POSTORDER
#define MATRIX_ADV_TREE_LEVELORDER GRAPH_ADV_TREE_LEVELORDER

typedef GraphAdvCycleKind MatrixAdvCycleKind;
#define MATRIX_ADV_CYCLE_NONE GRAPH_ADV_CYCLE_NONE
#define MATRIX_ADV_CYCLE_ZERO GRAPH_ADV_CYCLE_ZERO
#define MATRIX_ADV_CYCLE_NEGATIVE GRAPH_ADV_CYCLE_NEGATIVE

typedef GraphAdvCentralityKind MatrixAdvCentralityKind;
#define MATRIX_ADV_CENTRALITY_DEGREE GRAPH_ADV_CENTRALITY_DEGREE
#define MATRIX_ADV_CENTRALITY_CLOSENESS GRAPH_ADV_CENTRALITY_CLOSENESS
#define MATRIX_ADV_CENTRALITY_BETWEENNESS GRAPH_ADV_CENTRALITY_BETWEENNESS
#define MATRIX_ADV_CENTRALITY_HARMONIC GRAPH_ADV_CENTRALITY_HARMONIC
#define MATRIX_ADV_CENTRALITY_EIGENVECTOR GRAPH_ADV_CENTRALITY_EIGENVECTOR

typedef GraphAdvActivation MatrixAdvActivation;
#define MATRIX_ADV_ACT_NONE GRAPH_ADV_ACT_NONE
#define MATRIX_ADV_ACT_RELU GRAPH_ADV_ACT_RELU
#define MATRIX_ADV_ACT_TANH GRAPH_ADV_ACT_TANH

/** One graph-convolution layer: out = act(A_norm * in * W + bias). */
typedef GraphAdvGnnLayer MatrixAdvGnnLayer;

/** Heuristic estimate of the remaining cost from v to goal; must be >= 0. */
typedef GraphHeuristic MatrixAdvHeuristic;

/* ------------------------------------------------------------------ */
/*  Helpers                                                             */
/* ------------------------------------------------------------------ */

void         matrix_adv_path_free(MatrixAdvPath *path);
void         matrix_adv_cycles_free(MatrixAdvCycles *cycles);
MatrixStatus matrix_adv_path_from_parents(const int *parent, size_t n, int source, int target,
                                          double cost, MatrixAdvPath *out);
MatrixStatus matrix_adv_apsp_path(const double *dist, const int *pred, size_t n,
                                  int source, int target, MatrixAdvPath *out);

/* ------------------------------------------------------------------ */
/*  The algorithms                                                      */
/* ------------------------------------------------------------------ */

/*  1 */ MatrixStatus matrix_adv_longest_path(const MatrixGraph *g, int source,
                                              double *dist, int *parent);
/*  2 */ MatrixStatus matrix_adv_pagerank(const MatrixGraph *g, double damping, double tol,
                                          size_t max_iter, bool weighted,
                                          double *rank, size_t *iterations);
/*  3 */ MatrixStatus matrix_adv_reverse_delete(const MatrixGraph *g, MatrixEdge *edges,
                                                size_t *edge_count, double *total_weight,
                                                size_t *components);
/*  4 */ MatrixStatus matrix_adv_rocha_thatte(const MatrixGraph *g, size_t max_len,
                                              size_t max_cycles, MatrixAdvCycles *out);
/*  5 */ MatrixStatus matrix_adv_tree_traversal(const MatrixGraph *g, int root,
                                                MatrixAdvTreeOrder order, int *out, size_t *count);
/*  6 */ MatrixStatus matrix_adv_zero_weight_cycle(const MatrixGraph *g, double eps,
                                                   MatrixAdvCycleKind *kind,
                                                   int *cycle, size_t *cycle_len);
/*  7 */ MatrixStatus matrix_adv_subgraph_isomorphism(const MatrixGraph *pattern,
                                                      const MatrixGraph *target, bool induced,
                                                      bool match_weights, int *mapping,
                                                      bool *found);
/*  8 */ MatrixStatus matrix_adv_coffman_graham(const MatrixGraph *g, size_t processors,
                                                int *label, int *slot, size_t *makespan);
/*  9 */ MatrixStatus matrix_adv_crossing_number_heuristic(const MatrixGraph *g,
                                                           size_t max_passes,
                                                           long long *crossings, int *order);
/* 10 */ MatrixStatus matrix_adv_color_refinement(const MatrixGraph *g, int *color,
                                                  size_t *num_colors, size_t *rounds);
/* 11 */ MatrixStatus matrix_adv_color_coding_path(const MatrixGraph *g, size_t k,
                                                   size_t trials, uint64_t seed,
                                                   bool *found, int *path);
/* 12 */ MatrixStatus matrix_adv_centrality(const MatrixGraph *g, MatrixAdvCentralityKind kind,
                                            double *out);
         MatrixStatus matrix_adv_harmonic_centrality(const MatrixGraph *g, double *out);
         MatrixStatus matrix_adv_eigenvector_centrality(const MatrixGraph *g, size_t max_iter,
                                                        double tol, double *out,
                                                        size_t *iterations);
/* 13 */ MatrixStatus matrix_adv_astar(const MatrixGraph *g, int source, int target,
                                       MatrixAdvHeuristic h, void *ctx,
                                       MatrixAdvPath *path, size_t *expanded);
/* 14 */ MatrixStatus matrix_adv_flood(const MatrixGraph *g, int source, size_t ttl,
                                       int *round, int *parent, size_t *messages,
                                       size_t *reached);
/* 15 */ MatrixStatus matrix_adv_force_layout(const MatrixGraph *g, size_t iterations,
                                              double *x, double *y);
/* 16 */ MatrixStatus matrix_adv_floyd_warshall(const MatrixGraph *g, double *dist, int *pred);
/* 17 */ MatrixStatus matrix_adv_gnn_forward(const MatrixGraph *g, const double *features,
                                             size_t in_dim, const MatrixAdvGnnLayer *layers,
                                             size_t num_layers, double *output);
/* 18 */ MatrixStatus matrix_adv_graph_isomorphism(const MatrixGraph *a, const MatrixGraph *b,
                                                   bool match_weights, bool *isomorphic,
                                                   int *mapping);
/* 19 */ MatrixStatus matrix_adv_ida_star(const MatrixGraph *g, int source, int target,
                                          MatrixAdvHeuristic h, void *ctx, size_t node_limit,
                                          MatrixAdvPath *path);
/* 20 */ MatrixStatus matrix_adv_iddfs(const MatrixGraph *g, int source, int target,
                                       size_t max_depth, MatrixAdvPath *path);
/* 21 */ MatrixStatus matrix_adv_jump_point_search(const MatrixGraph *g, size_t width,
                                                   int source, int target, MatrixAdvPath *path);
/* 22 */ MatrixStatus matrix_adv_johnson(const MatrixGraph *g, double *dist, int *pred);
/* 23 */ MatrixStatus matrix_adv_k_shortest_paths(const MatrixGraph *g, int source, int target,
                                                  size_t k, MatrixAdvPath *paths, size_t *found);
/* 24 */ MatrixStatus matrix_adv_widest_path(const MatrixGraph *g, int source,
                                             double *width, int *parent);
/* 25 */ MatrixStatus matrix_adv_lex_bfs(const MatrixGraph *g, int start, int *order);

/* ------------------------------------------------------------------ */
/*  Additions: drawing metrics and dynamic connectivity                 */
/* ------------------------------------------------------------------ */

/** Incremental/decremental connectivity tracker.  Create with _create, free with _destroy. */
typedef struct MatrixAdvDynConn MatrixAdvDynConn;

/* 26 */ MatrixStatus matrix_adv_crossing_count(const MatrixGraph *g, const int *order,
                                                size_t order_len, long long *crossings);
/* 27 */ MatrixStatus matrix_adv_circle_layout(const MatrixGraph *g, double radius,
                                               double *x, double *y);
/* 28 */ MatrixAdvDynConn *matrix_adv_dynconn_create(const MatrixGraph *g);
         void              matrix_adv_dynconn_destroy(MatrixAdvDynConn *dc);
         MatrixStatus      matrix_adv_dynconn_connected(MatrixAdvDynConn *dc, int u, int v,
                                                        bool *connected);
         MatrixStatus      matrix_adv_dynconn_components(MatrixAdvDynConn *dc, size_t *count);

#endif /* MATRIX_ADVANCED_H */
