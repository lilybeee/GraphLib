/*
 * graph_results.h — result and callback types of the algorithm library, shared
 * by every representation.
 *
 * The list and matrix libraries each used to define their own copies; the
 * unified Graph facade (graph.h) needs one type that works for both, so they
 * live here and list_algorithms.h / matrix_algorithms.h alias them
 * (ListGroups == MatrixGroups == GraphGroups, and so on).
 */

#ifndef GRAPH_RESULTS_H
#define GRAPH_RESULTS_H

#include <stdbool.h>
#include <stddef.h>

/** A list of vertex sets (components, cliques, ...).  Free with graph_groups_free(). */
typedef struct {
    size_t  count;     /* number of groups                                        */
    size_t *offsets;   /* count + 1 entries; group i = members[offsets[i] .. offsets[i+1]) */
    int    *members;
} GraphGroups;

/** A path with its total weight.  Free with graph_path_free(). */
typedef struct {
    int    *vertices;
    size_t  length;    /* number of vertices */
    double  cost;
} GraphPath;

typedef struct {
    double density;      /* edges / possible edges (self-loops excluded)     */
    double reciprocity;  /* fraction of arcs with a reverse arc; 1.0 if undirected */
    size_t max_degree;
} GraphNetworkStats;

/** Called once per vertex as it is visited; return false to stop early. */
typedef bool (*GraphVisitFn)(int vertex, int parent, void *ctx);

/** Estimate of the remaining cost from @vertex to @target; must never overestimate. */
typedef double (*GraphHeuristic)(int vertex, int target, void *ctx);

typedef enum { GRAPH_FLOW_DINIC, GRAPH_FLOW_EDMONDS_KARP } GraphFlowAlgo;

/* ---- Types shared by the advanced algorithms (list_advanced / matrix_advanced /
 *      graph_advanced), so one result type serves both representations. ---- */

/** A list of cycles, each stored as a run of vertices in `members`.  Free with graph_adv_cycles_free(). */
typedef struct {
    size_t  count;
    size_t *offsets;    /* count + 1 entries; cycle i is members[offsets[i] .. offsets[i+1]) */
    int    *members;
    bool    truncated;  /* true if the search stopped at max_cycles */
} GraphAdvCycles;

typedef enum {
    GRAPH_ADV_TREE_PREORDER,
    GRAPH_ADV_TREE_INORDER,
    GRAPH_ADV_TREE_POSTORDER,
    GRAPH_ADV_TREE_LEVELORDER
} GraphAdvTreeOrder;

typedef enum {
    GRAPH_ADV_CYCLE_NONE,       /* no cycle of weight <= 0 */
    GRAPH_ADV_CYCLE_ZERO,       /* a cycle of weight 0 (within eps) */
    GRAPH_ADV_CYCLE_NEGATIVE    /* a negative cycle exists (witness returned) */
} GraphAdvCycleKind;

typedef enum {
    GRAPH_ADV_CENTRALITY_DEGREE,
    GRAPH_ADV_CENTRALITY_CLOSENESS,
    GRAPH_ADV_CENTRALITY_BETWEENNESS,
    GRAPH_ADV_CENTRALITY_HARMONIC,
    GRAPH_ADV_CENTRALITY_EIGENVECTOR
} GraphAdvCentralityKind;

typedef enum {
    GRAPH_ADV_ACT_NONE,
    GRAPH_ADV_ACT_RELU,
    GRAPH_ADV_ACT_TANH
} GraphAdvActivation;

/** One graph-convolution layer: out = act(A_norm * in * W + bias). */
typedef struct {
    size_t              in_dim;
    size_t              out_dim;
    const double       *weights;   /* in_dim x out_dim, row-major */
    const double       *bias;      /* out_dim entries, or NULL */
    GraphAdvActivation  activation;
} GraphAdvGnnLayer;

/** Frees the arrays inside @groups and zeroes it.  Safe with NULL. */
void graph_groups_free(GraphGroups *groups);

/** Frees the vertex array inside @path and zeroes it.  Safe with NULL. */
void graph_path_free(GraphPath *path);

#endif /* GRAPH_RESULTS_H */
