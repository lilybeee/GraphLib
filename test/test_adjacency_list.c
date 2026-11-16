/*
 * test_adjacency_list.c — unit tests for the adjacency-list graph backend.
 *
 * Each test mirrors the corresponding test in test_adjacency_matrix.c so that
 * both backends are validated against the same scenarios.
 *
 * Compile & run:
 *   cc -std=c11 -Wall -Wextra -Werror adjacency_list.c test_adjacency_list.c -lm -o list_test
 *   ./list_test
 */

#include "adjacency_list.h"
#include "backend.h" /* compiled here so the interface is checked by -Werror */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/*  1. Basic graph operations                                           */
/* ------------------------------------------------------------------ */

static void test_graph_operations(void)
{
    ListGraph *graph = list_graph_create(3, false);
    assert(graph != NULL);
    assert(list_graph_vertex_count(graph) == 3);
    assert(list_graph_edge_count(graph) == 0);
    assert(!list_graph_is_directed(graph));

    /* Valid edge */
    assert(list_graph_add_edge(graph, 0, 1, 2.5));
    /* Duplicate */
    assert(!list_graph_add_edge(graph, 0, 1, 9.0));
    /* Out-of-range */
    assert(!list_graph_add_edge(graph, -1, 1, 1.0));
    assert(!list_graph_add_edge(graph, 0, 3, 1.0));
    /* Non-finite weight */
    assert(!list_graph_add_edge(graph, 0, 2, INFINITY));
    assert(!list_graph_add_edge(graph, 0, 2, NAN));
    assert(list_graph_edge_count(graph) == 1);

    /* has_edge and symmetry for undirected */
    double weight = 0.0;
    assert(list_graph_has_edge(graph, 1, 0, &weight));
    assert(fabs(weight - 2.5) < 1e-9);

    /* Remove */
    assert(list_graph_remove_edge(graph, 1, 0));
    assert(!list_graph_has_edge(graph, 0, 1, NULL));
    assert(!list_graph_remove_edge(graph, 0, 1)); /* already gone */
    assert(list_graph_edge_count(graph) == 0);

    list_graph_destroy(graph);
}

/* ------------------------------------------------------------------ */
/*  2. Traversals and connected components                              */
/* ------------------------------------------------------------------ */

static void test_traversals_and_components(void)
{
    /*
     * Graph (undirected):
     *   0 -- 1
     *   |    |
     *   2 -- 3
     *   (isolated component: 4 -- 5)
     */
    ListGraph *graph = list_graph_create(6, false);
    assert(graph != NULL);
    assert(list_graph_add_edge(graph, 0, 1, 1.0));
    assert(list_graph_add_edge(graph, 0, 2, 1.0));
    assert(list_graph_add_edge(graph, 1, 3, 1.0));
    assert(list_graph_add_edge(graph, 2, 3, 1.0));
    assert(list_graph_add_edge(graph, 4, 5, 1.0));

    int order[6];
    int parent[6];

    /* BFS from 0 should visit exactly {0,1,2,3} */
    assert(list_bfs(graph, 0, order, parent) == 4);
    /* First vertex must be the start */
    assert(order[0] == 0);
    /* Parent of start is -1 */
    assert(parent[0] == -1);
    /* 1 and 2 are direct children of 0 */
    assert(parent[1] == 0);
    assert(parent[2] == 0);
    /* 3 is connected to both 1 and 2; either can be its BFS parent.
     * The exact parent depends on neighbour traversal order, which is
     * implementation-defined (we prepend edges, so order may vary). */
    assert(parent[3] == 1 || parent[3] == 2);

    /* DFS from 0 should also visit exactly {0,1,2,3} */
    assert(list_dfs(graph, 0, order, parent) == 4);
    assert(order[0] == 0);

    /* Connected components: 2 components */
    int component[6];
    assert(list_connected_components(graph, component) == 2);
    assert(component[0] == component[1]);
    assert(component[0] == component[2]);
    assert(component[0] == component[3]);
    assert(component[4] == component[5]);
    assert(component[0] != component[4]);

    /* The undirected 0-1-3-2-0 loop is a cycle */
    assert(list_has_cycle(graph));

    list_graph_destroy(graph);

    /* Weakly connected components of a directed graph */
    ListGraph *directed = list_graph_create(4, true);
    assert(directed != NULL);
    assert(list_graph_add_edge(directed, 0, 1, 1.0));
    assert(list_graph_add_edge(directed, 2, 1, 1.0));
    /* Vertices 3 is isolated; 0,1,2 share a weak component (0→1←2) */
    assert(list_connected_components(directed, component) == 2);
    list_graph_destroy(directed);
}

/* ------------------------------------------------------------------ */
/*  3. Cycle detection and topological sort                             */
/* ------------------------------------------------------------------ */

static void test_cycles_and_topological_sort(void)
{
    /*
     * DAG:
     *   5→2, 5→0, 4→0, 4→1, 2→3, 3→1
     */
    ListGraph *graph = list_graph_create(6, true);
    assert(graph != NULL);
    assert(list_graph_add_edge(graph, 5, 2, 1.0));
    assert(list_graph_add_edge(graph, 5, 0, 1.0));
    assert(list_graph_add_edge(graph, 4, 0, 1.0));
    assert(list_graph_add_edge(graph, 4, 1, 1.0));
    assert(list_graph_add_edge(graph, 2, 3, 1.0));
    assert(list_graph_add_edge(graph, 3, 1, 1.0));
    assert(!list_has_cycle(graph));

    int order[6], position[6];
    assert(list_topological_sort(graph, order));
    for (int i = 0; i < 6; ++i) {
        position[order[i]] = i;
    }
    /* Check all implied ordering constraints */
    assert(position[5] < position[2]);
    assert(position[5] < position[0]);
    assert(position[4] < position[0]);
    assert(position[4] < position[1]);
    assert(position[2] < position[3]);
    assert(position[3] < position[1]);

    /* Adding a back-edge creates a cycle */
    assert(list_graph_add_edge(graph, 1, 5, 1.0));
    assert(list_has_cycle(graph));
    assert(!list_topological_sort(graph, order));

    list_graph_destroy(graph);

    /* An undirected tree has no cycle */
    ListGraph *tree = list_graph_create(3, false);
    assert(tree != NULL);
    assert(list_graph_add_edge(tree, 0, 1, 1.0));
    assert(list_graph_add_edge(tree, 1, 2, 1.0));
    assert(!list_has_cycle(tree));
    /* Topological sort is only for directed graphs */
    assert(!list_topological_sort(tree, order));
    /* Adding the closing edge creates a cycle */
    assert(list_graph_add_edge(tree, 2, 0, 1.0));
    assert(list_has_cycle(tree));
    list_graph_destroy(tree);
}

/* ------------------------------------------------------------------ */
/*  4. Dijkstra shortest paths                                          */
/* ------------------------------------------------------------------ */

static void test_shortest_paths(void)
{
    /*
     * Directed weighted graph:
     *   0→1 (10), 0→2 (3), 2→1 (1), 1→3 (2), 2→3 (8)
     *
     *   Shortest paths from 0:
     *     0: 0
     *     1: 0→2→1 = 4
     *     2: 0→2   = 3
     *     3: 0→2→1→3 = 6
     *     4: ∞ (unreachable)
     */
    ListGraph *graph = list_graph_create(5, true);
    assert(graph != NULL);
    assert(list_graph_add_edge(graph, 0, 1, 10.0));
    assert(list_graph_add_edge(graph, 0, 2,  3.0));
    assert(list_graph_add_edge(graph, 2, 1,  1.0));
    assert(list_graph_add_edge(graph, 1, 3,  2.0));
    assert(list_graph_add_edge(graph, 2, 3,  8.0));

    double distance[5];
    int    par[5];
    assert(list_dijkstra(graph, 0, distance, par));
    assert(fabs(distance[0])       < 1e-9);
    assert(fabs(distance[1] - 4.0) < 1e-9);
    assert(fabs(distance[2] - 3.0) < 1e-9);
    assert(fabs(distance[3] - 6.0) < 1e-9);
    assert(isinf(distance[4]));
    assert(par[0] == -1);
    assert(par[1] == 2);
    assert(par[2] == 0);
    assert(par[3] == 1);
    assert(par[4] == -1);

    list_graph_destroy(graph);

    /* Negative weight → Dijkstra must refuse */
    ListGraph *negative = list_graph_create(2, true);
    assert(negative != NULL);
    assert(list_graph_add_edge(negative, 0, 1, -1.0));
    assert(!list_dijkstra(negative, 0, distance, par));
    list_graph_destroy(negative);
}

/* ------------------------------------------------------------------ */
/*  5. Kruskal's MST                                                    */
/* ------------------------------------------------------------------ */

static void test_minimum_spanning_tree(void)
{
    /*
     * Undirected weighted graph (4 vertices):
     *   0-1 (10), 0-2 (6), 0-3 (5), 1-3 (15), 2-3 (4)
     *
     *   MST edges: 2-3 (4), 0-3 (5), 0-1 (10) → total = 19
     */
    ListGraph *graph = list_graph_create(4, false);
    assert(graph != NULL);
    assert(list_graph_add_edge(graph, 0, 1, 10.0));
    assert(list_graph_add_edge(graph, 0, 2,  6.0));
    assert(list_graph_add_edge(graph, 0, 3,  5.0));
    assert(list_graph_add_edge(graph, 1, 3, 15.0));
    assert(list_graph_add_edge(graph, 2, 3,  4.0));

    ListEdge mst[3];
    size_t   edge_count = 0;
    double   total = list_mst_kruskal(graph, mst, &edge_count);
    assert(edge_count == 3);
    assert(fabs(total - 19.0) < 1e-9);

    list_graph_destroy(graph);

    /* Disconnected graph: no spanning tree */
    ListGraph *disconnected = list_graph_create(3, false);
    assert(disconnected != NULL);
    assert(list_graph_add_edge(disconnected, 0, 1, 1.0));
    /* vertex 2 is isolated */
    assert(isinf(list_mst_kruskal(disconnected, mst, &edge_count)));
    assert(edge_count == 1);
    list_graph_destroy(disconnected);
}

/* ------------------------------------------------------------------ */
/*  6. Edge cases                                                       */
/* ------------------------------------------------------------------ */

static void test_edge_cases(void)
{
    /* Empty graph */
    ListGraph *empty = list_graph_create(0, true);
    assert(empty != NULL);
    assert(list_graph_vertex_count(empty) == 0);
    assert(list_connected_components(empty, NULL) == 0);
    assert(list_topological_sort(empty, NULL));
    assert(!list_has_cycle(empty));
    list_graph_destroy(empty);

    /* Single vertex */
    ListGraph *single = list_graph_create(1, false);
    assert(single != NULL);
    assert(list_graph_vertex_count(single) == 1);

    int    order[1];
    double dist[1];
    assert(list_bfs(single, 0, order, NULL) == 1);
    assert(order[0] == 0);
    assert(list_dfs(single, 0, order, NULL) == 1);
    assert(list_dijkstra(single, 0, dist, NULL));
    assert(fabs(dist[0]) < 1e-9);

    ListEdge mst_edges[1];
    size_t   ec = 0;
    /* Single vertex: MST weight is 0 */
    assert(fabs(list_mst_kruskal(single, mst_edges, &ec)) < 1e-9);
    list_graph_destroy(single);

    /* Self-loop on a directed graph should be detected as a cycle */
    ListGraph *selfloop = list_graph_create(2, true);
    assert(selfloop != NULL);
    assert(list_graph_add_edge(selfloop, 0, 0, 1.0));
    assert(list_has_cycle(selfloop));
    list_graph_destroy(selfloop);

    /* NULL safety */
    list_graph_destroy(NULL);
    assert(list_graph_vertex_count(NULL) == 0);
    assert(list_graph_edge_count(NULL) == 0);
    assert(!list_graph_is_directed(NULL));
    assert(!list_graph_has_edge(NULL, 0, 0, NULL));
    assert(!list_graph_add_edge(NULL, 0, 1, 1.0));
    assert(!list_graph_remove_edge(NULL, 0, 1));
}

/* ------------------------------------------------------------------ */
/*  7. Dynamic vertices: add, lazy remove, ID reuse                     */
/* ------------------------------------------------------------------ */

static void test_dynamic_vertices(void)
{
    ListGraph *graph = list_graph_create(0, false);
    assert(graph != NULL);
    assert(list_graph_vertex_count(graph) == 0);
    assert(list_graph_id_bound(graph) == 0);

    int a = list_graph_add_vertex(graph);
    int b = list_graph_add_vertex(graph);
    int c = list_graph_add_vertex(graph);
    assert(a == 0 && b == 1 && c == 2);
    assert(list_graph_add_edge(graph, a, b, 1.0));
    assert(list_graph_add_edge(graph, b, c, 2.0));
    assert(list_graph_add_edge(graph, c, c, 3.0)); /* self-loop */
    assert(list_graph_edge_count(graph) == 3);

    uint64_t before = list_graph_version(graph);
    assert(list_graph_remove_vertex(graph, b));
    assert(list_graph_version(graph) > before);
    assert(!list_graph_remove_vertex(graph, b));  /* already dead */
    assert(!list_graph_is_alive(graph, b));
    assert(list_graph_vertex_count(graph) == 2);
    assert(list_graph_id_bound(graph) == 3);
    assert(list_graph_edge_count(graph) == 1);    /* only the self-loop left */
    assert(list_graph_degree(graph, a) == 0);
    assert(!list_graph_add_edge(graph, a, b, 1.0)); /* dead endpoint */
    assert(!list_graph_has_edge(graph, a, b, NULL));

    /* Algorithms ignore the dead slot. */
    int order[3], parent[3], component[3];
    assert(list_bfs(graph, a, order, parent) == 1);
    assert(!list_bfs(graph, b, order, parent));
    assert(list_connected_components(graph, component) == 2);
    assert(component[b] == -1);
    assert(component[a] != component[c]);

    /* The dead ID is reused before the ID space grows. */
    assert(list_graph_add_vertex(graph) == b);
    assert(list_graph_degree(graph, b) == 0);
    assert(list_graph_id_bound(graph) == 3);
    assert(list_graph_add_vertex(graph) == 3);
    list_graph_destroy(graph);

    /* Directed: removing a vertex must also delete its in-edges. */
    ListGraph *dg = list_graph_create(4, true);
    assert(dg != NULL);
    assert(list_graph_add_edge(dg, 0, 3, -1.0));
    assert(list_graph_add_edge(dg, 1, 3, 1.0));
    assert(list_graph_add_edge(dg, 3, 2, 1.0));
    assert(list_graph_add_edge(dg, 3, 3, 1.0));
    assert(list_graph_add_edge(dg, 0, 1, 1.0));
    double dist[4];
    assert(!list_dijkstra(dg, 0, dist, NULL));    /* negative edge present */
    assert(list_graph_remove_vertex(dg, 3));
    assert(list_graph_edge_count(dg) == 1);
    assert(list_graph_degree(dg, 0) == 1);
    assert(list_graph_degree(dg, 1) == 0);
    assert(list_dijkstra(dg, 0, dist, NULL));     /* negative edge gone with it */
    assert(fabs(dist[1] - 1.0) < 1e-9);
    list_graph_destroy(dg);
}

/* ------------------------------------------------------------------ */
/*  8. Compaction and remap                                             */
/* ------------------------------------------------------------------ */

static void test_compaction(void)
{
    ListGraph *graph = list_graph_create(8, false);
    assert(graph != NULL);
    /* Path 0-1-2-3-4-5-6-7 with weight = left endpoint + 1 */
    for (int i = 0; i < 7; ++i) {
        assert(list_graph_add_edge(graph, i, i + 1, i + 1.0));
    }
    assert(!list_graph_needs_compaction(graph));

    assert(list_graph_remove_vertex(graph, 1));
    assert(list_graph_remove_vertex(graph, 4));
    assert(list_graph_remove_vertex(graph, 5));
    assert(list_graph_needs_compaction(graph));   /* 3 of 8 dead > 25% */
    assert(list_graph_edge_count(graph) == 2);    /* 2-3 and 6-7 survive */

    int remap[8];
    assert(list_graph_compact(graph, remap));
    const int expect[8] = { 0, -1, 1, 2, -1, -1, 3, 4 };
    for (int i = 0; i < 8; ++i) {
        assert(remap[i] == expect[i]);
    }
    assert(list_graph_id_bound(graph) == 5);
    assert(list_graph_vertex_count(graph) == 5);
    assert(!list_graph_needs_compaction(graph));
    assert(list_graph_edge_count(graph) == 2);

    double weight = 0.0;
    assert(list_graph_has_edge(graph, 1, 2, &weight)); /* was 2-3, weight 3 */
    assert(fabs(weight - 3.0) < 1e-9);
    assert(list_graph_has_edge(graph, 4, 3, &weight)); /* was 7-6, weight 7 */
    assert(fabs(weight - 7.0) < 1e-9);
    assert(!list_graph_has_edge(graph, 0, 1, NULL));

    int component[5];
    assert(list_connected_components(graph, component) == 3);

    /* New vertices extend the compacted ID space. */
    assert(list_graph_add_vertex(graph) == 5);

    /* Compaction with nothing dead is a no-op that returns the identity. */
    int ident[6];
    assert(list_graph_compact(graph, ident));
    for (int i = 0; i < 6; ++i) {
        assert(ident[i] == i);
    }
    assert(list_graph_compact(graph, NULL));      /* remap is optional */
    list_graph_destroy(graph);

    /* Remove everything, then compact down to an empty graph. */
    ListGraph *all = list_graph_create(3, true);
    assert(all != NULL);
    for (int i = 0; i < 3; ++i) {
        assert(list_graph_remove_vertex(all, i));
    }
    assert(list_graph_compact(all, NULL));
    assert(list_graph_id_bound(all) == 0);
    assert(list_graph_add_vertex(all) == 0);
    list_graph_destroy(all);
}

/* ------------------------------------------------------------------ */
/*  9. Growth across segment boundaries                                 */
/* ------------------------------------------------------------------ */

static void test_segment_growth(void)
{
    /* Larger than any plausible segment, so the directory must grow. */
    const int start = 3000, extra = 5000, total = start + extra;
    ListGraph *graph = list_graph_create((size_t)start, false);
    assert(graph != NULL);
    for (int i = 0; i < extra; ++i) {
        assert(list_graph_add_vertex(graph) == start + i);
    }
    assert(list_graph_vertex_count(graph) == (size_t)total);

    for (int i = 0; i + 1 < total; ++i) {
        assert(list_graph_add_edge(graph, i, i + 1, 1.0));
    }
    /* A shortcut that jumps across several segments. */
    assert(list_graph_add_edge(graph, 5, total - 5, 1.0));

    int *order = malloc((size_t)total * sizeof(*order));
    int *parent = malloc((size_t)total * sizeof(*parent));
    assert(order != NULL && parent != NULL);
    assert(list_bfs(graph, 0, order, parent) == (size_t)total);
    assert(parent[total - 5] == 5);

    /* Free a whole segment's worth and compact: segments must be released. */
    size_t bytes_before = list_graph_memory_bytes(graph);
    for (int i = 1000; i < 5000; ++i) {
        assert(list_graph_remove_vertex(graph, i));
    }
    assert(list_graph_needs_compaction(graph));
    assert(list_graph_compact(graph, NULL));
    assert(list_graph_vertex_count(graph) == (size_t)(total - 4000));
    assert(list_graph_id_bound(graph) == (size_t)(total - 4000));
    assert(list_graph_memory_bytes(graph) < bytes_before);

    free(order);
    free(parent);
    list_graph_destroy(graph);
}

/* ------------------------------------------------------------------ */
/*  10. Large graphs: no recursion limit, no buffer overflow            */
/* ------------------------------------------------------------------ */

static void test_large_graphs(void)
{
    /* A 200k-vertex path: the old recursive cycle check overflowed the stack. */
    const int n = 200000;
    ListGraph *path = list_graph_create((size_t)n, false);
    assert(path != NULL);
    for (int i = 0; i + 1 < n; ++i) {
        assert(list_graph_add_edge(path, i, i + 1, 1.0));
    }
    assert(!list_has_cycle(path));

    int    *order = malloc((size_t)n * sizeof(*order));
    int    *parent = malloc((size_t)n * sizeof(*parent));
    double *dist  = malloc((size_t)n * sizeof(*dist));
    assert(order != NULL && parent != NULL && dist != NULL);

    assert(list_dfs(path, 0, order, parent) == (size_t)n);
    assert(order[n - 1] == n - 1);
    assert(list_dijkstra(path, 0, dist, parent));
    assert(fabs(dist[n - 1] - (n - 1)) < 1e-6);
    assert(parent[n - 1] == n - 2);

    assert(list_graph_add_edge(path, n - 1, 0, 1.0));
    assert(list_has_cycle(path));
    list_graph_destroy(path);

    /* Directed chain: weak components via union-find, in O(V + E). */
    ListGraph *chain = list_graph_create((size_t)n, true);
    assert(chain != NULL);
    for (int i = 0; i + 1 < n; ++i) {
        assert(list_graph_add_edge(chain, i, i + 1, 1.0));
    }
    assert(list_connected_components(chain, parent) == 1);
    assert(list_topological_sort(chain, order));
    assert(order[0] == 0 && order[n - 1] == n - 1);
    list_graph_destroy(chain);

    /* Complete directed graph: the old DFS pushed up to E entries onto an
     * n-entry stack. */
    const int k = 300;
    ListGraph *complete = list_graph_create((size_t)k, true);
    assert(complete != NULL);
    for (int u = 0; u < k; ++u) {
        for (int v = 0; v < k; ++v) {
            if (u != v) {
                assert(list_graph_add_edge(complete, u, v, 1.0));
            }
        }
    }
    assert(list_dfs(complete, 0, order, parent) == (size_t)k);
    assert(list_bfs(complete, 0, order, parent) == (size_t)k);
    assert(list_has_cycle(complete));
    list_graph_destroy(complete);

    free(order);
    free(parent);
    free(dist);
}

/* ------------------------------------------------------------------ */
/*  11. Randomised comparison against a reference model                 */
/* ------------------------------------------------------------------ */

#define MODEL_N 96

typedef struct {
    bool   alive[MODEL_N];
    bool   edge[MODEL_N][MODEL_N];
    double w[MODEL_N][MODEL_N];
    size_t bound;
    bool   directed;
} Model;

static unsigned long long rng_state = 88172645463325252ULL;

static unsigned rng(unsigned bound)
{
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return (unsigned)(rng_state % bound);
}

static size_t model_vertices(const Model *m)
{
    size_t c = 0;
    for (size_t i = 0; i < m->bound; ++i) {
        c += m->alive[i];
    }
    return c;
}

static size_t model_edges(const Model *m)
{
    size_t c = 0;
    for (size_t u = 0; u < m->bound; ++u) {
        for (size_t v = 0; v < m->bound; ++v) {
            if (m->edge[u][v] && (m->directed || u <= v)) {
                ++c;
            }
        }
    }
    return c;
}

/* Weak component count of the model (edges treated as undirected). */
static size_t model_components(const Model *m, int *comp)
{
    size_t count = 0;
    for (size_t i = 0; i < m->bound; ++i) {
        comp[i] = -1;
    }
    for (size_t s = 0; s < m->bound; ++s) {
        if (!m->alive[s] || comp[s] != -1) {
            continue;
        }
        int queue[MODEL_N];
        size_t head = 0, tail = 0;
        queue[tail++] = (int)s;
        comp[s] = (int)count;
        while (head < tail) {
            const int u = queue[head++];
            for (size_t v = 0; v < m->bound; ++v) {
                if ((m->edge[u][v] || m->edge[v][u]) && comp[v] == -1) {
                    comp[v] = (int)count;
                    queue[tail++] = (int)v;
                }
            }
        }
        ++count;
    }
    return count;
}

/* Directed reachability count from s (undirected models store both arcs). */
static size_t model_reach(const Model *m, int s)
{
    bool seen[MODEL_N] = { false };
    int  queue[MODEL_N];
    size_t head = 0, tail = 0;
    queue[tail++] = s;
    seen[s] = true;
    while (head < tail) {
        const int u = queue[head++];
        for (size_t v = 0; v < m->bound; ++v) {
            if (m->edge[u][v] && !seen[v]) {
                seen[v] = true;
                queue[tail++] = (int)v;
            }
        }
    }
    return tail;
}

static void check_against_model(const ListGraph *g, const Model *m)
{
    assert(list_graph_id_bound(g) == m->bound);
    assert(list_graph_vertex_count(g) == model_vertices(m));
    assert(list_graph_edge_count(g) == model_edges(m));

    for (size_t u = 0; u < m->bound; ++u) {
        assert(list_graph_is_alive(g, (int)u) == m->alive[u]);
        size_t degree = 0;
        for (size_t v = 0; v < m->bound; ++v) {
            double w = 0.0;
            const bool has = list_graph_has_edge(g, (int)u, (int)v, &w);
            assert(has == m->edge[u][v]);
            if (has) {
                assert(w == m->w[u][v]);
                ++degree;
            }
        }
        assert(list_graph_degree(g, (int)u) == degree);
    }

    int component[MODEL_N], expect[MODEL_N];
    assert(list_connected_components(g, component) == model_components(m, expect));
    for (size_t i = 0; i < m->bound; ++i) {
        assert(component[i] == expect[i]);
    }

    int order[MODEL_N], parent[MODEL_N];
    for (size_t s = 0; s < m->bound; ++s) {
        if (m->alive[s]) {
            assert(list_bfs(g, (int)s, order, parent) == model_reach(m, (int)s));
            assert(list_dfs(g, (int)s, order, parent) == model_reach(m, (int)s));
            break; /* one start vertex per check keeps the test fast */
        }
    }
}

static void model_compact(Model *m, const int *remap)
{
    Model next = *m;
    memset(next.alive, 0, sizeof(next.alive));
    memset(next.edge, 0, sizeof(next.edge));
    size_t live = 0;
    for (size_t u = 0; u < m->bound; ++u) {
        if (remap[u] < 0) {
            continue;
        }
        ++live;
        next.alive[remap[u]] = true;
        for (size_t v = 0; v < m->bound; ++v) {
            if (m->edge[u][v]) {
                next.edge[remap[u]][remap[v]] = true;
                next.w[remap[u]][remap[v]]    = m->w[u][v];
            }
        }
    }
    next.bound = live;
    *m = next;
}

static void run_random_model(bool directed, int steps)
{
    ListGraph *g = list_graph_create(8, directed);
    assert(g != NULL);
    Model m;
    memset(&m, 0, sizeof(m));
    m.directed = directed;
    m.bound    = 8;
    for (int i = 0; i < 8; ++i) {
        m.alive[i] = true;
    }

    for (int step = 0; step < steps; ++step) {
        const unsigned op = rng(100);
        const int u = (int)rng((unsigned)m.bound + 1); /* sometimes out of range */
        const int v = (int)rng((unsigned)m.bound + 1);
        const bool in_range = (size_t)u < m.bound && (size_t)v < m.bound;

        if (op < 12) { /* add vertex */
            size_t dead = m.bound - model_vertices(&m);
            if (dead == 0 && m.bound >= MODEL_N) {
                continue;
            }
            const int id = list_graph_add_vertex(g);
            assert(id >= 0 && (size_t)id < MODEL_N);
            assert((size_t)id >= m.bound || !m.alive[id]); /* new or recycled */
            if ((size_t)id == m.bound) {
                ++m.bound;
            }
            m.alive[id] = true;
        } else if (op < 22) { /* remove vertex */
            const bool ok = list_graph_remove_vertex(g, u);
            const bool should = (size_t)u < m.bound && m.alive[u];
            assert(ok == should);
            if (ok) {
                m.alive[u] = false;
                for (size_t i = 0; i < m.bound; ++i) {
                    m.edge[u][i] = m.edge[i][u] = false;
                }
            }
        } else if (op < 70) { /* add edge */
            const double w = (double)rng(50) / 4.0;
            const bool ok = list_graph_add_edge(g, u, v, w);
            const bool should = in_range && m.alive[u] && m.alive[v] && !m.edge[u][v];
            assert(ok == should);
            if (ok) {
                m.edge[u][v] = true;
                m.w[u][v]    = w;
                if (!directed) {
                    m.edge[v][u] = true;
                    m.w[v][u]    = w;
                }
            }
        } else if (op < 92) { /* remove edge */
            const bool ok = list_graph_remove_edge(g, u, v);
            const bool should = in_range && m.alive[u] && m.alive[v] && m.edge[u][v];
            assert(ok == should);
            if (ok) {
                m.edge[u][v] = false;
                if (!directed) {
                    m.edge[v][u] = false;
                }
            }
        } else if (list_graph_needs_compaction(g) || op >= 98) { /* compact */
            int remap[MODEL_N];
            assert(list_graph_compact(g, remap));
            model_compact(&m, remap);
        }

        if (step % 7 == 0) {
            check_against_model(g, &m);
        }
    }
    check_against_model(g, &m);
    list_graph_destroy(g);
}

static void test_random_against_model(void)
{
    run_random_model(false, 4000);
    run_random_model(true, 4000);
}

/* ------------------------------------------------------------------ */
/*  main                                                                */
/* ------------------------------------------------------------------ */

int main(void)
{
    test_graph_operations();
    test_traversals_and_components();
    test_cycles_and_topological_sort();
    test_shortest_paths();
    test_minimum_spanning_tree();
    test_edge_cases();
    test_dynamic_vertices();
    test_compaction();
    test_segment_growth();
    test_large_graphs();
    test_random_against_model();
    puts("All adjacency-list tests passed.");
    return 0;
}
