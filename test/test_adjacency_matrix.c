#include "adjacency_matrix.h"

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void test_graph_operations(void) {
    MatrixGraph *graph = matrix_graph_create(3, false);
    assert(graph != NULL);
    assert(matrix_graph_vertex_count(graph) == 3);
    assert(matrix_graph_edge_count(graph) == 0);
    assert(!matrix_graph_is_directed(graph));

    assert(matrix_graph_add_edge(graph, 0, 1, 2.5));
    assert(!matrix_graph_add_edge(graph, 0, 1, 9.0));
    assert(!matrix_graph_add_edge(graph, -1, 1, 1.0));
    assert(!matrix_graph_add_edge(graph, 0, 3, 1.0));
    assert(!matrix_graph_add_edge(graph, 0, 2, INFINITY));
    assert(!matrix_graph_add_edge(graph, 0, 2, NAN));
    assert(matrix_graph_edge_count(graph) == 1);

    double weight = 0.0;
    assert(matrix_graph_has_edge(graph, 1, 0, &weight));
    assert(fabs(weight - 2.5) < 1e-9);
    assert(matrix_graph_remove_edge(graph, 1, 0));
    assert(!matrix_graph_has_edge(graph, 0, 1, NULL));
    assert(!matrix_graph_remove_edge(graph, 0, 1));
    assert(matrix_graph_edge_count(graph) == 0);
    matrix_graph_destroy(graph);
}

static void test_traversals_and_components(void) {
    MatrixGraph *graph = matrix_graph_create(6, false);
    assert(graph != NULL);
    assert(matrix_graph_add_edge(graph, 0, 1, 1.0));
    assert(matrix_graph_add_edge(graph, 0, 2, 1.0));
    assert(matrix_graph_add_edge(graph, 1, 3, 1.0));
    assert(matrix_graph_add_edge(graph, 2, 3, 1.0));
    assert(matrix_graph_add_edge(graph, 4, 5, 1.0));

    int order[6];
    int parent[6];
    assert(matrix_bfs(graph, 0, order, parent) == 4);
    const int expected_bfs[] = {0, 1, 2, 3};
    for (size_t i = 0; i < 4; ++i) {
        assert(order[i] == expected_bfs[i]);
    }
    assert(parent[0] == -1 && parent[1] == 0 && parent[2] == 0 && parent[3] == 1);

    assert(matrix_dfs(graph, 0, order, parent) == 4);
    const int expected_dfs[] = {0, 1, 3, 2};
    for (size_t i = 0; i < 4; ++i) {
        assert(order[i] == expected_dfs[i]);
    }

    int component[6];
    assert(matrix_connected_components(graph, component) == 2);
    assert(component[0] == component[3]);
    assert(component[4] == component[5]);
    assert(component[0] != component[4]);
    assert(matrix_has_cycle(graph));
    matrix_graph_destroy(graph);

    MatrixGraph *directed = matrix_graph_create(4, true);
    assert(directed != NULL);
    assert(matrix_graph_add_edge(directed, 0, 1, 1.0));
    assert(matrix_graph_add_edge(directed, 2, 1, 1.0));
    assert(matrix_connected_components(directed, component) == 2);
    matrix_graph_destroy(directed);
}

static void test_cycles_and_topological_sort(void) {
    MatrixGraph *graph = matrix_graph_create(6, true);
    assert(graph != NULL);
    assert(matrix_graph_add_edge(graph, 5, 2, 1.0));
    assert(matrix_graph_add_edge(graph, 5, 0, 1.0));
    assert(matrix_graph_add_edge(graph, 4, 0, 1.0));
    assert(matrix_graph_add_edge(graph, 4, 1, 1.0));
    assert(matrix_graph_add_edge(graph, 2, 3, 1.0));
    assert(matrix_graph_add_edge(graph, 3, 1, 1.0));
    assert(!matrix_has_cycle(graph));

    int order[6];
    int position[6];
    assert(matrix_topological_sort(graph, order));
    for (int i = 0; i < 6; ++i) {
        position[order[i]] = i;
    }
    assert(position[5] < position[2]);
    assert(position[5] < position[0]);
    assert(position[4] < position[0]);
    assert(position[4] < position[1]);
    assert(position[2] < position[3]);
    assert(position[3] < position[1]);

    assert(matrix_graph_add_edge(graph, 1, 5, 1.0));
    assert(matrix_has_cycle(graph));
    assert(!matrix_topological_sort(graph, order));
    matrix_graph_destroy(graph);

    MatrixGraph *tree = matrix_graph_create(3, false);
    assert(tree != NULL);
    assert(matrix_graph_add_edge(tree, 0, 1, 1.0));
    assert(matrix_graph_add_edge(tree, 1, 2, 1.0));
    assert(!matrix_has_cycle(tree));
    assert(!matrix_topological_sort(tree, order));
    assert(matrix_graph_add_edge(tree, 2, 0, 1.0));
    assert(matrix_has_cycle(tree));
    matrix_graph_destroy(tree);
}

static void test_shortest_paths(void) {
    MatrixGraph *graph = matrix_graph_create(5, true);
    assert(graph != NULL);
    assert(matrix_graph_add_edge(graph, 0, 1, 10.0));
    assert(matrix_graph_add_edge(graph, 0, 2, 3.0));
    assert(matrix_graph_add_edge(graph, 2, 1, 1.0));
    assert(matrix_graph_add_edge(graph, 1, 3, 2.0));
    assert(matrix_graph_add_edge(graph, 2, 3, 8.0));

    double distance[5];
    int parent[5];
    assert(matrix_dijkstra(graph, 0, distance, parent));
    assert(fabs(distance[0]) < 1e-9);
    assert(fabs(distance[1] - 4.0) < 1e-9);
    assert(fabs(distance[2] - 3.0) < 1e-9);
    assert(fabs(distance[3] - 6.0) < 1e-9);
    assert(isinf(distance[4]));
    assert(parent[1] == 2 && parent[2] == 0 && parent[3] == 1 && parent[4] == -1);
    matrix_graph_destroy(graph);

    MatrixGraph *negative = matrix_graph_create(2, true);
    assert(negative != NULL);
    assert(matrix_graph_add_edge(negative, 0, 1, -1.0));
    assert(!matrix_dijkstra(negative, 0, distance, parent));
    matrix_graph_destroy(negative);
}

static void test_minimum_spanning_tree(void) {
    MatrixGraph *graph = matrix_graph_create(4, false);
    assert(graph != NULL);
    assert(matrix_graph_add_edge(graph, 0, 1, 10.0));
    assert(matrix_graph_add_edge(graph, 0, 2, 6.0));
    assert(matrix_graph_add_edge(graph, 0, 3, 5.0));
    assert(matrix_graph_add_edge(graph, 1, 3, 15.0));
    assert(matrix_graph_add_edge(graph, 2, 3, 4.0));

    MatrixEdge edges[3];
    size_t edge_count = 0;
    const double total = matrix_mst_kruskal(graph, edges, &edge_count);
    assert(edge_count == 3);
    assert(fabs(total - 19.0) < 1e-9);
    matrix_graph_destroy(graph);

    MatrixGraph *disconnected = matrix_graph_create(3, false);
    assert(disconnected != NULL);
    assert(matrix_graph_add_edge(disconnected, 0, 1, 1.0));
    assert(isinf(matrix_mst_kruskal(disconnected, edges, &edge_count)));
    assert(edge_count == 1);
    matrix_graph_destroy(disconnected);
}

static void test_empty_graph(void) {
    MatrixGraph *graph = matrix_graph_create(0, true);
    assert(graph != NULL);
    assert(matrix_connected_components(graph, NULL) == 0);
    assert(matrix_topological_sort(graph, NULL));
    assert(!matrix_has_cycle(graph));
    matrix_graph_destroy(graph);
}

/* A real DFS never has a cross edge in an undirected graph, so on a triangle
 * the last vertex must hang off the previous one, not off the root. */
static void test_dfs_builds_a_real_tree(void) {
    MatrixGraph *graph = matrix_graph_create(3, false);
    assert(graph != NULL);
    assert(matrix_graph_add_edge(graph, 0, 1, 1.0));
    assert(matrix_graph_add_edge(graph, 0, 2, 1.0));
    assert(matrix_graph_add_edge(graph, 1, 2, 1.0));

    int order[3];
    int parent[3];
    assert(matrix_dfs(graph, 0, order, parent) == 3);
    assert(order[0] == 0 && order[1] == 1 && order[2] == 2);
    assert(parent[0] == -1 && parent[1] == 0 && parent[2] == 1);
    matrix_graph_destroy(graph);
}

static void test_dynamic_vertices(void) {
    MatrixGraph *graph = matrix_graph_create(0, false);
    assert(graph != NULL);
    assert(matrix_graph_vertex_count(graph) == 0);
    assert(matrix_graph_id_bound(graph) == 0);

    const int a = matrix_graph_add_vertex(graph);
    const int b = matrix_graph_add_vertex(graph);
    const int c = matrix_graph_add_vertex(graph);
    assert(a == 0 && b == 1 && c == 2);
    assert(matrix_graph_add_edge(graph, a, b, 1.0));
    assert(matrix_graph_add_edge(graph, b, c, 2.0));
    assert(matrix_graph_add_edge(graph, c, c, 3.0));   /* self-loop */
    assert(matrix_graph_edge_count(graph) == 3);

    const uint64_t before = matrix_graph_version(graph);
    assert(matrix_graph_remove_vertex(graph, b));
    assert(matrix_graph_version(graph) > before);
    assert(!matrix_graph_remove_vertex(graph, b));
    assert(!matrix_graph_is_alive(graph, b));
    assert(matrix_graph_vertex_count(graph) == 2);
    assert(matrix_graph_id_bound(graph) == 3);
    assert(matrix_graph_edge_count(graph) == 1);
    assert(matrix_graph_degree(graph, a) == 0);
    assert(!matrix_graph_add_edge(graph, a, b, 1.0));
    assert(!matrix_graph_has_edge(graph, a, b, NULL));

    int order[3];
    int parent[3];
    int component[3];
    assert(matrix_bfs(graph, a, order, parent) == 1);
    assert(matrix_bfs(graph, b, order, parent) == 0);
    assert(matrix_connected_components(graph, component) == 2);
    assert(component[b] == -1 && component[a] != component[c]);

    /* The dead ID is reused before the ID space grows. */
    assert(matrix_graph_add_vertex(graph) == b);
    assert(matrix_graph_degree(graph, b) == 0);
    assert(matrix_graph_id_bound(graph) == 3);
    assert(matrix_graph_add_vertex(graph) == 3);
    matrix_graph_destroy(graph);

    /* Directed: removing a vertex must also delete its in-edges. */
    MatrixGraph *directed = matrix_graph_create(4, true);
    assert(directed != NULL);
    assert(matrix_graph_add_edge(directed, 0, 3, -1.0));
    assert(matrix_graph_add_edge(directed, 1, 3, 1.0));
    assert(matrix_graph_add_edge(directed, 3, 2, 1.0));
    assert(matrix_graph_add_edge(directed, 3, 3, 1.0));
    assert(matrix_graph_add_edge(directed, 0, 1, 1.0));
    double distance[4];
    assert(!matrix_dijkstra(directed, 0, distance, NULL));   /* negative edge */
    assert(matrix_graph_remove_vertex(directed, 3));
    assert(matrix_graph_edge_count(directed) == 1);
    assert(matrix_graph_degree(directed, 0) == 1);
    assert(matrix_graph_degree(directed, 1) == 0);
    assert(matrix_dijkstra(directed, 0, distance, NULL));    /* gone with it */
    assert(fabs(distance[1] - 1.0) < 1e-9);
    matrix_graph_destroy(directed);
}

static void test_compaction(void) {
    MatrixGraph *graph = matrix_graph_create(8, false);
    assert(graph != NULL);
    for (int i = 0; i < 7; ++i) {
        assert(matrix_graph_add_edge(graph, i, i + 1, i + 1.0));
    }
    assert(!matrix_graph_needs_compaction(graph));

    assert(matrix_graph_remove_vertex(graph, 1));
    assert(matrix_graph_remove_vertex(graph, 4));
    assert(matrix_graph_remove_vertex(graph, 5));
    assert(matrix_graph_needs_compaction(graph));       /* 3 of 8 dead > 25% */
    assert(matrix_graph_edge_count(graph) == 2);        /* 2-3 and 6-7 survive */

    int remap[8];
    assert(matrix_graph_compact(graph, remap));
    const int expect[8] = {0, -1, 1, 2, -1, -1, 3, 4};
    for (int i = 0; i < 8; ++i) {
        assert(remap[i] == expect[i]);
    }
    assert(matrix_graph_id_bound(graph) == 5);
    assert(matrix_graph_vertex_count(graph) == 5);
    assert(!matrix_graph_needs_compaction(graph));
    assert(matrix_graph_edge_count(graph) == 2);

    double weight = 0.0;
    assert(matrix_graph_has_edge(graph, 1, 2, &weight));   /* was 2-3, weight 3 */
    assert(fabs(weight - 3.0) < 1e-9);
    assert(matrix_graph_has_edge(graph, 4, 3, &weight));   /* was 7-6, weight 7 */
    assert(fabs(weight - 7.0) < 1e-9);
    assert(!matrix_graph_has_edge(graph, 0, 1, NULL));
    assert(matrix_graph_degree(graph, 1) == 1 && matrix_graph_degree(graph, 0) == 0);

    int component[5];
    assert(matrix_connected_components(graph, component) == 3);
    assert(matrix_graph_add_vertex(graph) == 5);

    int identity[6];
    assert(matrix_graph_compact(graph, identity));          /* nothing dead */
    for (int i = 0; i < 6; ++i) {
        assert(identity[i] == i);
    }
    assert(matrix_graph_compact(graph, NULL));              /* remap optional */
    matrix_graph_destroy(graph);

    MatrixGraph *all = matrix_graph_create(3, true);
    assert(all != NULL);
    for (int i = 0; i < 3; ++i) {
        assert(matrix_graph_remove_vertex(all, i));
    }
    assert(matrix_graph_compact(all, NULL));
    assert(matrix_graph_id_bound(all) == 0);
    assert(matrix_graph_add_vertex(all) == 0);
    matrix_graph_destroy(all);
}

/* Growing past 64, 128 and 256 IDs reallocates the tile directory; edges that
 * were already stored must survive every step. */
static void test_growth_across_tiles(void) {
    MatrixGraph *graph = matrix_graph_create(60, false);
    assert(graph != NULL);
    assert(matrix_graph_add_edge(graph, 0, 59, 1.5));
    assert(matrix_graph_add_edge(graph, 10, 10, 2.5));

    const int total = 300;
    for (int i = 60; i < total; ++i) {
        assert(matrix_graph_add_vertex(graph) == i);
    }
    assert(matrix_graph_vertex_count(graph) == (size_t)total);

    double weight = 0.0;
    assert(matrix_graph_has_edge(graph, 59, 0, &weight) && weight == 1.5);
    assert(matrix_graph_has_edge(graph, 10, 10, &weight) && weight == 2.5);

    /* A chain that crosses every tile boundary, plus a long shortcut. */
    for (int i = 59; i + 1 < total; ++i) {
        assert(matrix_graph_add_edge(graph, i, i + 1, 1.0));
    }
    assert(matrix_graph_add_edge(graph, 0, total - 1, 1.0));

    int *order = malloc((size_t)total * sizeof(*order));
    int *parent = malloc((size_t)total * sizeof(*parent));
    assert(order != NULL && parent != NULL);
    /* Vertex 0 reaches itself and 59 .. total-1; vertices 1 .. 58 are isolated. */
    assert(matrix_bfs(graph, 0, order, parent) == (size_t)(1 + total - 59));
    assert(parent[59] == 0 && parent[total - 1] == 0);   /* the shortcut wins */
    assert(parent[100] == 99);
    free(order);
    free(parent);
    matrix_graph_destroy(graph);
}

/* Deep and wide graphs: no recursion limit, no buffer overflow, and a sparse
 * graph on a large ID space pays only for the tiles it touches. */
static void test_large_graphs(void) {
    const int n = 20000;
    MatrixGraph *path = matrix_graph_create((size_t)n, false);
    assert(path != NULL);
    for (int i = 0; i + 1 < n; ++i) {
        assert(matrix_graph_add_edge(path, i, i + 1, 1.0));
    }

    /* A dense n x n matrix would need n * n * 9 bytes (3.6 GB here).  A path
     * pays n^2 / 8 bytes for the presence bits (50 MB) plus only the weight
     * tiles along the diagonal band (about 30 MB). */
    assert(matrix_graph_memory_bytes(path) < (size_t)128 * 1024 * 1024);
    assert(matrix_graph_memory_bytes(path) < (size_t)n * (size_t)n * 9 / 20);
    assert(!matrix_has_cycle(path));    /* the recursive version overflowed here */

    int *order = malloc((size_t)n * sizeof(*order));
    int *parent = malloc((size_t)n * sizeof(*parent));
    double *distance = malloc((size_t)n * sizeof(*distance));
    assert(order != NULL && parent != NULL && distance != NULL);

    assert(matrix_dfs(path, 0, order, parent) == (size_t)n);
    assert(order[n - 1] == n - 1 && parent[n - 1] == n - 2);
    assert(matrix_bfs(path, 0, order, parent) == (size_t)n);
    assert(matrix_dijkstra(path, 0, distance, parent));
    assert(fabs(distance[n - 1] - (n - 1)) < 1e-6);
    assert(matrix_connected_components(path, parent) == 1);

    assert(matrix_graph_add_edge(path, n - 1, 0, 1.0));
    assert(matrix_has_cycle(path));
    matrix_graph_destroy(path);

    MatrixGraph *chain = matrix_graph_create((size_t)n, true);
    assert(chain != NULL);
    for (int i = 0; i + 1 < n; ++i) {
        assert(matrix_graph_add_edge(chain, i, i + 1, 1.0));
    }
    assert(matrix_connected_components(chain, parent) == 1);
    assert(matrix_topological_sort(chain, order));
    assert(order[0] == 0 && order[n - 1] == n - 1);
    matrix_graph_destroy(chain);

    const int k = 300;
    MatrixGraph *complete = matrix_graph_create((size_t)k, true);
    assert(complete != NULL);
    for (int u = 0; u < k; ++u) {
        for (int v = 0; v < k; ++v) {
            if (u != v) {
                assert(matrix_graph_add_edge(complete, u, v, 1.0));
            }
        }
    }
    assert(matrix_graph_edge_count(complete) == (size_t)k * (k - 1));
    assert(matrix_dfs(complete, 0, order, parent) == (size_t)k);
    assert(matrix_bfs(complete, 0, order, parent) == (size_t)k);
    assert(matrix_has_cycle(complete));
    matrix_graph_destroy(complete);

    free(order);
    free(parent);
    free(distance);
}

/* Neighbour stepping visits exactly the out-neighbours, in ascending order,
 * across tile boundaries. */
static void test_neighbor_iteration(void) {
    MatrixGraph *graph = matrix_graph_create(200, true);
    assert(graph != NULL);
    const int targets[] = {3, 63, 64, 65, 127, 128, 199};
    for (size_t i = 0; i < sizeof(targets) / sizeof(targets[0]); ++i) {
        assert(matrix_graph_add_edge(graph, 5, targets[i], targets[i] + 0.5));
    }
    size_t cursor = 0;
    int to = -1;
    double weight = 0.0;
    for (size_t i = 0; i < sizeof(targets) / sizeof(targets[0]); ++i) {
        assert(matrix_graph_next_neighbor(graph, 5, &cursor, &to, &weight));
        assert(to == targets[i] && weight == targets[i] + 0.5);
    }
    assert(!matrix_graph_next_neighbor(graph, 5, &cursor, &to, &weight));
    assert(matrix_graph_degree(graph, 5) == 7);

    cursor = 0;
    assert(!matrix_graph_next_neighbor(graph, 6, &cursor, &to, &weight));   /* no edges */
    assert(!matrix_graph_next_neighbor(graph, 500, &cursor, &to, &weight)); /* bad id */
    matrix_graph_destroy(graph);
}

/* ---- Randomised comparison against a reference model ---- */

#define MODEL_N 160   /* crosses the 64 and 128 tile boundaries */

typedef struct {
    bool alive[MODEL_N];
    bool edge[MODEL_N][MODEL_N];
    double w[MODEL_N][MODEL_N];
    size_t bound;
    bool directed;
} Model;

static unsigned long long rng_state = 88172645463325252ULL;

static unsigned rng(unsigned bound) {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return (unsigned)(rng_state % bound);
}

static size_t model_vertices(const Model *m) {
    size_t count = 0;
    for (size_t i = 0; i < m->bound; ++i) {
        count += m->alive[i];
    }
    return count;
}

static size_t model_edges(const Model *m) {
    size_t count = 0;
    for (size_t u = 0; u < m->bound; ++u) {
        for (size_t v = 0; v < m->bound; ++v) {
            if (m->edge[u][v] && (m->directed || u <= v)) {
                ++count;
            }
        }
    }
    return count;
}

/* Weak component count of the model (edges treated as undirected). */
static size_t model_components(const Model *m, int *comp) {
    size_t count = 0;
    for (size_t i = 0; i < m->bound; ++i) {
        comp[i] = -1;
    }
    for (size_t s = 0; s < m->bound; ++s) {
        if (!m->alive[s] || comp[s] != -1) {
            continue;
        }
        int queue[MODEL_N];
        size_t head = 0;
        size_t tail = 0;
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

static size_t model_reach(const Model *m, int s) {
    bool seen[MODEL_N] = {false};
    int queue[MODEL_N];
    size_t head = 0;
    size_t tail = 0;
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

static void check_against_model(const MatrixGraph *g, const Model *m) {
    assert(matrix_graph_id_bound(g) == m->bound);
    assert(matrix_graph_vertex_count(g) == model_vertices(m));
    assert(matrix_graph_edge_count(g) == model_edges(m));

    for (size_t u = 0; u < m->bound; ++u) {
        assert(matrix_graph_is_alive(g, (int)u) == m->alive[u]);
        size_t degree = 0;
        for (size_t v = 0; v < m->bound; ++v) {
            double w = 0.0;
            const bool has = matrix_graph_has_edge(g, (int)u, (int)v, &w);
            assert(has == m->edge[u][v]);
            if (has) {
                assert(w == m->w[u][v]);
                ++degree;
            }
        }
        assert(matrix_graph_degree(g, (int)u) == degree);
    }

    int component[MODEL_N];
    int expect[MODEL_N];
    assert(matrix_connected_components(g, component) == model_components(m, expect));
    for (size_t i = 0; i < m->bound; ++i) {
        assert(component[i] == expect[i]);
    }

    int order[MODEL_N];
    int parent[MODEL_N];
    for (size_t s = 0; s < m->bound; ++s) {
        if (m->alive[s]) {
            assert(matrix_bfs(g, (int)s, order, parent) == model_reach(m, (int)s));
            assert(matrix_dfs(g, (int)s, order, parent) == model_reach(m, (int)s));
            break;   /* one start vertex per check keeps the test fast */
        }
    }
}

static void model_compact(Model *m, const int *remap) {
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
                next.w[remap[u]][remap[v]] = m->w[u][v];
            }
        }
    }
    next.bound = live;
    *m = next;
}

static void run_random_model(bool directed, int steps) {
    MatrixGraph *g = matrix_graph_create(8, directed);
    assert(g != NULL);
    Model m;
    memset(&m, 0, sizeof(m));
    m.directed = directed;
    m.bound = 8;
    for (int i = 0; i < 8; ++i) {
        m.alive[i] = true;
    }

    for (int step = 0; step < steps; ++step) {
        const unsigned op = rng(100);
        const int u = (int)rng((unsigned)m.bound + 1);   /* sometimes out of range */
        const int v = (int)rng((unsigned)m.bound + 1);
        const bool in_range = (size_t)u < m.bound && (size_t)v < m.bound;

        if (op < 12) {   /* add vertex */
            const size_t dead = m.bound - model_vertices(&m);
            if (dead == 0 && m.bound >= MODEL_N) {
                continue;
            }
            const int id = matrix_graph_add_vertex(g);
            assert(id >= 0 && (size_t)id < MODEL_N);
            assert((size_t)id >= m.bound || !m.alive[id]);   /* new or recycled */
            if ((size_t)id == m.bound) {
                ++m.bound;
            }
            m.alive[id] = true;
        } else if (op < 22) {   /* remove vertex */
            const bool ok = matrix_graph_remove_vertex(g, u);
            const bool should = (size_t)u < m.bound && m.alive[u];
            assert(ok == should);
            if (ok) {
                m.alive[u] = false;
                for (size_t i = 0; i < m.bound; ++i) {
                    m.edge[u][i] = m.edge[i][u] = false;
                }
            }
        } else if (op < 70) {   /* add edge */
            const double w = (double)rng(50) / 4.0 - (rng(10) == 0 ? 3.0 : 0.0);
            const bool ok = matrix_graph_add_edge(g, u, v, w);
            const bool should = in_range && m.alive[u] && m.alive[v] && !m.edge[u][v];
            assert(ok == should);
            if (ok) {
                m.edge[u][v] = true;
                m.w[u][v] = w;
                if (!directed) {
                    m.edge[v][u] = true;
                    m.w[v][u] = w;
                }
            }
        } else if (op < 92) {   /* remove edge */
            const bool ok = matrix_graph_remove_edge(g, u, v);
            const bool should = in_range && m.alive[u] && m.alive[v] && m.edge[u][v];
            assert(ok == should);
            if (ok) {
                m.edge[u][v] = false;
                if (!directed) {
                    m.edge[v][u] = false;
                }
            }
        } else if (matrix_graph_needs_compaction(g) || op >= 98) {   /* compact */
            int remap[MODEL_N];
            assert(matrix_graph_compact(g, remap));
            model_compact(&m, remap);
        }

        if (step % 7 == 0) {
            check_against_model(g, &m);
        }
    }
    check_against_model(g, &m);

    /* Dijkstra must agree with the model about whether a negative edge exists
     * (the graph tracks that with a counter, not a scan). */
    bool negative = false;
    for (size_t u = 0; u < m.bound; ++u) {
        for (size_t v = 0; v < m.bound; ++v) {
            negative = negative || (m.edge[u][v] && m.w[u][v] < 0.0);
        }
    }
    for (size_t s = 0; s < m.bound; ++s) {
        if (m.alive[s]) {
            double distance[MODEL_N];
            assert(matrix_dijkstra(g, (int)s, distance, NULL) == !negative);
            break;
        }
    }
    matrix_graph_destroy(g);
}

static void test_random_against_model(void) {
    run_random_model(false, 4000);
    run_random_model(true, 4000);
}

int main(void) {
    test_graph_operations();
    test_traversals_and_components();
    test_cycles_and_topological_sort();
    test_shortest_paths();
    test_minimum_spanning_tree();
    test_empty_graph();
    test_dfs_builds_a_real_tree();
    test_dynamic_vertices();
    test_compaction();
    test_growth_across_tiles();
    test_large_graphs();
    test_neighbor_iteration();
    test_random_against_model();
    puts("All adjacency-matrix tests passed.");
    return 0;
}
