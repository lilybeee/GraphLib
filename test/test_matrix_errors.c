/* Error-code tests: every failing core operation must record the right code. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "adjacency_matrix.h"
#include "graph_status.h"
#include "matrix_algorithms.h"
#include "test_util_matrix.h"

#define EXPECT_ERR(expr, code) \
    do { graph_clear_error(); (void)(expr); CHECK(graph_last_error() == (code)); } while (0)

int main(void)
{
    MatrixGraph *g = matrix_graph_create(4, false);
    MatrixGraph *d = matrix_graph_create(4, true);
    CHECK(g != NULL && d != NULL);
    CHECK(matrix_graph_add_edge(g, 0, 1, 1.0));
    CHECK(graph_last_error() == GRAPH_OK);

    EXPECT_ERR(matrix_graph_add_edge(g, 0, 1, 1.0), GRAPH_ERR_DUPLICATE_EDGE);
    EXPECT_ERR(matrix_graph_add_edge(g, 0, 9, 1.0), GRAPH_ERR_INVALID_VERTEX);
    EXPECT_ERR(matrix_graph_add_edge(NULL, 0, 1, 1.0), GRAPH_ERR_INVALID_ARGUMENT);
    EXPECT_ERR(matrix_graph_add_edge(g, 1, 2, NAN), GRAPH_ERR_INVALID_WEIGHT);
    EXPECT_ERR(matrix_graph_add_edge(g, 1, 2, INFINITY), GRAPH_ERR_INVALID_WEIGHT);
    EXPECT_ERR(matrix_graph_remove_edge(g, 2, 3), GRAPH_ERR_EDGE_NOT_FOUND);
    EXPECT_ERR(matrix_graph_remove_edge(g, -1, 3), GRAPH_ERR_INVALID_VERTEX);
    EXPECT_ERR(matrix_graph_remove_vertex(g, 77), GRAPH_ERR_INVALID_VERTEX);
    EXPECT_ERR(matrix_graph_add_vertex(NULL), GRAPH_ERR_INVALID_ARGUMENT);

    int order[8], parent[8];
    EXPECT_ERR(matrix_bfs(g, 50, order, parent), GRAPH_ERR_INVALID_VERTEX);
    EXPECT_ERR(matrix_bfs(g, 0, NULL, parent), GRAPH_ERR_INVALID_ARGUMENT);
    EXPECT_ERR(matrix_dfs(g, 50, order, parent), GRAPH_ERR_INVALID_VERTEX);
    EXPECT_ERR(matrix_topological_sort(g, order), GRAPH_ERR_UNSUPPORTED);
    EXPECT_ERR(matrix_has_cycle(NULL), GRAPH_ERR_INVALID_ARGUMENT);

    /* Cycle in a directed graph -> GRAPH_ERR_CYCLE. */
    CHECK(matrix_graph_add_edge(d, 0, 1, 1.0));
    CHECK(matrix_graph_add_edge(d, 1, 2, 1.0));
    CHECK(matrix_graph_add_edge(d, 2, 0, 1.0));
    EXPECT_ERR(matrix_topological_sort(d, order), GRAPH_ERR_CYCLE);

    /* Negative weight -> Dijkstra refuses with a specific code. */
    CHECK(matrix_graph_add_edge(g, 2, 3, -2.0));
    double dist[8];
    EXPECT_ERR(matrix_dijkstra(g, 0, dist, parent), GRAPH_ERR_NEGATIVE_WEIGHT);

    /* Disconnected graph -> Kruskal reports NOT_FOUND. */
    MatrixGraph *h = matrix_graph_create(4, false);
    MatrixEdge edges[4];
    size_t n_edges;
    CHECK(matrix_graph_add_edge(h, 0, 1, 1.0));
    graph_clear_error();
    CHECK(isinf(matrix_mst_kruskal(h, edges, &n_edges)));
    CHECK(graph_last_error() == GRAPH_ERR_NOT_FOUND);
    graph_clear_error();
    CHECK(isinf(matrix_mst_kruskal(d, edges, &n_edges)));
    CHECK(graph_last_error() == GRAPH_ERR_UNSUPPORTED);

    /* Status-returning library: codes come back directly. */
    long long tri;
    CHECK(matrix_triangle_count(NULL, &tri) == GRAPH_ERR_INVALID_ARGUMENT);

    /* Every code has a distinct non-empty name. */
    for (int a = 0; a < GRAPH_ERR_COUNT; ++a) {
        const char *sa = graph_status_string((GraphStatus)a);
        CHECK(sa != NULL && sa[0] != '\0');
        for (int b = a + 1; b < GRAPH_ERR_COUNT; ++b) {
            CHECK(strcmp(sa, graph_status_string((GraphStatus)b)) != 0);
        }
    }
    CHECK(graph_status_string((GraphStatus)999) != NULL);

    matrix_graph_destroy(g);
    matrix_graph_destroy(d);
    matrix_graph_destroy(h);
    puts("All matrix error-code tests passed.");
    return 0;
}
