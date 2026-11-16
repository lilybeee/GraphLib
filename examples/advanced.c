/*
 * advanced.c — a tour of the algorithms in graph_advanced.h (graph_adv_*): the ones beyond the
 * everyday set.  Each section is a few lines and prints a result you can check by hand.
 *
 *   make examples
 */

#include <stdio.h>

#include "graph_advanced.h"
#include "graph_algorithms.h"

int main(void)
{
    /* A road map: weights are distances in km.  0-1 (4), 0-2 (2), 1-2 (1), 1-3 (5), 2-3 (8), 3-4 (3). */
    Graph *roads = graph_create_sized(5, false);
    const struct { int a, b; double km; } r[] = {{0, 1, 4}, {0, 2, 2}, {1, 2, 1}, {1, 3, 5}, {2, 3, 8}, {3, 4, 3}};
    for (size_t i = 0; i < sizeof r / sizeof *r; i++) graph_add_edge(roads, r[i].a, r[i].b, r[i].km);

    /* Widest path: the route whose NARROWEST road is as wide as possible (here, weights are
     * capacities in tonnes).  Useful for "what is the heaviest truck that can get there?". */
    double width[5];
    graph_adv_widest_path(roads, 0, width, NULL);
    printf("heaviest load from 0 to 4: %.0f t\n", width[4]);

    /* Tree traversal: visit a tree in the order you need. */
    Graph *tree = graph_create_sized(7, false);       /* 0 is the root; 1,2 its children; 3-6 theirs */
    const int t[][2] = {{0, 1}, {0, 2}, {1, 3}, {1, 4}, {2, 5}, {2, 6}};
    for (size_t i = 0; i < 6; i++) graph_add_edge(tree, t[i][0], t[i][1], 1.0);
    int order[7];
    size_t n;
    const char *names[] = {"pre-order", "post-order", "level-order"};
    const GraphAdvTreeOrder kinds[] = {GRAPH_ADV_TREE_PREORDER, GRAPH_ADV_TREE_POSTORDER, GRAPH_ADV_TREE_LEVELORDER};
    for (int k = 0; k < 3; k++) {
        graph_adv_tree_traversal(tree, 0, kinds[k], order, &n);
        printf("%-11s:", names[k]);
        for (size_t i = 0; i < n; i++) printf(" %d", order[i]);
        printf("\n");
    }

    /* Isomorphism: are two networks the same shape, just numbered differently?  Optionally the
     * weights must match too.  Colour refinement makes this fast on large graphs. */
    Graph *a = graph_create_sized(5, false), *b = graph_create_sized(5, false);
    for (int i = 0; i < 5; i++) {
        graph_add_edge(a, i, (i + 1) % 5, 1.0);                 /* a 5-cycle ...                   */
        graph_add_edge(b, (i * 2) % 5, ((i + 1) * 2) % 5, 1.0); /* ... and the same cycle relabelled */
    }
    bool same;
    graph_adv_graph_isomorphism(a, b, false, &same, NULL);
    printf("5-cycle vs relabelled 5-cycle isomorphic: %s\n", same ? "yes" : "no");

    /* Subgraph isomorphism: does the big graph contain a triangle? */
    Graph *triangle = graph_create_sized(3, false);
    graph_add_edge(triangle, 0, 1, 1.0); graph_add_edge(triangle, 1, 2, 1.0); graph_add_edge(triangle, 0, 2, 1.0);
    bool found;
    graph_adv_subgraph_isomorphism(triangle, roads, true, false, NULL, &found);
    printf("roads contain a triangle (0-1-2): %s\n", found ? "yes" : "no");

    /* All-pairs shortest paths with Johnson's algorithm, then rebuild one route from `pred`. */
    double dist[25];
    int pred[25];
    graph_adv_johnson(roads, dist, pred);
    GraphPath route;
    graph_adv_apsp_path(dist, pred, 5, 0, 4, &route);
    printf("shortest road 0 -> 4: %.0f km via", route.cost);
    for (size_t i = 0; i < route.length; i++) printf(" %d", route.vertices[i]);
    printf("\n");
    graph_path_free(&route);

    /* Dynamic connectivity: ask "are u and v connected?" while the graph changes. */
    GraphAdvDynConn *dc = graph_adv_dynconn_create(roads);
    bool connected;
    graph_adv_dynconn_connected(dc, 0, 4, &connected);
    printf("\nconnected 0-4 now: %s\n", connected ? "yes" : "no");
    graph_remove_edge(roads, 3, 4);                       /* the only road to town 4 closes */
    graph_adv_dynconn_connected(dc, 0, 4, &connected);
    printf("connected 0-4 after closing road 3-4: %s\n", connected ? "yes" : "no");
    graph_adv_dynconn_destroy(dc);

    graph_destroy(roads); graph_destroy(tree); graph_destroy(a); graph_destroy(b); graph_destroy(triangle);
    return 0;
}
