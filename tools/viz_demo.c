/*
 * viz_demo — draw a generated graph.
 *   viz_demo <kind> <n> <out.png|out.svg> [matrix]
 * kinds: random (Erdos-Renyi), scalefree (Barabasi-Albert), grid, tree, clusters
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "graph_algorithms.h"
#include "graph_viz.h"

static Graph *make(const char *kind, size_t n)
{
    if (strcmp(kind, "random") == 0)    { return graph_gen_erdos_renyi(GRAPH_LIST, n, 3.0 / (double)n, 7, false, 1, 1); }
    if (strcmp(kind, "scalefree") == 0) { return graph_gen_barabasi_albert(GRAPH_LIST, n, 2, 7); }
    if (strcmp(kind, "tree") == 0)      { return graph_gen_random_tree(GRAPH_LIST, n, 7); }
    if (strcmp(kind, "grid") == 0) {
        size_t side = 1;
        while (side * side < n) { ++side; }
        return graph_gen_grid(GRAPH_LIST, side, side);
    }
    if (strcmp(kind, "clusters") == 0) {            /* dense groups of ~50 joined by a few bridges */
        Graph *g = graph_create_sized(n, false);
        const size_t size = 50;
        unsigned long long r = 12345;
        for (size_t u = 0; u < n; ++u) {
            for (int e = 0; e < 4; ++e) {
                r = r * 6364136223846793005ULL + 1442695040888963407ULL;
                const size_t base = u / size * size, span = base + size <= n ? size : n - base;
                const size_t v = base + (size_t)(r >> 33) % span;
                if (v != u) { (void)graph_add_edge(g, (int)u, (int)v, 1.0); }
            }
            if (u % size == 0 && u > 0) { (void)graph_add_edge(g, (int)u, (int)(u - size), 1.0); }
        }
        return g;
    }
    return NULL;
}

int main(int argc, char **argv)
{
    if (argc < 4) {
        fprintf(stderr, "usage: %s <random|scalefree|grid|tree|clusters> <n> <out.png|out.svg> [matrix]\n", argv[0]);
        return 2;
    }
    Graph *g = make(argv[1], (size_t)strtoul(argv[2], NULL, 10));
    if (g == NULL) {
        fprintf(stderr, "could not build graph: %s\n", graph_status_string(graph_last_error()));
        return 1;
    }
    VizOptions o = graph_viz_default_options();
    o.title = argv[1];
    if (argc > 4 && strcmp(argv[4], "matrix") == 0) { o.layout = VIZ_LAYOUT_MATRIX; }
    const GraphStatus st = graph_visualize(g, argv[3], &o);
    printf("%zu vertices, %zu edges -> %s (%s)\n", graph_vertex_count(g), graph_edge_count(g),
           argv[3], graph_status_string(st));
    graph_destroy(g);
    return st == GRAPH_OK ? 0 : 1;
}
