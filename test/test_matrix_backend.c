/*
 * test_matrix_backend.c — the matrix exposed through its GraphBackend vtable,
 * and the shape of its cost model.
 *
 * Compile & run:  make matrix_backend_test && ./matrix_backend_test
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "matrix_backend.h"

static void test_vtable(void) {
    const GraphBackend *b = matrix_backend();
    assert(strcmp(b->name, "matrix") == 0);

    /* Everything the observer and migration rely on is present. */
    assert(b->create && b->destroy && b->iter_begin && b->iter_next);
    assert(b->estimate_cost && b->estimate_memory && b->estimate_build_cost);
    assert(b->bfs && b->dfs && b->dijkstra && b->mst_kruskal);
    assert(b->add_new_edge == NULL);   /* the duplicate check is already O(1) */

    void *g = b->create(4, false);
    assert(g != NULL);
    assert(b->add_edge(g, 0, 1, 1.5));
    assert(b->add_edge(g, 0, 2, 2.5));
    assert(b->add_edge(g, 2, 3, 0.5));
    assert(b->edge_count(g) == 3 && b->degree(g, 0) == 2);

    /* Vertex 0's neighbourhood, in ascending order: {1: 1.5, 2: 2.5}. */
    GraphIter it;
    int v;
    double w;
    b->iter_begin(g, 0, &it);
    assert(b->iter_next(&it, &v, &w) && v == 1 && w == 1.5);
    assert(b->iter_next(&it, &v, &w) && v == 2 && w == 2.5);
    assert(!b->iter_next(&it, &v, &w));

    /* Dead and out-of-range vertices iterate as empty. */
    assert(b->remove_vertex(g, 1));
    b->iter_begin(g, 1, &it);
    assert(!b->iter_next(&it, &v, &w));
    b->iter_begin(g, 99, &it);
    assert(!b->iter_next(&it, &v, &w));

    int order[4];
    int parent[4];
    assert(b->bfs(g, 0, order, parent) == 3);   /* 0, 2, 3 */
    b->destroy(g);
}

/* The formulas must rank the representations the way their structure says. */
static void test_cost_model_shape(void) {
    const GraphBackend *b = matrix_backend();

    GraphShape sparse = {1000, 1000, 2000, false, false};
    GraphShape dense = {1000, 1000, 200000, false, false};

    Workload lookups;
    workload_reset(&lookups);
    lookups.n[OP_HAS_EDGE] = 1000;
    /* Point queries do not depend on density at all. */
    assert(b->estimate_cost(&sparse, &lookups) == b->estimate_cost(&dense, &lookups));

    Workload algos;
    workload_reset(&algos);
    algos.n[OP_ALGO_RUN] = 10;
    assert(b->estimate_cost(&dense, &algos) > b->estimate_cost(&sparse, &algos));

    /* A wider row costs more to walk. */
    GraphShape narrow = {1000, 1000, 5000, false, false};
    GraphShape wide = {1000, 100000, 5000, false, false};
    Workload scans;
    workload_reset(&scans);
    scans.n[OP_NEIGHBOR_SCAN] = 100;
    assert(b->estimate_cost(&wide, &scans) > b->estimate_cost(&narrow, &scans));

    /* Memory is bounded by the tiles an edge set can touch, not by n^2: a
     * sparse graph on a big ID space is far cheaper than a dense one. */
    GraphShape huge_sparse = {50000, 50000, 100000, false, false};
    GraphShape huge_dense = {50000, 50000, 50000000, false, false};
    assert(b->estimate_memory(&huge_sparse) < b->estimate_memory(&huge_dense));
    assert(b->estimate_memory(&dense) > b->estimate_memory(&sparse));

    /* The estimate tracks reality: build a graph, compare to memory_bytes. */
    void *g = b->create(1000, false);
    assert(g != NULL);
    unsigned long long state = 12345;
    for (int i = 0; i < 20000; ++i) {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        const int u = (int)((state >> 33) % 1000);
        const int w = (int)((state >> 13) % 1000);
        (void)b->add_edge(g, u, w, 1.0);
    }
    GraphShape actual = {b->vertex_count(g), b->id_bound(g), b->edge_count(g), false, false};
    const double estimated = (double)b->estimate_memory(&actual);
    const double measured = (double)b->memory_bytes(g);
    assert(estimated > 0.5 * measured && estimated < 2.0 * measured);
    b->destroy(g);
}

int main(void) {
    test_vtable();
    test_cost_model_shape();
    puts("All matrix-backend tests passed.");
    return 0;
}
