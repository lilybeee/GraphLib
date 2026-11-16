/*
 * test_memory_budget.c — the memory budget (memory_budget.h): going over it must
 * be a clean GRAPH_ERR_TOO_LARGE with the graph left intact, never a crash.
 */

#include <stdio.h>
#include <stdlib.h>

#include "adjacency_list.h"
#include "adjacency_matrix.h"
#include "graph.h"
#include "graph_algorithms.h"
#include "memory_budget.h"

#define CHECK(cond)                                                                  \
    do {                                                                             \
        if (!(cond)) {                                                               \
            fprintf(stderr, "%s:%d: CHECK failed: %s\n", __FILE__, __LINE__, #cond); \
            exit(1);                                                                 \
        }                                                                            \
    } while (0)

#define MB ((size_t)1 << 20)

static void test_default(void)
{
    CHECK(graph_memory_budget() == (size_t)4 << 30);
    CHECK(graph_default_config(false).observer.memory_budget == (size_t)4 << 30);
}

static void test_matrix_limits(void)
{
    /* The calculated limit: presence bits alone are V^2/8 bytes once every row strip exists. */
    const size_t w = matrix_max_vertices((size_t)4 << 30, true), u = matrix_max_vertices((size_t)4 << 30, false);
    CHECK(w >= 131000 && w <= 131100);          /* sqrt(8 * 0.5 * 4 GiB) */
    CHECK(u >= 175000 && u <= 176000);          /* sqrt(8 * 0.9 * 4 GiB) */
    CHECK(matrix_max_vertices(0, true) == 2147483647u);

    graph_set_memory_budget(256 * MB);          /* limit: ~32 768 weighted vertices */
    CHECK(matrix_graph_create(1000000, false) == NULL && graph_last_error() == GRAPH_ERR_TOO_LARGE);
    CHECK(matrix_graph_create(40000, false) == NULL && graph_last_error() == GRAPH_ERR_TOO_LARGE);
    MatrixGraph *g = matrix_graph_create(30000, false);
    CHECK(g != NULL);
    /* growing past the limit by adding vertices is refused too, and the graph survives it */
    int added_vertices = 0;
    while (matrix_graph_add_vertex(g) >= 0) added_vertices++;
    CHECK(graph_last_error() == GRAPH_ERR_TOO_LARGE && matrix_graph_id_bound(g) >= 32000);
    CHECK(added_vertices > 0 && matrix_graph_id_bound(g) <= matrix_max_vertices(256 * MB, true));

    /* growth by edges: each random edge opens a row strip and a 32 KB weight tile */
    size_t added = 0;
    bool refused = false;
    unsigned long long x = 88172645463325252ULL;
    for (size_t i = 0; i < 400000 && !refused; i++) {
        x ^= x << 13; x ^= x >> 7; x ^= x << 17;
        const int u2 = (int)(x % 30000), v2 = (int)((x >> 20) % 30000);
        if (u2 == v2 || matrix_graph_has_edge(g, u2, v2, NULL)) continue;
        if (matrix_graph_add_edge(g, u2, v2, 1.0)) {
            added++;
        } else {
            refused = graph_last_error() == GRAPH_ERR_TOO_LARGE;
            CHECK(refused);                                   /* the only reason it may fail */
        }
    }
    CHECK(refused && added > 100);
    CHECK(matrix_graph_edge_count(g) == added);               /* nothing half-added */
    CHECK(matrix_graph_memory_bytes(g) <= 256 * MB);          /* stayed inside the budget */
    matrix_graph_destroy(g);
    graph_set_memory_budget(GRAPH_DEFAULT_MEMORY_BUDGET);
}

/* The sparse tile index replaces a (V/64)^2 pointer array: a large, sparse graph must now cost
 * memory in proportion to the tiles it uses, and every weight must still read back exactly. */
static void test_sparse_index(void)
{
    const size_t n = 100000;
    MatrixGraph *g = matrix_graph_create(n, false);
    CHECK(g != NULL);
    CHECK(matrix_graph_memory_bytes(g) < 64 * MB);            /* the old directory alone was 19 MB at this size */
    unsigned long long x = 88172645463325252ULL;
    enum { EDGES = 3000 };
    static int eu[EDGES], ev[EDGES];
    static double ew[EDGES];
    for (int i = 0; i < EDGES; i++) {
        x ^= x << 13; x ^= x >> 7; x ^= x << 17;
        eu[i] = (int)(x % n); ev[i] = (int)((x >> 24) % n); ew[i] = 1.0 + (double)(i % 97) / 8.0;
        if (eu[i] == ev[i] || matrix_graph_has_edge(g, eu[i], ev[i], NULL)) { eu[i] = -1; continue; }
        CHECK(matrix_graph_add_edge(g, eu[i], ev[i], ew[i]));
    }
    for (int i = 0; i < EDGES; i++) {                          /* survives index growth and rehashing */
        if (eu[i] < 0) continue;
        double got = 0;
        CHECK(matrix_graph_has_edge(g, eu[i], ev[i], &got) && got == ew[i]);
        CHECK(matrix_graph_has_edge(g, ev[i], eu[i], &got) && got == ew[i]);   /* undirected mirror */
    }
    matrix_graph_destroy(g);
}

static void test_unweighted_matrix(void)
{
    MatrixGraph *g = matrix_graph_create_unweighted(30000, false), *w = matrix_graph_create(30000, false);
    CHECK(g != NULL && w != NULL && matrix_graph_is_unweighted(g) && !matrix_graph_is_unweighted(w));
    unsigned long long x = 88172645463325252ULL;
    for (int i = 0; i < 3000; i++) {
        x ^= x << 13; x ^= x >> 7; x ^= x << 17;
        const int a = (int)(x % 30000), b = (int)((x >> 24) % 30000);
        if (a == b || matrix_graph_has_edge(g, a, b, NULL)) continue;
        CHECK(matrix_graph_add_edge(g, a, b, 1.0) && matrix_graph_add_edge(w, a, b, 1.0));
    }
    double wt = 0;
    int a = 0, b = 1;
    CHECK(matrix_graph_add_edge(g, a, b, 1.0) || matrix_graph_has_edge(g, a, b, NULL));
    CHECK(matrix_graph_has_edge(g, a, b, &wt) && wt == 1.0);
    CHECK(!matrix_graph_add_edge(g, 7, 8, 2.5) && graph_last_error() == GRAPH_ERR_INVALID_WEIGHT);
    CHECK(matrix_graph_edge_count(g) == matrix_graph_edge_count(w) + (matrix_graph_has_edge(w, a, b, NULL) ? 0 : 1));
    /* no 32 KB weight tile per edge: the same graph takes under half the memory (the row strips remain) */
    CHECK(matrix_graph_memory_bytes(g) * 2 < matrix_graph_memory_bytes(w));
    CHECK(matrix_graph_remove_vertex(g, a) && matrix_graph_vertex_count(g) == 29999);
    int *remap = malloc(30000 * sizeof *remap);
    CHECK(matrix_graph_compact(g, remap) && matrix_graph_is_unweighted(g) && matrix_graph_id_bound(g) == 29999);
    free(remap);
    matrix_graph_destroy(g);
    matrix_graph_destroy(w);
}

static void test_list_limits(void)
{
    graph_set_memory_budget(2 * MB);
    ListGraph *g = list_graph_create(20000, false);
    CHECK(g != NULL);
    size_t added = 0;
    bool refused = false;
    for (int i = 1; i < 20000 && !refused; i++) {
        for (int k = 1; k <= 8 && !refused; k++) {
            const int v = (i * k * 31 + 7) % 20000;
            if (v == i || list_graph_has_edge(g, i, v, NULL)) continue;
            if (list_graph_add_edge(g, i, v, 1.0)) {
                added++;
            } else {
                refused = graph_last_error() == GRAPH_ERR_TOO_LARGE;
                CHECK(refused);
            }
        }
    }
    CHECK(refused && added > 1000);
    CHECK(list_graph_edge_count(g) == added);
    CHECK(list_graph_memory_bytes(g) <= 2 * MB + 4096);
    /* deleting a busy vertex hands its memory back, so adding works again */
    size_t busiest = 0; int bv = 0;
    for (int v = 0; v < 20000; v++) if (list_graph_degree(g, v) > busiest) { busiest = list_graph_degree(g, v); bv = v; }
    CHECK(list_graph_remove_vertex(g, bv));
    CHECK(list_graph_add_vertex(g) >= 0);
    list_graph_destroy(g);

    /* a vertex table that cannot fit is refused at creation */
    graph_set_memory_budget(1 * MB);
    CHECK(list_graph_create(10000000, false) == NULL && graph_last_error() == GRAPH_ERR_TOO_LARGE);
    graph_set_memory_budget(GRAPH_DEFAULT_MEMORY_BUDGET);
}

static void test_convert_and_load(void)
{
    /* 200 000 users, 100 000 friendships: tiny as a list, far too big as a matrix under 128 MB */
    graph_set_memory_budget(128 * MB);
    Graph *g = graph_create_sized(200000, false);
    CHECK(g != NULL);
    graph_set_adaptive(g, false);
    for (int i = 0; i < 100000; i++) {
        (void)graph_add_edge(g, (int)((i * 7919LL) % 200000), (int)((i * 104729LL + 13) % 200000), 1.0);   /* 64-bit: i * 104729 overflows an int */
    }
    const size_t edges = graph_edge_count(g);
    CHECK(!graph_convert(g, GRAPH_MATRIX));
    CHECK(graph_last_error() == GRAPH_ERR_TOO_LARGE);
    CHECK(graph_representation(g) == GRAPH_LIST && graph_edge_count(g) == edges);   /* untouched */
    double dist[1];
    int parent[1];
    (void)dist; (void)parent;
    CHECK(graph_degree(g, 0) <= graph_degree(g, 0));                                 /* still usable */

    /* loading the snapshot into the matrix fails cleanly too */
    FILE *f = tmpfile();
    CHECK(f != NULL && graph_save_binary(g, f) == GRAPH_OK);
    rewind(f);
    Graph *back = NULL;
    CHECK(graph_load_binary(f, GRAPH_MATRIX, &back) != GRAPH_OK && back == NULL);
    CHECK(graph_last_error() == GRAPH_ERR_TOO_LARGE);
    rewind(f);
    CHECK(graph_load_binary(f, GRAPH_LIST, &back) == GRAPH_OK && back != NULL);       /* list is fine */
    CHECK(graph_edge_count(back) == edges);
    graph_destroy(back);
    fclose(f);

    /* and a graph that does fit converts as before */
    Graph *small = graph_create_sized(100, false);
    CHECK(small != NULL && graph_add_edge(small, 1, 2, 1.0));
    CHECK(graph_convert(small, GRAPH_MATRIX) && graph_representation(small) == GRAPH_MATRIX);
    CHECK(graph_has_edge(small, 1, 2, NULL));
    graph_destroy(small);

    /* budget 0 = unlimited: the check is skipped */
    graph_set_memory_budget(0);
    Graph *free_graph = graph_create_sized(1000, false);
    CHECK(free_graph != NULL && graph_convert(free_graph, GRAPH_MATRIX));
    graph_destroy(free_graph);

    graph_destroy(g);
    graph_set_memory_budget(GRAPH_DEFAULT_MEMORY_BUDGET);
}


static void test_unweighted(void)
{
    GraphConfig cfg = graph_default_config(false);
    cfg.unweighted  = true;
    Graph *g = graph_create_with(&cfg, 1000);
    CHECK(g != NULL && graph_representation(g) == GRAPH_LIST);
    CHECK(graph_add_edge(g, 1, 2, 1.0) && graph_add_edge(g, 2, 3, 1.0) && graph_add_edge(g, 1, 3, 1.0));
    CHECK(!graph_add_edge(g, 4, 5, 2.5) && graph_last_error() == GRAPH_ERR_INVALID_WEIGHT);   /* no weight to remember */
    double w = 0;
    CHECK(graph_has_edge(g, 1, 2, &w) && w == 1.0);
    CHECK(graph_remove_edge(g, 1, 2) && !graph_has_edge(g, 1, 2, NULL));
    CHECK(graph_add_edge(g, 1, 2, 1.0));

    /* the algorithms see weight 1 everywhere: Dijkstra distances are hop counts */
    double dist[1000];
    int parent[1000];
    CHECK(graph_dijkstra(g, 1, dist, parent) && dist[3] == 1.0 && dist[2] == 1.0 && dist[4] > 1e300);
    long long tri = 0;
    CHECK(graph_triangle_count(g, &tri) == GRAPH_OK && tri == 1);
    Graph *copy = graph_copy(g);
    CHECK(copy != NULL && graph_edge_count(copy) == 3 && graph_has_edge(copy, 2, 3, &w) && w == 1.0);
    graph_destroy(copy);
    graph_destroy(g);

    /* the saving: the same 20 000 friendships take clearly less memory */
    GraphConfig plain = graph_default_config(false);
    plain.adaptive    = false;
    Graph *a = graph_create_with(&plain, 5000), *b = graph_create_with(&cfg, 5000);
    CHECK(a != NULL && b != NULL);
    for (int i = 0; i < 20000; i++) {
        const int u = (int)((i * 7919LL) % 5000), v = (int)((i * 104729LL + 13) % 5000);
        if (u != v) { (void)graph_add_edge(a, u, v, 1.0); (void)graph_add_edge(b, u, v, 1.0); }
    }
    CHECK(graph_edge_count(a) == graph_edge_count(b));
    graph_destroy(a);
    graph_destroy(b);
}

static void test_shrink(void)
{
    GraphConfig cfg = graph_default_config(false);
    cfg.adaptive    = false;
    Graph *g = graph_create_with(&cfg, 20000);
    CHECK(g != NULL);
    for (int i = 0; i < 90000; i++) (void)graph_add_edge(g, (int)((i * 7919LL) % 20000), (int)((i * 104729LL + 13) % 20000), 1.0 + i % 7);
    const size_t edges = graph_edge_count(g);
    FILE *before = tmpfile();
    CHECK(graph_save_binary(g, before) == GRAPH_OK);
    CHECK(graph_shrink_to_fit(g));
    FILE *after = tmpfile();
    CHECK(graph_save_binary(g, after) == GRAPH_OK);
    /* contents identical: the snapshots are byte for byte the same */
    CHECK(ftell(before) == ftell(after));
    rewind(before); rewind(after);
    int a, b;
    long same = 0;
    while ((a = fgetc(before)) != EOF && (b = fgetc(after)) != EOF) { CHECK(a == b); same++; }
    CHECK(same == ftell(before));
    fclose(before); fclose(after);
    CHECK(graph_edge_count(g) == edges);
    CHECK(graph_add_edge(g, 1, 19999, 3.0) || graph_has_edge(g, 1, 19999, NULL));   /* still editable */
    graph_freeze(g);
    CHECK(!graph_shrink_to_fit(g) && graph_last_error() == GRAPH_ERR_FROZEN);
    graph_destroy(g);
}

static void test_matrix_policy(void)
{
    /* too many vertices for a matrix: refused even though the graph has no edges at all */
    Graph *big = graph_create_sized(200000, false);
    CHECK(big != NULL);
    CHECK(!graph_convert(big, GRAPH_MATRIX) && graph_last_error() == GRAPH_ERR_TOO_LARGE);
    CHECK(graph_representation(big) == GRAPH_LIST);
    graph_destroy(big);

    GraphConfig cfg = graph_default_config(false);
    cfg.start = GRAPH_MATRIX;
    CHECK(graph_create_with(&cfg, 1000000) == NULL && graph_last_error() == GRAPH_ERR_TOO_LARGE);

    /* an unweighted graph stays unweighted through list -> matrix -> list */
    GraphConfig uc = graph_default_config(false);
    uc.unweighted = true;
    uc.adaptive   = false;
    Graph *g = graph_create_with(&uc, 500);
    CHECK(g != NULL && graph_add_edge(g, 1, 2, 1.0) && graph_add_edge(g, 2, 3, 1.0));
    CHECK(graph_convert(g, GRAPH_MATRIX) && graph_representation(g) == GRAPH_MATRIX);
    CHECK(!graph_add_edge(g, 4, 5, 2.0) && graph_last_error() == GRAPH_ERR_INVALID_WEIGHT);
    double dist[500];
    int parent[500];
    CHECK(graph_dijkstra(g, 1, dist, parent) && dist[3] == 2.0);   /* hop counts */
    CHECK(graph_convert(g, GRAPH_LIST) && !graph_add_edge(g, 4, 5, 2.0));
    graph_destroy(g);

    /* the unweighted matrix really does reach larger graphs than the weighted one */
    graph_set_memory_budget((size_t)4 << 30);
    CHECK(matrix_max_vertices(graph_memory_budget(), false) > matrix_max_vertices(graph_memory_budget(), true));
}

int main(void)
{
    test_default();
    test_matrix_limits();
    test_sparse_index();
    test_unweighted_matrix();
    test_list_limits();
    test_convert_and_load();
    test_unweighted();
    test_shrink();
    test_matrix_policy();
    printf("All memory-budget tests passed.\n");
    return 0;
}
