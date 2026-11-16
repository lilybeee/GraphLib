/*
 * test_graph_shard.c — the sharding engine (graph_shard.h).
 *
 *   1. Every shard algorithm gives the same answer as the in-memory algorithm, on a graph forced
 *      into many shards with a cache too small to hold them all (so shards really are evicted).
 *   2. The mode switch (graph_open_edge_list), error cases, and cleanup of the shard files.
 *   3. A million users, five million friendships, processed under a hard memory cap.
 */

#define _POSIX_C_SOURCE 200809L

#include <dirent.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <unistd.h>

#include "graph.h"
#include "graph_algorithms.h"
#include "graph_shard.h"
#include "graph_stream.h"

#define CHECK(cond)                                                                  \
    do {                                                                             \
        if (!(cond)) {                                                               \
            fprintf(stderr, "%s:%d: CHECK failed: %s\n", __FILE__, __LINE__, #cond); \
            exit(1);                                                                 \
        }                                                                            \
    } while (0)

#define MB ((size_t)1 << 20)

static bool same_partition(const int *a, const int *b, size_t n)
{
    int *map = malloc(n * sizeof *map), *back = malloc(n * sizeof *back);
    for (size_t i = 0; i < n; i++) map[i] = back[i] = -1;
    bool ok = true;
    for (size_t v = 0; v < n && ok; v++) {
        if (map[a[v]] == -1 && back[b[v]] == -1) { map[a[v]] = b[v]; back[b[v]] = a[v]; }
        else if (map[a[v]] != b[v] || back[b[v]] != a[v]) ok = false;
    }
    free(map); free(back);
    return ok;
}

static size_t files_in(const char *dir)
{
    size_t n = 0;
    DIR *d = opendir(dir);
    if (d == NULL) return 0;
    for (struct dirent *e; (e = readdir(d)) != NULL;) if (e->d_name[0] != '.') n++;
    closedir(d);
    return n;
}

/* A social-like file: hubs from preferential attachment, small islands, and users with no friends. */
static FILE *social_file(size_t n_connected, size_t m, Graph **keep)
{
    Graph *g = graph_gen_barabasi_albert(GRAPH_LIST, n_connected, m, 11);
    CHECK(g != NULL);
    graph_set_adaptive(g, false);
    for (size_t i = 0; i < n_connected / 2; i++) (void)graph_add_vertex(g);
    for (int i = 0; i < 300; i++) {
        const int a = (int)(n_connected + (size_t)i * 5);
        (void)graph_add_edge(g, a, a + 1, 1.0);
    }
    FILE *f = tmpfile();
    CHECK(f != NULL && graph_write_edge_list(g, f) == GRAPH_OK);
    fprintf(f, "# vertices %zu\n", graph_id_bound(g));      /* so trailing silent users exist */
    fflush(f);
    rewind(f);
    if (keep != NULL) *keep = g; else graph_destroy(g);
    return f;
}

static void test_matches_in_memory(void)
{
    FILE *f = social_file(40000, 5, NULL);
    ShardConfig cfg = shard_default_config();
    cfg.memory_budget = 3 * MB;                              /* forces ~20 shards and a small cache */
    ShardedGraph *sg = NULL;
    CHECK(shard_build(f, &cfg, &sg) == GRAPH_OK && sg != NULL);
    ShardStats st0 = shard_stats(sg);
    CHECK(st0.shards > 8 && st0.cache_capacity < st0.shards);   /* the cache cannot hold the graph */

    rewind(f);
    Graph *g = NULL;
    CHECK(graph_read_edge_list(f, false, GRAPH_LIST, &g) == GRAPH_OK);
    const size_t n = graph_id_bound(g);
    CHECK(shard_vertex_count(sg) == n && shard_edge_count(sg) == graph_edge_count(g));

    /* degrees */
    size_t *deg = malloc(n * sizeof *deg);
    CHECK(shard_degrees(sg, deg) == GRAPH_OK);
    for (size_t v = 0; v < n; v++) CHECK(deg[v] == graph_degree(g, (int)v) && deg[v] == shard_degree(sg, (int)v));

    /* components */
    int *c1 = malloc(n * sizeof *c1), *c2 = malloc(n * sizeof *c2);
    size_t k1 = 0;
    CHECK(shard_components(sg, c1, &k1) == GRAPH_OK);
    const size_t k2 = graph_connected_components(g, c2);
    CHECK(k1 == k2 && same_partition(c1, c2, n));

    /* PageRank: same definition, so the same numbers */
    double *r1 = malloc(n * sizeof *r1), *r2 = malloc(n * sizeof *r2);
    size_t it1 = 0, it2 = 0;
    const size_t loads_before = shard_stats(sg).loads;
    CHECK(shard_pagerank(sg, 0.85, 1e-12, 500, r1, &it1) == GRAPH_OK);
    CHECK(graph_pagerank(g, 0.85, 1e-12, 500, r2, &it2) == GRAPH_OK);
    double worst = 0, sum = 0;
    for (size_t v = 0; v < n; v++) { worst = fmax(worst, fabs(r1[v] - r2[v])); sum += r1[v]; }
    CHECK(worst < 1e-9 && fabs(sum - 1.0) < 1e-9);
    /* each round sweeps every shard, and the cache is smaller than the graph, so every round after
     * the first must reload at least the shards that do not fit; alternating the sweep direction
     * keeps it from reloading more than that plus the shards it just left behind */
    const size_t loaded = shard_stats(sg).loads - loads_before, spill = st0.shards - st0.cache_capacity;
    CHECK(loaded >= st0.shards + (it1 - 1) * spill);
    CHECK(loaded < it1 * st0.shards);                        /* better than a full reload every round */

    /* BFS: hop distances equal Dijkstra on unit weights; and depth-limited */
    double *dist = malloc(n * sizeof *dist);
    int *parent = malloc(n * sizeof *parent), *level = malloc(n * sizeof *level);
    for (int src = 0; src < 3; src++) {
        const int source = src * 7919 % 40000;
        CHECK(graph_dijkstra(g, source, dist, parent));
        size_t reached = 0, want = 0;
        CHECK(shard_bfs(sg, source, SIZE_MAX, level, &reached) == GRAPH_OK);
        for (size_t v = 0; v < n; v++) {
            CHECK((dist[v] > 1e300 && level[v] == -1) || (int)dist[v] == level[v]);
            if (level[v] >= 0) want++;
        }
        CHECK(reached == want);
        CHECK(shard_bfs(sg, source, 2, level, &reached) == GRAPH_OK);
        for (size_t v = 0; v < n; v++) CHECK(dist[v] <= 2.0 ? level[v] == (int)dist[v] : level[v] == -1);
    }
    CHECK(shard_bfs(sg, (int)n, 2, level, NULL) == GRAPH_ERR_INVALID_VERTEX);

    free(deg); free(c1); free(c2); free(r1); free(r2); free(dist); free(parent); free(level);
    graph_destroy(g);
    shard_destroy(sg);
    fclose(f);
}

static void test_mode_switch_and_errors(void)
{
    FILE *f = social_file(5000, 3, NULL);

    /* a generous budget loads an ordinary graph... */
    ShardConfig cfg = shard_default_config();
    cfg.memory_budget = (size_t)512 * MB;
    Graph *g = NULL;
    ShardedGraph *sg = NULL;
    CHECK(graph_open_edge_list(f, &cfg, &g, &sg) == GRAPH_OK && g != NULL && sg == NULL);
    const size_t n = graph_id_bound(g), edges = graph_edge_count(g);
    graph_destroy(g);

    /* ...unless the user switches sharding on... */
    cfg.force = true;
    CHECK(graph_open_edge_list(f, &cfg, &g, &sg) == GRAPH_OK && g == NULL && sg != NULL);
    CHECK(shard_vertex_count(sg) == n && shard_edge_count(sg) == edges);
    shard_destroy(sg);

    /* ...and a budget too small for the in-memory list switches it on by itself */
    cfg.force = false;
    cfg.memory_budget = 512 * 1024;                          /* the list would need ~790 KB */
    CHECK(graph_open_edge_list(f, &cfg, &g, &sg) == GRAPH_OK && g == NULL && sg != NULL);
    shard_destroy(sg);

    /* a budget that cannot even hold the per-user state is refused, and leaves nothing behind */
    char dir[] = "/tmp/shard_test_XXXXXX";
    CHECK(mkdtemp(dir) != NULL);
    cfg.memory_budget = 64 * 1024;                           /* 64 KB: 12 000 users need ~290 KB */
    cfg.dir = dir;
    CHECK(shard_build(f, &cfg, &sg) == GRAPH_ERR_TOO_LARGE && sg == NULL);
    CHECK(files_in(dir) == 0);

    /* a normal build writes files into the chosen directory and removes every one of them */
    cfg.memory_budget = 4 * MB;
    CHECK(shard_build(f, &cfg, &sg) == GRAPH_OK);
    CHECK(files_in(dir) == shard_stats(sg).shards && files_in(dir) > 0);
    shard_destroy(sg);
    CHECK(files_in(dir) == 0);
    rmdir(dir);

    /* malformed input, bad arguments */
    FILE *bad = tmpfile();
    fputs("0 1\nnot an edge\n", bad);
    cfg.dir = NULL;
    CHECK(shard_build(bad, &cfg, &sg) == GRAPH_ERR_IO && sg == NULL);
    fclose(bad);
    CHECK(shard_build(NULL, &cfg, &sg) == GRAPH_ERR_INVALID_ARGUMENT);
    CHECK(shard_degree(NULL, 0) == 0 && shard_vertex_count(NULL) == 0);
    shard_destroy(NULL);

    /* repeated edges are merged exactly as the in-memory loader merges them */
    FILE *dup = tmpfile();
    fputs("0 1\n1 0\n0 1\n2 2\n2 2\n3 1\n", dup);
    CHECK(shard_build(dup, &cfg, &sg) == GRAPH_OK);
    rewind(dup);
    Graph *dg = NULL;
    CHECK(graph_read_edge_list(dup, false, GRAPH_LIST, &dg) == GRAPH_OK);
    CHECK(shard_edge_count(sg) == graph_edge_count(dg) && shard_edge_count(sg) == 3);   /* 0-1, 2-2, 3-1 */
    for (int v = 0; v < 4; v++) CHECK(shard_degree(sg, v) == graph_degree(dg, v));
    graph_destroy(dg);
    shard_destroy(sg);
    fclose(dup);

    /* the empty graph and a graph with a single edge */
    FILE *empty = tmpfile();
    CHECK(shard_build(empty, &cfg, &sg) == GRAPH_OK && shard_vertex_count(sg) == 0);
    CHECK(shard_components(sg, NULL, NULL) == GRAPH_OK);
    shard_destroy(sg);
    fclose(empty);
    FILE *one = tmpfile();
    fputs("# directed\n0 3\n", one);
    CHECK(shard_build(one, &cfg, &sg) == GRAPH_OK && shard_vertex_count(sg) == 4 && shard_is_directed(sg));
    int level[4];
    size_t reached = 0;
    CHECK(shard_bfs(sg, 0, SIZE_MAX, level, &reached) == GRAPH_OK && reached == 2 && level[3] == 1 && level[1] == -1);
    CHECK(shard_bfs(sg, 3, SIZE_MAX, level, &reached) == GRAPH_OK && reached == 1);   /* arcs only go one way */
    shard_destroy(sg);
    fclose(one);
    fclose(f);
}

/* A million users and five million friendships on a "device" with 200 MB of address space and a
 * 64 MB budget.  Holding that graph as a list needs ~230 MB, so it could not be loaded at all. */
static void test_million_users_small_device(void)
{
    const size_t n = 1000000, m = 5000000;
    FILE *f = tmpfile();
    CHECK(f != NULL);
    uint64_t x = 88172645463325252ULL;
    size_t loops = 0;
    for (size_t i = 0; i < m; i++) {
        x ^= x << 13; x ^= x >> 7; x ^= x << 17;
        const unsigned u = (unsigned)(x % (n * 8 / 10)), v = (unsigned)((x >> 24) % (n * 8 / 10));
        loops += u == v;
        fprintf(f, "%u %u\n", u, v);
    }
    fprintf(f, "# vertices %zu\n", n);
    fflush(f);
    rewind(f);

    struct rlimit lim = {(rlim_t)200 * MB, (rlim_t)200 * MB};
    CHECK(setrlimit(RLIMIT_AS, &lim) == 0);

    ShardConfig cfg = shard_default_config();
    cfg.memory_budget = 64 * MB;
    Graph *g = NULL;
    ShardedGraph *sg = NULL;
    CHECK(graph_open_edge_list(f, &cfg, &g, &sg) == GRAPH_OK);
    CHECK(g == NULL && sg != NULL);                           /* too big for the budget: sharded */
    CHECK(shard_vertex_count(sg) == n);
    CHECK(shard_edge_count(sg) <= m && shard_edge_count(sg) > m - m / 100);   /* repeats merged: a few lost */

    size_t *deg = malloc(n * sizeof *deg);
    CHECK(deg != NULL && shard_degrees(sg, deg) == GRAPH_OK);
    size_t sum = 0;
    for (size_t v = 0; v < n; v++) sum += deg[v];
    CHECK(sum <= 2 * shard_edge_count(sg) && sum + loops >= 2 * shard_edge_count(sg));   /* a loop counts once */
    free(deg);

    int *comp = malloc(n * sizeof *comp);
    size_t count = 0;
    CHECK(comp != NULL && shard_components(sg, comp, &count) == GRAPH_OK && count > n / 5);   /* 20% have none */
    rewind(f);
    unsigned a, b;
    size_t checked = 0;
    while (fscanf(f, "%u %u", &a, &b) == 2) { CHECK(comp[a] == comp[b]); checked++; }
    CHECK(checked == m);

    double *rank = malloc(n * sizeof *rank);
    size_t iters = 0;
    CHECK(rank != NULL && shard_pagerank(sg, 0.85, 1e-8, 100, rank, &iters) == GRAPH_OK && iters > 3);
    double total = 0, top = 0;
    for (size_t v = 0; v < n; v++) { total += rank[v]; top = fmax(top, rank[v]); }
    CHECK(fabs(total - 1.0) < 1e-6 && top > 1.0 / (double)n);

    int *level = malloc(n * sizeof *level);
    size_t reached = 0;
    CHECK(level != NULL && shard_bfs(sg, 0, 3, level, &reached) == GRAPH_OK && reached > 1);
    for (size_t v = 0; v < n; v++) CHECK(level[v] <= 3);

    const ShardStats st = shard_stats(sg);
    CHECK(st.loads > st.shards);                              /* shards were evicted and reloaded */
    printf("  1,000,000 users on a 64 MB budget: %zu shards, cache holds %zu, %zu loads, %.0f MB read\n",
           st.shards, st.cache_capacity, st.loads, (double)st.bytes_read / (double)MB);
    free(comp); free(rank); free(level);
    shard_destroy(sg);
    fclose(f);
}

int main(int argc, char **argv)
{
    test_matches_in_memory();
    test_mode_switch_and_errors();
    if (argc < 2 || strcmp(argv[1], "quick") != 0) {       /* "quick" skips the million-user run (for valgrind) */
        test_million_users_small_device();
    }
    printf("All graph shard tests passed.\n");
    return 0;
}
