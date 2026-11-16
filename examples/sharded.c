/*
 * sharded.c — running on a device with very little memory.
 *
 * We pretend to be a small machine: the memory budget is lowered to 10 MiB, then a 200 000-user
 * social network (an edge-list FILE) is opened.  The ordinary in-memory graph would need about
 * 35 MiB, so graph_open_edge_list switches the sharding engine on by itself: the friendships stay
 * on disk in shard files and every algorithm becomes a sweep over them.  The same file is then
 * loaded the ordinary way (budget lifted) to show that both give the same answers, and what the
 * low-memory mode costs in time.
 *
 * Build and run:  make build/examples/sharded && ./build/examples/sharded
 */

#define _POSIX_C_SOURCE 200809L

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "graph.h"
#include "graph_algorithms.h"
#include "graph_shard.h"

static double now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec + (double)t.tv_nsec * 1e-9;
}

/* A made-up social network written as an edge list: popular users attract friends, 20% of the
 * users never connect to anyone. */
static FILE *make_network(size_t users)
{
    FILE *f = tmpfile();
    uint64_t x = 88172645463325252ULL;
    int *ends = malloc(users * 24 * sizeof *ends);
    size_t en = 0;
    for (size_t u = 0; u < users; u++) {
        if (u % 5 == 4) continue;                                  /* the silent fifth */
        x ^= x << 13; x ^= x >> 7; x ^= x << 17;
        const size_t links = 1 + (size_t)(x % 8);
        for (size_t k = 0; k < links && en > 0; k++) {
            x ^= x << 13; x ^= x >> 7; x ^= x << 17;
            const int friend_id = (x % 10 < 8) ? ends[x / 10 % en] : (int)(x / 10 % u);
            if (friend_id != (int)u && friend_id != (int)u + 0 && (size_t)friend_id % 5 != 4) {
                fprintf(f, "%zu %d\n", u, friend_id);
                ends[en++] = (int)u;
                ends[en++] = friend_id;
            }
        }
        if (en == 0) { ends[en++] = (int)u; ends[en++] = (int)u; }
    }
    fprintf(f, "# vertices %zu\n", users);
    free(ends);
    fflush(f);
    rewind(f);
    return f;
}

int main(void)
{
    const size_t users = 200000;
    FILE *file = make_network(users);

    /* ---- the low-memory mode ---- */
    graph_set_memory_budget((size_t)10 << 20);                     /* a small device: 10 MiB */
    ShardConfig cfg = shard_default_config();                      /* budget = the one just set */
    Graph *g = NULL;
    ShardedGraph *sg = NULL;
    if (graph_open_edge_list(file, &cfg, &g, &sg) != GRAPH_OK) {
        printf("open failed: %s\n", graph_status_string(graph_last_error()));
        return 1;
    }
    printf("budget 10 MiB -> %s\n", sg != NULL ? "sharded (the graph would not fit in memory)"
                                               : "loaded in memory");
    ShardStats st = shard_stats(sg);
    printf("%zu users, %zu friendships in %zu shard files; the cache holds %zu at a time\n\n",
           shard_vertex_count(sg), shard_edge_count(sg), st.shards, st.cache_capacity);

    int *component = malloc(users * sizeof *component);
    double *rank = malloc(users * sizeof *rank);
    int *level = malloc(users * sizeof *level);
    size_t groups = 0, rounds = 0, reached = 0;

    double t = now();
    shard_components(sg, component, &groups);
    printf("connected groups: %zu  (%.2f s)\n", groups, now() - t);

    t = now();
    shard_pagerank(sg, 0.85, 1e-9, 100, rank, &rounds);
    size_t best = 0;
    for (size_t v = 1; v < users; v++) if (rank[v] > rank[best]) best = v;
    printf("most influential user: %zu (PageRank %.5f, %zu friends), %zu rounds  (%.2f s)\n", best, rank[best],
           shard_degree(sg, (int)best), rounds, now() - t);

    t = now();
    shard_bfs(sg, (int)best, 2, level, &reached);
    printf("people within 2 hops of them: %zu  (%.2f s)\n", reached - 1, now() - t);

    st = shard_stats(sg);
    printf("disk traffic: %zu shard loads, %.0f MB read\n\n", st.loads, (double)st.bytes_read / 1048576.0);

    /* ---- the same file, the ordinary way: lift the budget and compare ---- */
    graph_set_memory_budget((size_t)4 << 30);
    cfg.memory_budget = (size_t)4 << 30;
    Graph *mem = NULL;
    ShardedGraph *none = NULL;
    graph_open_edge_list(file, &cfg, &mem, &none);
    printf("budget 4 GiB  -> %s\n", mem != NULL ? "loaded in memory" : "sharded");

    double *rank2 = malloc(users * sizeof *rank2);
    size_t rounds2 = 0;
    t = now();
    graph_pagerank(mem, 0.85, 1e-9, 100, rank2, &rounds2);
    double worst = 0;
    for (size_t v = 0; v < users; v++) worst = fmax(worst, fabs(rank[v] - rank2[v]));
    printf("in-memory PageRank: %zu rounds (%.2f s); largest difference from the sharded result: %.1e\n",
           rounds2, now() - t, worst);

    free(component); free(rank); free(rank2); free(level);
    graph_destroy(mem);
    shard_destroy(sg);
    fclose(file);
    return 0;
}
