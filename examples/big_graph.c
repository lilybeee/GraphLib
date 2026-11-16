/*
 * big_graph.c — working with big graphs on a machine with limited memory.
 *
 * Shows the four tools the library gives you for that, in the order you would reach for them:
 *   1. the memory budget (a request that cannot fit is refused, it does not crash),
 *   2. unweighted graphs (half the memory when edges carry no weight),
 *   3. shrink-to-fit (give back growth slack once the graph is built),
 *   4. streaming (analyse an edge-list file without ever holding the graph).
 *
 *   make examples
 */

#include <stdio.h>
#include <stdlib.h>

#include "graph_algorithms.h"
#include "graph_stream.h"

static unsigned long long seed = 88172645463325252ULL;
static unsigned rnd(unsigned bound) { seed ^= seed << 13; seed ^= seed >> 7; seed ^= seed << 17; return (unsigned)(seed % bound); }

/* A social network: every user befriends a few random users; every fifth user has no friends. */
static void befriend(Graph *g, size_t users)
{
    for (size_t u = 0; u < users; u++) {
        if (u % 5 == 4) continue;
        for (int k = 0; k < 5; k++) {
            graph_add_edge(g, (int)u, (int)rnd((unsigned)users), 1.0);
        }
    }
}

int main(void)
{
    const size_t users = 200000;

    /* 1. The budget.  Every graph may use at most graph_memory_budget() bytes (4 GiB by default;
     *    graph_set_memory_budget changes it for graphs created afterwards).  Asking for more is
     *    an error you can handle. */
    printf("memory budget: %zu MiB\n", graph_memory_budget() >> 20);
    Graph *g = graph_create_sized(users, false);
    befriend(g, users);
    printf("%zu users, %zu friendships, stored as a %s\n", graph_vertex_count(g), graph_edge_count(g),
           graph_rep_name(graph_representation(g)));
    if (!graph_convert(g, GRAPH_MATRIX)) {
        printf("convert to matrix refused: %s\n", graph_status_string(graph_last_error()));
    }
    printf("still a %s, still usable: user 0 has %zu friends\n",
           graph_rep_name(graph_representation(g)), graph_degree(g, 0));
    printf("largest matrix the budget allows: %zu vertices weighted, %zu unweighted\n",
           graph_matrix_max_vertices(true), graph_matrix_max_vertices(false));

    /* 2. Unweighted: friendships carry no weight, so do not store one.  Every edge reads back
     *    as 1.0, and any other weight is refused instead of silently dropped. */
    GraphConfig cfg = graph_default_config(false);
    cfg.unweighted  = true;
    Graph *lean = graph_create_with(&cfg, users);
    seed = 88172645463325252ULL;
    befriend(lean, users);
    printf("\nunweighted copy: %zu friendships\n", graph_edge_count(lean));
    if (!graph_add_edge(lean, 1, 2, 2.5)) {
        printf("add_edge with weight 2.5: %s\n", graph_status_string(graph_last_error()));
    }

    /* 3. Shrink-to-fit: once the graph is built, release the spare room each friend list keeps
     *    for growth.  Contents are unchanged. */
    graph_shrink_to_fit(g);
    graph_shrink_to_fit(lean);

    /* 4. Streaming: write the graph as an edge list, then analyse the FILE.  Only a few integers
     *    per user are kept; no edge is ever stored, so the file may be far larger than memory. */
    FILE *file = tmpfile();
    graph_write_edge_list(g, file);
    rewind(file);
    GraphStreamInfo info;
    graph_stream_scan(file, &info);
    printf("\nstreamed file: %zu users, %zu friendships\n", info.vertices, info.edges);

    rewind(file);
    int *component = malloc(info.vertices * sizeof *component);
    size_t groups;
    graph_stream_components(file, info.vertices, component, &groups);
    printf("connected groups found while streaming: %zu\n", groups);

    rewind(file);
    size_t *degree = malloc(info.vertices * sizeof *degree);
    graph_stream_degrees(file, info.vertices, false, degree, NULL);
    size_t isolated = 0, busiest = 0;
    for (size_t v = 0; v < info.vertices; v++) {
        if (degree[v] == 0) isolated++;
        if (degree[v] > busiest) busiest = degree[v];
    }
    printf("users with no friends: %zu, most friends of any user: %zu\n", isolated, busiest);

    free(component);
    free(degree);
    fclose(file);
    graph_destroy(g);
    graph_destroy(lean);
    return 0;
}
