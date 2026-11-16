/*
 * test_graph_stream.c — graph_stream.h: analyses of an edge-list file that never
 * build the graph.  Checked against the in-memory algorithms, on malformed input,
 * and on a million-user file processed under a hard memory cap.
 */

#define _POSIX_C_SOURCE 200809L

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>

#include "graph.h"
#include "graph_algorithms.h"
#include "graph_stream.h"

#define CHECK(cond)                                                                  \
    do {                                                                             \
        if (!(cond)) {                                                               \
            fprintf(stderr, "%s:%d: CHECK failed: %s\n", __FILE__, __LINE__, #cond); \
            exit(1);                                                                 \
        }                                                                            \
    } while (0)

static FILE *file_with(const char *text)
{
    FILE *f = tmpfile();
    CHECK(f != NULL && fputs(text, f) >= 0);
    rewind(f);
    return f;
}

static void test_small(void)
{
    FILE *f = file_with("# undirected\n# vertices 8\n\n0 1\n1 2 2.5\n# a comment\n3 4\n2 0\n");
    GraphStreamInfo info;
    CHECK(graph_stream_scan(f, &info) == GRAPH_OK);
    CHECK(info.vertices == 8 && info.edges == 4 && !info.directed);   /* 5,6,7 exist, in no edge */
    rewind(f);
    size_t deg[8];
    CHECK(graph_stream_degrees(f, 8, false, deg, NULL) == GRAPH_OK);
    CHECK(deg[0] == 2 && deg[1] == 2 && deg[2] == 2 && deg[3] == 1 && deg[4] == 1 && deg[5] == 0 && deg[7] == 0);
    rewind(f);
    int comp[8];
    size_t count = 0;
    CHECK(graph_stream_components(f, 8, comp, &count) == GRAPH_OK);
    CHECK(count == 5);                                                 /* {0,1,2} {3,4} {5} {6} {7} */
    CHECK(comp[0] == 0 && comp[1] == 0 && comp[2] == 0 && comp[3] == 1 && comp[4] == 1 && comp[5] == 2 &&
          comp[6] == 3 && comp[7] == 4);
    fclose(f);

    /* directed: out- and in-degree */
    f = file_with("# directed\n0 1\n0 2\n2 1\n");
    CHECK(graph_stream_scan(f, &info) == GRAPH_OK && info.directed && info.vertices == 3);
    rewind(f);
    size_t out[3], in[3];
    CHECK(graph_stream_degrees(f, 3, true, out, in) == GRAPH_OK);
    CHECK(out[0] == 2 && out[1] == 0 && out[2] == 1 && in[0] == 0 && in[1] == 2 && in[2] == 1);
    fclose(f);

    /* malformed input is reported, not guessed at */
    f = file_with("0 1\nhello world\n");
    CHECK(graph_stream_scan(f, &info) == GRAPH_ERR_IO);
    fclose(f);
    f = file_with("0 1\n1 9\n");
    CHECK(graph_stream_degrees(f, 4, false, deg, NULL) == GRAPH_ERR_INVALID_VERTEX);   /* user 9 >= n */
    fclose(f);
    char longline[2000];
    memset(longline, '1', sizeof longline - 1);
    longline[sizeof longline - 1] = '\0';
    f = file_with(longline);
    CHECK(graph_stream_scan(f, &info) == GRAPH_ERR_IO);
    fclose(f);
    CHECK(graph_stream_scan(NULL, &info) == GRAPH_ERR_INVALID_ARGUMENT);
}

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

/* The same answers as the in-memory algorithms, on a graph with hubs and unconnected users. */
static void test_matches_in_memory(void)
{
    const size_t n = 60000;
    Graph *g = graph_gen_barabasi_albert(GRAPH_LIST, n / 2, 2, 5);
    CHECK(g != NULL);
    graph_set_adaptive(g, false);
    for (size_t i = 0; i < n / 2; i++) (void)graph_add_vertex(g);        /* half the users know nobody */
    for (int i = 0; i < 400; i++) {                                       /* a few small islands */
        const int a = (int)(n / 2 + (size_t)i * 7), b = a + 1;
        (void)graph_add_edge(g, a, b, 1.0);
    }
    CHECK(graph_id_bound(g) == n);

    FILE *f = tmpfile();
    CHECK(f != NULL && graph_write_edge_list(g, f) == GRAPH_OK);
    rewind(f);
    GraphStreamInfo info;
    CHECK(graph_stream_scan(f, &info) == GRAPH_OK && info.edges == graph_edge_count(g));
    CHECK(info.vertices <= n);

    rewind(f);
    size_t *deg = calloc(n, sizeof *deg);
    CHECK(graph_stream_degrees(f, n, false, deg, NULL) == GRAPH_OK);
    for (size_t v = 0; v < n; v++) CHECK(deg[v] == graph_degree(g, (int)v));

    rewind(f);
    int *streamed = malloc(n * sizeof *streamed), *direct = malloc(n * sizeof *direct);
    size_t sc = 0;
    CHECK(graph_stream_components(f, n, streamed, &sc) == GRAPH_OK);
    const size_t dc = graph_connected_components(g, direct);
    CHECK(sc == dc && same_partition(streamed, direct, n));
    free(deg); free(streamed); free(direct);
    fclose(f);
    graph_destroy(g);
}

/* A million users and five million friendships, written straight to a file (no graph is ever
 * built), then analysed with the address space capped at 200 MB.  Holding that graph as a list
 * takes ~230 MB, so this would not fit; streaming needs about 20 MB. */
static void test_million_users_capped(void)
{
    const size_t n = 1000000, m = 5000000;
    FILE *f = tmpfile();
    CHECK(f != NULL);
    uint64_t x = 88172645463325252ULL;
    size_t loops = 0;
    for (size_t i = 0; i < m; i++) {
        x ^= x << 13; x ^= x >> 7; x ^= x << 17;
        /* only the first 80% of users have friends; the rest are silent */
        const unsigned u = (unsigned)(x % (n * 8 / 10)), v = (unsigned)((x >> 24) % (n * 8 / 10));
        loops += u == v;
        fprintf(f, "%u %u\n", u, v);
    }
    fflush(f);
    rewind(f);

    struct rlimit lim = {(rlim_t)200 << 20, (rlim_t)200 << 20};
    CHECK(setrlimit(RLIMIT_AS, &lim) == 0);

    GraphStreamInfo info;
    CHECK(graph_stream_scan(f, &info) == GRAPH_OK && info.edges == m);
    CHECK(info.vertices > n * 8 / 10 - 1000 && info.vertices <= n * 8 / 10);

    rewind(f);
    size_t *deg = malloc(info.vertices * sizeof *deg);
    CHECK(deg != NULL && graph_stream_degrees(f, info.vertices, false, deg, NULL) == GRAPH_OK);
    size_t sum = 0, maxd = 0;
    for (size_t v = 0; v < info.vertices; v++) { sum += deg[v]; if (deg[v] > maxd) maxd = deg[v]; }
    CHECK(sum == 2 * m - loops);                          /* both ends of every edge; a self-loop counts once */
    CHECK(maxd > 0);
    free(deg);

    rewind(f);
    int *comp = malloc(info.vertices * sizeof *comp);
    size_t count = 0;
    CHECK(comp != NULL && graph_stream_components(f, info.vertices, comp, &count) == GRAPH_OK);
    CHECK(count >= 1 && count < info.vertices);

    /* a final pass proves the labelling is consistent: both ends of every edge share a label */
    rewind(f);
    unsigned a, b;
    size_t checked = 0;
    while (fscanf(f, "%u %u", &a, &b) == 2) {
        CHECK(comp[a] == comp[b]);
        checked++;
    }
    CHECK(checked == m);
    free(comp);
    fclose(f);
}

int main(void)
{
    test_small();
    test_matches_in_memory();
    test_million_users_capped();
    printf("All graph stream tests passed.\n");
    return 0;
}
