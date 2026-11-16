/*
 * test_graph_social.c — the per-user social queries (graph_social.h), checked
 * against brute-force answers on many random graphs in both representations.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "graph.h"
#include "graph_algorithms.h"
#include "graph_social.h"

#define CHECK(cond)                                                                  \
    do {                                                                             \
        if (!(cond)) {                                                               \
            fprintf(stderr, "%s:%d: CHECK failed: %s\n", __FILE__, __LINE__, #cond); \
            exit(1);                                                                 \
        }                                                                            \
    } while (0)

static unsigned long long rs = 88172645463325252ULL;
static unsigned rnd(unsigned b) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return (unsigned)(rs % b); }

static const GraphRep REPS[2] = {GRAPH_LIST, GRAPH_MATRIX};

static Graph *make(GraphRep rep, size_t n, bool directed, const int (*e)[2], size_t m)
{
    Graph *g = graph_create_sized(n, directed);
    CHECK(g != NULL);
    graph_set_adaptive(g, false);
    for (size_t i = 0; i < m; i++) CHECK(graph_add_edge(g, e[i][0], e[i][1], 1.0));
    CHECK(graph_convert(g, rep));
    return g;
}

static void test_known(void)
{
    /* 0-1-2-3-4 is a path; 5 knows nobody; 0-2 shortcut */
    const int e[][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {0, 2}};
    for (int r = 0; r < 2; r++) {
        Graph *g = make(REPS[r], 6, false, e, 5);
        int users[8], hops[8];
        size_t count = 99;
        bool trunc = true;
        CHECK(graph_social_within_hops(g, 0, 2, 8, users, hops, &count, &trunc) == GRAPH_OK);
        CHECK(count == 3 && !trunc);                      /* 1, 2 at one hop; 3 at two */
        for (size_t i = 0; i < count; i++) {
            CHECK(hops[i] == (users[i] == 3 ? 2 : 1));
            CHECK(users[i] != 0);                         /* the source is not in its own results */
        }
        CHECK(hops[0] <= hops[1] && hops[1] <= hops[2]);  /* nearest first */
        CHECK(graph_social_within_hops(g, 0, 10, 2, users, hops, &count, &trunc) == GRAPH_OK);
        CHECK(count == 2 && trunc);                       /* capped, and says so */
        CHECK(graph_social_within_hops(g, 5, 3, 8, users, hops, &count, &trunc) == GRAPH_OK);
        CHECK(count == 0 && !trunc);                      /* a user with no connections */
        CHECK(graph_social_within_hops(g, 0, 0, 8, users, hops, &count, NULL) == GRAPH_OK && count == 0);
        CHECK(graph_social_within_hops(g, 99, 2, 8, users, hops, &count, NULL) == GRAPH_ERR_INVALID_VERTEX);
        CHECK(graph_social_within_hops(NULL, 0, 2, 8, users, hops, &count, NULL) == GRAPH_ERR_INVALID_ARGUMENT);

        GraphPath p = {0};
        CHECK(graph_social_path(g, 0, 4, SIZE_MAX, &p) == GRAPH_OK);
        CHECK(p.length == 4 && p.cost == 3.0 && p.vertices[0] == 0 && p.vertices[1] == 2 &&
              p.vertices[2] == 3 && p.vertices[3] == 4);   /* 0-2-3-4 uses the shortcut */
        graph_path_free(&p);
        CHECK(graph_social_path(g, 0, 4, 2, &p) == GRAPH_ERR_NOT_FOUND);      /* too far for 2 hops */
        CHECK(graph_social_path(g, 0, 5, SIZE_MAX, &p) == GRAPH_ERR_NOT_FOUND); /* 5 knows nobody */
        CHECK(graph_social_path(g, 3, 3, SIZE_MAX, &p) == GRAPH_OK && p.length == 1 && p.cost == 0);
        graph_path_free(&p);
        graph_destroy(g);

        /* directed graphs have no in-neighbours to search backwards through */
        Graph *d = make(REPS[r], 3, true, e, 2);
        CHECK(graph_social_path(d, 0, 2, SIZE_MAX, &p) == GRAPH_ERR_UNSUPPORTED);
        graph_destroy(d);
    }

    /* suggestions: 3 shares friends 1 and 2 with user 0, 4 shares only 1 */
    const int f[][2] = {{0, 1}, {0, 2}, {1, 3}, {2, 3}, {1, 4}, {6, 7}};
    for (int r = 0; r < 2; r++) {
        Graph *g = make(REPS[r], 8, false, f, 6);
        GraphSuggestion s[4];
        size_t count = 0;
        CHECK(graph_social_suggest_friends(g, 0, 4, s, &count) == GRAPH_OK);
        CHECK(count == 2 && s[0].user == 3 && s[0].mutual == 2 && s[1].user == 4 && s[1].mutual == 1);
        CHECK(graph_social_suggest_friends(g, 0, 1, s, &count) == GRAPH_OK && count == 1 && s[0].user == 3);
        CHECK(graph_social_suggest_friends(g, 5, 4, s, &count) == GRAPH_OK && count == 0);   /* no friends */
        CHECK(graph_social_suggest_friends(g, 99, 4, s, &count) == GRAPH_ERR_INVALID_VERTEX);
        graph_destroy(g);
    }
}

/* ---- random cross-check against brute force ---- */

static void bfs_dist(Graph *g, int src, int *dist)
{
    const size_t n = graph_id_bound(g);
    int *queue = malloc(n * sizeof *queue), *nb = malloc(n * sizeof *nb);
    for (size_t v = 0; v < n; v++) dist[v] = -1;
    size_t head = 0, tail = 0;
    dist[src] = 0;
    queue[tail++] = src;
    while (head < tail) {
        const int u = queue[head++];
        const size_t d = graph_neighbors(g, u, nb, NULL, n);
        for (size_t i = 0; i < d; i++) if (dist[nb[i]] < 0) { dist[nb[i]] = dist[u] + 1; queue[tail++] = nb[i]; }
    }
    free(queue); free(nb);
}

static void test_random(void)
{
    for (int trial = 0; trial < 150; trial++) {
        const size_t n = 20 + rnd(60);
        const bool directed = trial % 5 == 4;
        Graph *g = graph_create_sized(n, directed);
        graph_set_adaptive(g, false);
        for (size_t i = 0; i < n * 2; i++) {            /* sparse, so some users end up with nobody */
            const unsigned a = rnd((unsigned)n), b = rnd((unsigned)n);
            if (a != b) (void)graph_add_edge(g, (int)a, (int)b, 1.0);
        }
        CHECK(graph_convert(g, REPS[trial % 2]));
        int *dist = malloc(n * sizeof *dist), *users = malloc(n * sizeof *users), *hops = malloc(n * sizeof *hops);
        int *nbrs = malloc(n * sizeof *nbrs), *other = malloc(n * sizeof *other);

        for (int q = 0; q < 10; q++) {
            const int u = (int)rnd((unsigned)n);
            bfs_dist(g, u, dist);

            /* within k hops == the BFS ball of radius k, nearest first */
            const size_t k = 1 + rnd(4);
            size_t count = 0, want = 0;
            CHECK(graph_social_within_hops(g, u, k, n, users, hops, &count, NULL) == GRAPH_OK);
            for (size_t v = 0; v < n; v++) if ((int)v != u && dist[v] >= 1 && (size_t)dist[v] <= k) want++;
            CHECK(count == want);
            for (size_t i = 0; i < count; i++) {
                CHECK(hops[i] == dist[users[i]]);
                CHECK(i == 0 || hops[i - 1] <= hops[i]);
            }

            /* suggestions == brute-force mutual-friend counts */
            GraphSuggestion sug[8];
            size_t sc = 0;
            CHECK(graph_social_suggest_friends(g, u, 8, sug, &sc) == GRAPH_OK);
            const size_t ud = graph_neighbors(g, u, nbrs, NULL, n);
            size_t found = 0;
            int best_mutual = 0;
            for (size_t v = 0; v < n; v++) {
                if ((int)v == u || graph_has_edge(g, u, (int)v, NULL)) continue;
                int mutual = 0;
                for (size_t i = 0; i < ud; i++) if (graph_has_edge(g, nbrs[i], (int)v, NULL)) mutual++;
                if (mutual > 0) {
                    found++;
                    if (mutual > best_mutual) best_mutual = mutual;
                }
            }
            CHECK(sc == (found < 8 ? found : 8));
            for (size_t i = 0; i < sc; i++) {
                int mutual = 0;
                for (size_t j = 0; j < ud; j++) if (graph_has_edge(g, nbrs[j], sug[i].user, NULL)) mutual++;
                CHECK((size_t)mutual == sug[i].mutual);
                CHECK(sug[i].user != u && !graph_has_edge(g, u, sug[i].user, NULL));
                CHECK(i == 0 || sug[i - 1].mutual >= sug[i].mutual);
            }
            if (sc > 0) CHECK((int)sug[0].mutual == best_mutual);

            /* shortest path == BFS distance (undirected only) */
            if (!directed) {
                const int t = (int)rnd((unsigned)n);
                GraphPath p = {0};
                const GraphStatus st = graph_social_path(g, u, t, SIZE_MAX, &p);
                if (dist[t] < 0) {
                    CHECK(st == GRAPH_ERR_NOT_FOUND);
                } else {
                    CHECK(st == GRAPH_OK && p.length == (size_t)dist[t] + 1 && p.vertices[0] == u &&
                          p.vertices[p.length - 1] == t);
                    for (size_t i = 0; i + 1 < p.length; i++) CHECK(graph_has_edge(g, p.vertices[i], p.vertices[i + 1], NULL));
                    graph_path_free(&p);
                    if (dist[t] > 1) {                          /* one hop short is refused */
                        CHECK(graph_social_path(g, u, t, (size_t)dist[t] - 1, &p) == GRAPH_ERR_NOT_FOUND);
                    }
                }
            }
        }
        free(dist); free(users); free(hops); free(nbrs); free(other);
        graph_destroy(g);
    }
}

int main(void)
{
    test_known();
    test_random();
    printf("All graph social-query tests passed.\n");
    return 0;
}
