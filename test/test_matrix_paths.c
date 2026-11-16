/*
 * test_matrix_paths.c — Bellman-Ford, Floyd-Warshall, Johnson, A*, multi-source
 * Dijkstra, 0-1 BFS and Yen's k shortest paths, all cross-checked against each
 * other and against brute-force path enumeration.
 */

#include "test_util_matrix.h"

#define N 12

static double path_weight(const MatrixGraph *g, const int *p, size_t len)
{
    double c = 0.0;
    for (size_t i = 0; i + 1 < len; ++i) {
        double w;
        CHECK(matrix_graph_has_edge(g, p[i], p[i + 1], &w));
        c += w;
    }
    return c;
}

/* ---- reconstruction helpers ---- */

static void test_reconstruction(void)
{
    const int parent[] = { -1, 0, 1, 1, -1 };
    int path[8];
    CHECK(matrix_path_from_parent(parent, 5, 0, 2, path, 8) == 3);
    CHECK(path[0] == 0 && path[1] == 1 && path[2] == 2);
    CHECK(matrix_path_from_parent(parent, 5, 0, 0, path, 8) == 1 && path[0] == 0);
    CHECK(matrix_path_from_parent(parent, 5, 0, 4, path, 8) == 0);        /* unreachable   */
    CHECK(matrix_path_from_parent(parent, 5, 0, 2, path, 2) == 0);        /* cap too small */
    CHECK(matrix_path_from_parent(parent, 5, 0, 9, path, 8) == 0);        /* out of range  */
    CHECK(matrix_path_from_parent(NULL, 5, 0, 2, path, 8) == 0);

    const int loopy[] = { 1, 0 };                                         /* a cycle must not hang */
    CHECK(matrix_path_from_parent(loopy, 2, 5 - 4, 0, path, 8) <= 2);

    /* next[u*3+v]: 0 -> 1 -> 2 */
    const int nx[9] = { -1, 1, 1, -1, -1, 2, -1, -1, -1 };
    CHECK(matrix_path_from_next(nx, 3, 0, 2, path, 8) == 3);
    CHECK(path[0] == 0 && path[1] == 1 && path[2] == 2);
    CHECK(matrix_path_from_next(nx, 3, 2, 0, path, 8) == 0);
    CHECK(matrix_path_from_next(nx, 3, 1, 1, path, 8) == 1);
}

/* ---- positive weights: every algorithm must agree ---- */

static void test_positive_weights(void)
{
    for (int seed = 1; seed <= 60; ++seed) {
        for (int directed = 0; directed < 2; ++directed) {
            MatrixGraph *g = random_graph(N, 0.2, directed, 9, (unsigned)seed, (size_t)(seed % 3));
            const size_t n = matrix_graph_id_bound(g);
            double fw[N * N], jo[N * N];
            int fnext[N * N], jnext[N * N];
            CHECK_OK(matrix_floyd_warshall(g, fw, fnext));
            CHECK_OK(matrix_johnson(g, jo, jnext));
            for (size_t i = 0; i < n * n; ++i) {
                CHECK(near(fw[i], jo[i]));
            }
            for (size_t s = 0; s < n; ++s) {
                if (!matrix_graph_is_alive(g, (int)s)) {
                    for (size_t v = 0; v < n; ++v) {
                        CHECK(isinf(fw[s * n + v]) && isinf(jo[s * n + v]));
                    }
                    continue;
                }
                double dj[N], bf[N];
                int pj[N], pb[N];
                CHECK(matrix_dijkstra(g, (int)s, dj, pj));
                CHECK_OK(matrix_bellman_ford(g, (int)s, bf, pb));
                for (size_t v = 0; v < n; ++v) {
                    CHECK(near(dj[v], fw[s * n + v]));
                    CHECK(near(bf[v], fw[s * n + v]));
                    if (isinf(dj[v])) {
                        continue;
                    }
                    int path[N];
                    size_t len = matrix_path_from_parent(pb, n, (int)s, (int)v, path, N);
                    CHECK(len >= 1 && near(path_weight(g, path, len), dj[v]));
                    len = matrix_path_from_next(fnext, n, (int)s, (int)v, path, N);
                    CHECK(len >= 1 && near(path_weight(g, path, len), dj[v]));
                    len = matrix_path_from_next(jnext, n, (int)s, (int)v, path, N);
                    CHECK(len >= 1 && near(path_weight(g, path, len), dj[v]));
                }
            }
            matrix_graph_destroy(g);
        }
    }
}

/* ---- negative weights and negative cycles ---- */

static void test_negative_weights(void)
{
    for (int seed = 1; seed <= 60; ++seed) {
        /* A DAG has no cycles at all, so negative weights cannot make a negative cycle. */
        MatrixGraph *g = matrix_gen_random_dag(N, 0.3, (unsigned)seed, -6.0, 6.0);
        double fw[N * N], jo[N * N], bf[N];
        CHECK_OK(matrix_floyd_warshall(g, fw, NULL));
        CHECK_OK(matrix_johnson(g, jo, NULL));
        for (size_t i = 0; i < N * N; ++i) {
            CHECK(near(fw[i], jo[i]));
        }
        for (int s = 0; s < N; ++s) {
            CHECK_OK(matrix_bellman_ford(g, s, bf, NULL));
            for (int v = 0; v < N; ++v) {
                CHECK(near(bf[v], fw[s * N + v]));
            }
        }
        matrix_graph_destroy(g);
    }

    const double cyc[] = { 0,1,1,  1,2,-3,  2,0,1,  2,3,5 };            /* cycle weight -1 */
    MatrixGraph *neg = build(5, true, cyc, 4);
    double d[5], all[25];
    CHECK(matrix_bellman_ford(neg, 0, d, NULL) == GRAPH_ERR_NEGATIVE_CYCLE);
    CHECK(matrix_floyd_warshall(neg, all, NULL) == GRAPH_ERR_NEGATIVE_CYCLE);
    CHECK(matrix_johnson(neg, all, NULL) == GRAPH_ERR_NEGATIVE_CYCLE);
    /* Bellman-Ford only reports cycles reachable from its source. */
    CHECK_OK(matrix_bellman_ford(neg, 4, d, NULL));
    CHECK(d[4] == 0.0 && isinf(d[0]));
    matrix_graph_destroy(neg);

    MatrixGraph *und = matrix_graph_create(3, false);             /* one negative undirected edge */
    CHECK(matrix_graph_add_edge(und, 0, 1, -1.0));
    CHECK(matrix_bellman_ford(und, 0, d, NULL) == GRAPH_ERR_NEGATIVE_CYCLE);
    CHECK(matrix_floyd_warshall(und, all, NULL) == GRAPH_ERR_NEGATIVE_CYCLE);
    CHECK(matrix_johnson(und, all, NULL) == GRAPH_ERR_NEGATIVE_CYCLE);
    matrix_graph_destroy(und);

    MatrixGraph *self = matrix_graph_create(2, true);
    CHECK(matrix_graph_add_edge(self, 0, 0, -2.0) && matrix_graph_add_edge(self, 0, 1, 1.0));
    CHECK(matrix_bellman_ford(self, 0, d, NULL) == GRAPH_ERR_NEGATIVE_CYCLE);
    CHECK(matrix_floyd_warshall(self, all, NULL) == GRAPH_ERR_NEGATIVE_CYCLE);
    matrix_graph_destroy(self);

    MatrixGraph *g = matrix_gen_path(3, true);
    CHECK(matrix_bellman_ford(g, 7, d, NULL) == GRAPH_ERR_INVALID_ARGUMENT);
    CHECK(matrix_bellman_ford(NULL, 0, d, NULL) == GRAPH_ERR_INVALID_ARGUMENT);
    matrix_graph_destroy(g);
}

/* ---- multi-source Dijkstra ---- */

static void test_multi_source(void)
{
    for (int seed = 1; seed <= 40; ++seed) {
        MatrixGraph *g = random_graph(N, 0.22, seed % 2, 9, (unsigned)seed, (size_t)(seed % 3));
        const size_t n = matrix_graph_id_bound(g);
        double fw[N * N], dist[N];
        int origin[N], parent[N], src[3], ns = 0;
        CHECK_OK(matrix_floyd_warshall(g, fw, NULL));
        for (int v = 0; v < (int)n && ns < 3; v += 3) {
            if (matrix_graph_is_alive(g, v)) {
                src[ns++] = v;
            }
        }
        if (ns == 0) {
            matrix_graph_destroy(g);
            continue;
        }
        CHECK_OK(matrix_dijkstra_multi(g, src, (size_t)ns, dist, parent, origin));
        for (size_t v = 0; v < n; ++v) {
            double best = INFINITY;
            for (int i = 0; i < ns; ++i) {
                best = fmin(best, fw[(size_t)src[i] * n + v]);
            }
            CHECK(near(dist[v], best));
            if (!isinf(best)) {
                CHECK(near(fw[(size_t)origin[v] * n + v], best));
                bool is_source = false;
                for (int i = 0; i < ns; ++i) {
                    is_source = is_source || origin[v] == src[i];
                }
                CHECK(is_source);
            } else {
                CHECK(origin[v] == -1 && parent[v] == -1);
            }
        }
        matrix_graph_destroy(g);
    }
    MatrixGraph *g = matrix_gen_path(4, false);
    double d[4];
    const int dead[] = { 9 }, ok[] = { 0 };
    CHECK(matrix_dijkstra_multi(g, dead, 1, d, NULL, NULL) == GRAPH_ERR_INVALID_ARGUMENT);
    CHECK(matrix_dijkstra_multi(g, ok, 0, d, NULL, NULL) == GRAPH_ERR_INVALID_ARGUMENT);
    CHECK(matrix_graph_add_edge(g, 0, 3, -1.0) == true);
    CHECK(matrix_dijkstra_multi(g, ok, 1, d, NULL, NULL) == GRAPH_ERR_INVALID_ARGUMENT);
    matrix_graph_destroy(g);
}

/* ---- 0-1 BFS ---- */

static void test_zero_one_bfs(void)
{
    for (int seed = 1; seed <= 50; ++seed) {
        t_seed((unsigned)seed + 1000u);
        const bool directed = seed % 2;
        MatrixGraph *g = matrix_graph_create(N, directed);
        for (int u = 0; u < N; ++u) {
            for (int v = directed ? 0 : u + 1; v < N; ++v) {
                if (u != v && t_rand(100) < 25) {
                    CHECK(matrix_graph_add_edge(g, u, v, (double)t_rand(2)));
                }
            }
        }
        (void)matrix_graph_remove_vertex(g, (int)t_rand(N));
        double fw[N * N], d[N];
        int parent[N];
        CHECK_OK(matrix_floyd_warshall(g, fw, NULL));
        for (int s = 0; s < N; ++s) {
            if (!matrix_graph_is_alive(g, s)) {
                CHECK(matrix_zero_one_bfs(g, s, d, NULL) == GRAPH_ERR_INVALID_ARGUMENT);
                continue;
            }
            CHECK_OK(matrix_zero_one_bfs(g, s, d, parent));
            for (int v = 0; v < N; ++v) {
                CHECK(near(d[v], fw[s * N + v]));
                if (!isinf(d[v])) {
                    int path[N];
                    const size_t len = matrix_path_from_parent(parent, N, s, v, path, N);
                    CHECK(len >= 1 && near(path_weight(g, path, len), d[v]));
                }
            }
        }
        matrix_graph_destroy(g);
    }
    MatrixGraph *g = matrix_gen_path(3, true);
    CHECK(matrix_graph_add_edge(g, 0, 2, 2.0));
    double d[3];
    CHECK(matrix_zero_one_bfs(g, 0, d, NULL) == GRAPH_ERR_INVALID_ARGUMENT);   /* weight 2 */
    matrix_graph_destroy(g);
}

/* ---- A* ---- */

typedef struct {
    const double *dist;      /* true distances, row-major N x N */
    size_t        n;
    unsigned      salt;
} HCtx;

/** Admissible but wildly inconsistent: a random fraction of the true distance. */
static double h_inconsistent(int v, int target, void *ctx)
{
    const HCtx *c = ctx;
    const double d = c->dist[(size_t)v * c->n + (size_t)target];
    if (isinf(d)) {
        return 0.0;
    }
    const unsigned x = ((unsigned)v * 2654435761u) ^ c->salt;
    return d * (0.2 + 0.8 * (double)(x % 1000u) / 1000.0);
}

static double h_zero(int v, int target, void *ctx)
{
    (void)v; (void)target; (void)ctx;
    return 0.0;
}

static void test_a_star(void)
{
    for (int seed = 1; seed <= 50; ++seed) {
        for (int directed = 0; directed < 2; ++directed) {
            MatrixGraph *g = random_graph(N, 0.2, directed, 9, (unsigned)seed, (size_t)(seed % 3));
            const size_t n = matrix_graph_id_bound(g);
            double fw[N * N];
            CHECK_OK(matrix_floyd_warshall(g, fw, NULL));
            HCtx ctx = { fw, n, (unsigned)seed * 7919u };
            for (int s = 0; s < (int)n; ++s) {
                for (int t = 0; t < (int)n; ++t) {
                    if (!matrix_graph_is_alive(g, s) || !matrix_graph_is_alive(g, t)) {
                        continue;
                    }
                    MatrixHeuristic hs[3] = { NULL, h_zero, h_inconsistent };
                    for (int k = 0; k < 3; ++k) {
                        MatrixPath p = { NULL, 0, 0.0 };
                        const MatrixStatus st = matrix_a_star(g, s, t, hs[k], &ctx, &p);
                        const double want = fw[(size_t)s * n + (size_t)t];
                        if (isinf(want)) {
                            CHECK(st == GRAPH_ERR_NOT_FOUND);
                            continue;
                        }
                        CHECK(st == GRAPH_OK);
                        CHECK(near(p.cost, want));
                        CHECK(p.length >= 1 && p.vertices[0] == s && p.vertices[p.length - 1] == t);
                        CHECK(near(path_weight(g, p.vertices, p.length), want));
                        matrix_path_free(&p);
                    }
                }
            }
            matrix_graph_destroy(g);
        }
    }
    MatrixGraph *g = matrix_gen_path(4, true);
    MatrixPath p = { NULL, 0, 0.0 };
    CHECK_OK(matrix_a_star(g, 2, 2, NULL, NULL, &p));
    CHECK(p.length == 1 && p.cost == 0.0);
    matrix_path_free(&p);
    CHECK(matrix_a_star(g, 3, 0, NULL, NULL, &p) == GRAPH_ERR_NOT_FOUND);
    CHECK(matrix_a_star(g, 9, 0, NULL, NULL, &p) == GRAPH_ERR_INVALID_ARGUMENT);
    CHECK(matrix_graph_add_edge(g, 0, 3, -2.0));
    CHECK(matrix_a_star(g, 0, 3, NULL, NULL, &p) == GRAPH_ERR_INVALID_ARGUMENT);
    matrix_graph_destroy(g);
}

/* ---- Yen's k shortest paths ---- */

typedef struct {
    const MatrixGraph *g;
    int    target;
    bool   on_path[16];
    int    stack[16];
    double costs[4096];
    size_t count;
} Enum;

static void enumerate(Enum *e, int u, size_t depth, double cost)
{
    if (u == e->target) {
        CHECK(e->count < 4096);
        e->costs[e->count++] = cost;
        return;
    }
    for (int v = 0; v < 16; ++v) {
        double w;
        if (v < (int)matrix_graph_id_bound(e->g) && !e->on_path[v] &&
            matrix_graph_has_edge(e->g, u, v, &w)) {
            e->on_path[v] = true;
            e->stack[depth] = v;
            enumerate(e, v, depth + 1, cost + w);
            e->on_path[v] = false;
        }
    }
}

static int cmp_double(const void *a, const void *b)
{
    const double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

static void test_k_shortest(void)
{
    for (int seed = 1; seed <= 80; ++seed) {
        for (int directed = 0; directed < 2; ++directed) {
            MatrixGraph *g = random_graph(8, 0.4, directed, 4, (unsigned)seed, (size_t)(seed % 2));
            for (int s = 0; s < 8; s += 3) {
                for (int t = 1; t < 8; t += 3) {
                    if (!matrix_graph_is_alive(g, s) || !matrix_graph_is_alive(g, t)) {
                        continue;
                    }
                    Enum e;
                    memset(&e, 0, sizeof(e));
                    e.g = g;
                    e.target = t;
                    e.on_path[s] = true;
                    enumerate(&e, s, 0, 0.0);
                    qsort(e.costs, e.count, sizeof(double), cmp_double);

                    MatrixPath paths[7];
                    size_t found = 99;
                    const MatrixStatus st = matrix_k_shortest_paths(g, s, t, 7, paths, &found);
                    if (e.count == 0) {
                        CHECK(st == GRAPH_ERR_NOT_FOUND && found == 0);
                        continue;
                    }
                    CHECK(st == GRAPH_OK);
                    CHECK(found == (e.count < 7 ? e.count : 7));
                    for (size_t i = 0; i < found; ++i) {
                        CHECK(near(paths[i].cost, e.costs[i]));
                        CHECK(near(path_weight(g, paths[i].vertices, paths[i].length), paths[i].cost));
                        CHECK(paths[i].vertices[0] == s && paths[i].vertices[paths[i].length - 1] == t);
                        for (size_t a = 0; a < paths[i].length; ++a) {         /* loopless */
                            for (size_t b = a + 1; b < paths[i].length; ++b) {
                                CHECK(paths[i].vertices[a] != paths[i].vertices[b]);
                            }
                        }
                        for (size_t j = 0; j < i; ++j) {                       /* distinct */
                            CHECK(paths[i].length != paths[j].length ||
                                  memcmp(paths[i].vertices, paths[j].vertices,
                                         paths[i].length * sizeof(int)) != 0);
                        }
                    }
                    for (size_t i = 0; i < found; ++i) {
                        matrix_path_free(&paths[i]);
                    }
                }
            }
            matrix_graph_destroy(g);
        }
    }
    MatrixGraph *g = matrix_gen_path(3, true);
    MatrixPath paths[2];
    size_t found;
    CHECK(matrix_k_shortest_paths(g, 0, 2, 0, paths, &found) == GRAPH_ERR_INVALID_ARGUMENT);
    CHECK_OK(matrix_k_shortest_paths(g, 1, 1, 2, paths, &found));
    CHECK(found == 1 && paths[0].length == 1);
    matrix_path_free(&paths[0]);
    CHECK(matrix_k_shortest_paths(g, 2, 0, 2, paths, &found) == GRAPH_ERR_NOT_FOUND);
    matrix_graph_destroy(g);
}

/* ---- all-pairs size limit ---- */

static void test_too_large(void)
{
    MatrixGraph *g = matrix_graph_create(8193, false);
    double d[1];
    CHECK(matrix_floyd_warshall(g, d, NULL) == GRAPH_ERR_TOO_LARGE);
    CHECK(matrix_johnson(g, d, NULL) == GRAPH_ERR_TOO_LARGE);
    matrix_graph_destroy(g);
}

int main(void)
{
    test_reconstruction();
    test_positive_weights();
    test_negative_weights();
    test_multi_source();
    test_zero_one_bfs();
    test_a_star();
    test_k_shortest();
    test_too_large();
    printf("All matrix path tests passed.\n");
    return 0;
}
