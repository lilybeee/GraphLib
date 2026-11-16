/*
 * test_list_analysis.c — traversal callbacks, centrality, clustering,
 * bipartite, eccentricity, transpose, similarity, statistics.
 */

#include "test_util_list.h"

/* ---- traversal callbacks ---- */

typedef struct {
    int    order[64];
    int    parent[64];
    size_t n;
    size_t stop_after;
} Trace;

static bool trace_fn(int v, int parent, void *ctx)
{
    Trace *t = ctx;
    t->order[t->n]  = v;
    t->parent[t->n] = parent;
    ++t->n;
    return t->stop_after == 0 || t->n < t->stop_after;
}

static void test_visit(void)
{
    for (int seed = 1; seed <= 20; ++seed) {
        for (int directed = 0; directed < 2; ++directed) {
            ListGraph *g = random_graph(20, 0.12, directed, 5, (unsigned)seed, 2);
            int start = 0;
            while (!list_graph_is_alive(g, start)) {
                ++start;
            }
            int order[64], parent[64];
            Trace t = { {0}, {0}, 0, 0 };

            size_t want = list_bfs(g, start, order, parent);
            CHECK(list_bfs_visit(g, start, trace_fn, &t) == want && t.n == want);
            CHECK(memcmp(order, t.order, want * sizeof(int)) == 0);
            for (size_t i = 0; i < want; ++i) {
                CHECK(t.parent[i] == parent[t.order[i]]);
            }

            t.n  = 0;
            want = list_dfs(g, start, order, parent);
            CHECK(list_dfs_visit(g, start, trace_fn, &t) == want && t.n == want);
            CHECK(memcmp(order, t.order, want * sizeof(int)) == 0);
            for (size_t i = 0; i < want; ++i) {
                CHECK(t.parent[i] == parent[t.order[i]]);
            }
            list_graph_destroy(g);
        }
    }
    /* Early stop, bad arguments. */
    ListGraph *p = list_gen_path(10, false);
    Trace t = { {0}, {0}, 0, 3 };
    CHECK(list_bfs_visit(p, 0, trace_fn, &t) == 3 && t.n == 3);
    t.n = 0;
    CHECK(list_dfs_visit(p, 0, trace_fn, &t) == 3 && t.n == 3);
    CHECK(list_bfs_visit(p, 99, trace_fn, &t) == 0);
    CHECK(list_dfs_visit(p, 0, NULL, &t) == 0);
    list_graph_destroy(p);
}

/* ---- centrality ---- */

static void test_degree_closeness_known(void)
{
    ListGraph *star = list_gen_star(5);
    double c[5];
    CHECK_OK(list_degree_centrality(star, c));
    CHECK(near(c[0], 1.0) && near(c[1], 0.25));

    ListGraph *path = list_gen_path(5, false);
    CHECK_OK(list_closeness_centrality(path, c));
    CHECK(near(c[2], 4.0 / 6.0));                       /* distances 1,1,2,2 */
    CHECK(near(c[0], 4.0 / 10.0));                      /* 1+2+3+4 = 10      */
    CHECK(c[2] > c[1] && c[1] > c[0]);

    ListGraph *isolated = list_graph_create(3, false);
    CHECK_OK(list_closeness_centrality(isolated, c));
    CHECK(c[0] == 0.0 && c[1] == 0.0);
    list_graph_destroy(isolated);
    list_graph_destroy(path);
    list_graph_destroy(star);
}

/** Brute-force betweenness: sigma via hop distances and path counting. */
static void brute_betweenness(const ListGraph *g, double *bc)
{
    const size_t n = list_graph_id_bound(g);
    int *D = malloc(n * n * sizeof(int));
    double *sigma = calloc(n * n, sizeof(double));
    hops(g, D);
    for (size_t s = 0; s < n; ++s) {
        if (!list_graph_is_alive(g, (int)s)) {
            continue;
        }
        /* Count shortest paths from s by processing vertices in distance order. */
        sigma[s * n + s] = 1.0;
        for (int d = 1; d < (int)n; ++d) {
            for (size_t v = 0; v < n; ++v) {
                if (D[s * n + v] != d) {
                    continue;
                }
                for (size_t u = 0; u < n; ++u) {
                    if (D[s * n + u] == d - 1 && list_graph_has_edge(g, (int)u, (int)v, NULL)) {
                        sigma[s * n + v] += sigma[s * n + u];
                    }
                }
            }
        }
    }
    for (size_t v = 0; v < n; ++v) {
        bc[v] = 0.0;
        for (size_t s = 0; s < n; ++s) {
            for (size_t t = 0; t < n; ++t) {
                if (s == v || t == v || s == t || D[s * n + t] >= HOP_INF || D[s * n + v] >= HOP_INF ||
                    D[v * n + t] >= HOP_INF) {
                    continue;
                }
                if (D[s * n + v] + D[v * n + t] == D[s * n + t]) {
                    bc[v] += sigma[s * n + v] * sigma[v * n + t] / sigma[s * n + t];
                }
            }
        }
        if (!list_graph_is_directed(g)) {
            bc[v] /= 2.0;
        }
    }
    free(D);
    free(sigma);
}

static void test_centrality_random(void)
{
    for (int seed = 1; seed <= 40; ++seed) {
        for (int directed = 0; directed < 2; ++directed) {
            ListGraph *g = random_graph(14, 0.15 + 0.02 * (seed % 5), directed, 5, (unsigned)seed,
                                          (size_t)(seed % 3));
            const size_t n = list_graph_id_bound(g);
            double got[14], want[14], close[14];
            CHECK_OK(list_betweenness_centrality(g, got));
            brute_betweenness(g, want);
            for (size_t v = 0; v < n; ++v) {
                if (list_graph_is_alive(g, (int)v)) {
                    CHECK(near(got[v], want[v]));
                } else {
                    CHECK(got[v] == 0.0);
                }
            }
            /* Closeness against hop distances. */
            int D[14 * 14];
            hops(g, D);
            CHECK_OK(list_closeness_centrality(g, close));
            const size_t vc = list_graph_vertex_count(g);
            for (size_t s = 0; s < n; ++s) {
                if (!list_graph_is_alive(g, (int)s)) {
                    CHECK(close[s] == 0.0);
                    continue;
                }
                size_t r = 0;
                double sum = 0;
                for (size_t t = 0; t < n; ++t) {
                    if (t != s && D[s * n + t] < HOP_INF) {
                        ++r;
                        sum += D[s * n + t];
                    }
                }
                const double expect = (r > 0 && vc > 1)
                    ? ((double)r / (double)(vc - 1)) * ((double)r / sum) : 0.0;
                CHECK(near(close[s], expect));
            }
            list_graph_destroy(g);
        }
    }
    ListGraph *path = list_gen_path(5, false);
    double b[5];
    CHECK_OK(list_betweenness_centrality(path, b));
    CHECK(near(b[0], 0) && near(b[1], 3) && near(b[2], 4) && near(b[3], 3) && near(b[4], 0));
    list_graph_destroy(path);
    ListGraph *star = list_gen_star(6);
    double bs[6];
    CHECK_OK(list_betweenness_centrality(star, bs));
    CHECK(near(bs[0], 10.0) && near(bs[1], 0.0));       /* C(5,2) pairs of leaves */
    list_graph_destroy(star);
}

/* ---- triangles and clustering ---- */

static void test_triangles_clustering(void)
{
    ListGraph *k4 = list_gen_complete(4, false);
    long long t;
    CHECK_OK(list_triangle_count(k4, &t));
    CHECK(t == 4);
    double local[4], avg;
    CHECK_OK(list_clustering_coefficient(k4, local, &avg));
    CHECK(near(avg, 1.0) && near(local[0], 1.0));
    list_graph_destroy(k4);

    ListGraph *tree = list_gen_random_tree(20, 5);
    CHECK_OK(list_triangle_count(tree, &t));
    CHECK(t == 0);
    CHECK_OK(list_clustering_coefficient(tree, NULL, &avg));
    CHECK(near(avg, 0.0));
    list_graph_destroy(tree);

    ListGraph *d = list_gen_complete(3, true);
    CHECK(list_triangle_count(d, &t) == GRAPH_ERR_UNSUPPORTED);
    CHECK(list_clustering_coefficient(d, local, &avg) == GRAPH_ERR_UNSUPPORTED);
    list_graph_destroy(d);

    for (int seed = 1; seed <= 40; ++seed) {
        ListGraph *g = random_graph(13, 0.3, false, 3, (unsigned)seed, (size_t)(seed % 3));
        const size_t n = list_graph_id_bound(g);
        long long want = 0;
        double lc[13], sum = 0.0;
        size_t alive = 0;
        for (size_t a = 0; a < n; ++a) {
            for (size_t b = a + 1; b < n; ++b) {
                for (size_t c = b + 1; c < n; ++c) {
                    want += list_graph_has_edge(g, (int)a, (int)b, NULL) &&
                            list_graph_has_edge(g, (int)b, (int)c, NULL) &&
                            list_graph_has_edge(g, (int)a, (int)c, NULL);
                }
            }
        }
        for (size_t v = 0; v < n; ++v) {
            double c = 0.0;
            if (list_graph_is_alive(g, (int)v)) {
                ++alive;
                size_t deg = 0, links = 0;
                for (size_t a = 0; a < n; ++a) {
                    if (a == v || !list_graph_has_edge(g, (int)v, (int)a, NULL)) {
                        continue;
                    }
                    ++deg;
                    for (size_t b = a + 1; b < n; ++b) {
                        links += b != v && list_graph_has_edge(g, (int)v, (int)b, NULL) &&
                                 list_graph_has_edge(g, (int)a, (int)b, NULL);
                    }
                }
                c = deg >= 2 ? 2.0 * (double)links / ((double)deg * (double)(deg - 1)) : 0.0;
                sum += c;
            }
            lc[v] = c;
        }
        double got_local[13], got_avg;
        CHECK_OK(list_triangle_count(g, &t));
        CHECK(t == want);
        CHECK_OK(list_clustering_coefficient(g, got_local, &got_avg));
        for (size_t v = 0; v < n; ++v) {
            CHECK(near(got_local[v], lc[v]));
        }
        CHECK(near(got_avg, alive ? sum / (double)alive : 0.0));
        list_graph_destroy(g);
    }
}

/* ---- bipartite ---- */

static void test_bipartite(void)
{
    bool r;
    int color[32];
    ListGraph *even = list_gen_cycle(6, false), *odd = list_gen_cycle(5, false);
    CHECK_OK(list_is_bipartite(even, &r, color));
    CHECK(r);
    for (int u = 0; u < 6; ++u) {
        CHECK(color[u] == 0 || color[u] == 1);
        CHECK(color[u] != color[(u + 1) % 6]);
    }
    CHECK_OK(list_is_bipartite(odd, &r, NULL));
    CHECK(!r);
    list_graph_destroy(even);
    list_graph_destroy(odd);

    ListGraph *grid = list_gen_grid(4, 5);
    CHECK_OK(list_is_bipartite(grid, &r, NULL));
    CHECK(r);
    list_graph_destroy(grid);

    ListGraph *loop = list_graph_create(3, false);
    CHECK(list_graph_add_edge(loop, 1, 1, 1.0));
    CHECK_OK(list_is_bipartite(loop, &r, NULL));
    CHECK(!r);
    list_graph_destroy(loop);

    /* A directed graph is judged by its underlying undirected graph. */
    ListGraph *dtri = list_gen_cycle(3, true), *dsq = list_gen_cycle(4, true);
    CHECK_OK(list_is_bipartite(dtri, &r, NULL));
    CHECK(!r);
    CHECK_OK(list_is_bipartite(dsq, &r, NULL));
    CHECK(r);
    list_graph_destroy(dtri);
    list_graph_destroy(dsq);

    for (int seed = 1; seed <= 30; ++seed) {      /* random graphs vs. 2^n colourings */
        ListGraph *g = random_graph(10, 0.2, seed % 2, 3, (unsigned)seed, (size_t)(seed % 2));
        const size_t n = list_graph_id_bound(g);
        bool exists = false;
        for (unsigned mask = 0; mask < (1u << n) && !exists; ++mask) {
            bool ok = true;
            for (size_t u = 0; u < n && ok; ++u) {
                for (size_t v = 0; v < n && ok; ++v) {
                    if (list_graph_is_alive(g, (int)u) && list_graph_is_alive(g, (int)v) &&
                        list_graph_has_edge(g, (int)u, (int)v, NULL)) {
                        ok = ((mask >> u) & 1u) != ((mask >> v) & 1u);
                    }
                }
            }
            exists = ok;
        }
        CHECK_OK(list_is_bipartite(g, &r, color));
        CHECK(r == exists);
        if (r) {
            for (size_t u = 0; u < n; ++u) {
                for (size_t v = 0; v < n; ++v) {
                    if (list_graph_is_alive(g, (int)u) && list_graph_has_edge(g, (int)u, (int)v, NULL)) {
                        CHECK(color[u] != color[v]);
                    }
                }
            }
        }
        list_graph_destroy(g);
    }
}

/* ---- eccentricity, diameter, radius ---- */

static void test_eccentricity(void)
{
    ListGraph *p = list_gen_path(6, false);
    int ecc[6], diam, rad;
    CHECK_OK(list_eccentricity(p, ecc));
    const int want[6] = { 5, 4, 3, 3, 4, 5 };
    CHECK(memcmp(ecc, want, sizeof(want)) == 0);
    CHECK_OK(list_diameter_radius(p, &diam, &rad));
    CHECK(diam == 5 && rad == 3);
    list_graph_destroy(p);

    ListGraph *s = list_gen_star(7);
    CHECK_OK(list_diameter_radius(s, &diam, &rad));
    CHECK(diam == 2 && rad == 1);
    list_graph_destroy(s);

    ListGraph *two = list_graph_create(4, false);
    CHECK(list_graph_add_edge(two, 0, 1, 1) && list_graph_add_edge(two, 2, 3, 1));
    CHECK(list_diameter_radius(two, &diam, &rad) == GRAPH_ERR_UNSUPPORTED);
    list_graph_destroy(two);

    ListGraph *empty = list_graph_create(0, false);
    CHECK(list_diameter_radius(empty, &diam, &rad) == GRAPH_ERR_INVALID_ARGUMENT);
    list_graph_destroy(empty);

    ListGraph *dc = list_gen_cycle(5, true);        /* strongly connected */
    CHECK_OK(list_diameter_radius(dc, &diam, &rad));
    CHECK(diam == 4 && rad == 4);
    list_graph_destroy(dc);
    ListGraph *dp = list_gen_path(4, true);         /* not strongly connected */
    CHECK(list_diameter_radius(dp, &diam, &rad) == GRAPH_ERR_UNSUPPORTED);
    list_graph_destroy(dp);

    for (int seed = 1; seed <= 30; ++seed) {            /* against hop distances */
        ListGraph *g = random_graph(12, 0.18, seed % 2, 3, (unsigned)seed, (size_t)(seed % 3));
        const size_t n = list_graph_id_bound(g);
        int D[144], e[12];
        hops(g, D);
        CHECK_OK(list_eccentricity(g, e));
        for (size_t u = 0; u < n; ++u) {
            if (!list_graph_is_alive(g, (int)u)) {
                CHECK(e[u] == -1);
                continue;
            }
            int m = 0;
            for (size_t v = 0; v < n; ++v) {
                if (D[u * n + v] < HOP_INF && D[u * n + v] > m) {
                    m = D[u * n + v];
                }
            }
            CHECK(e[u] == m);
        }
        list_graph_destroy(g);
    }
}

/* ---- transpose, neighbours, similarity, statistics ---- */

static void test_transpose(void)
{
    for (int seed = 1; seed <= 20; ++seed) {
        ListGraph *g = random_graph(12, 0.2, true, 7, (unsigned)seed, (size_t)(seed % 3));
        ListGraph *t = NULL, *tt = NULL;
        CHECK_OK(list_transpose(g, &t));
        const size_t n = list_graph_id_bound(g);
        CHECK(list_graph_id_bound(t) == n && list_graph_edge_count(t) == list_graph_edge_count(g));
        for (size_t u = 0; u < n; ++u) {
            CHECK(list_graph_is_alive(t, (int)u) == list_graph_is_alive(g, (int)u));
            for (size_t v = 0; v < n; ++v) {
                CHECK(wt(t, (int)v, (int)u) == wt(g, (int)u, (int)v));
                CHECK(list_graph_has_edge(t, (int)v, (int)u, NULL) ==
                      list_graph_has_edge(g, (int)u, (int)v, NULL));
            }
        }
        CHECK_OK(list_transpose(t, &tt));
        CHECK(same_graph(g, tt));
        list_graph_destroy(g);
        list_graph_destroy(t);
        list_graph_destroy(tt);
    }
    ListGraph *u = random_graph(9, 0.3, false, 4, 3, 1), *tu = NULL;
    CHECK_OK(list_transpose(u, &tu));
    CHECK(same_graph(u, tu));
    list_graph_destroy(u);
    list_graph_destroy(tu);
}

static void test_similarity_and_stats(void)
{
    ListGraph *k4 = list_gen_complete(4, false);
    int common[4];
    size_t cnt;
    double j;
    CHECK_OK(list_common_neighbors(k4, 0, 1, common, &cnt));
    CHECK(cnt == 2 && common[0] == 2 && common[1] == 3);
    CHECK_OK(list_jaccard(k4, 0, 1, &j));
    CHECK(near(j, 0.5));
    ListNetworkStats st;
    CHECK_OK(list_network_stats(k4, &st));
    CHECK(near(st.density, 1.0) && near(st.reciprocity, 1.0) && st.max_degree == 3);
    CHECK(list_common_neighbors(k4, 0, 9, common, &cnt) == GRAPH_ERR_INVALID_ARGUMENT);
    list_graph_destroy(k4);

    ListGraph *iso = list_graph_create(3, false);
    CHECK_OK(list_jaccard(iso, 0, 1, &j));
    CHECK(j == 0.0);
    list_graph_destroy(iso);

    ListGraph *c3 = list_gen_cycle(3, true);
    CHECK_OK(list_network_stats(c3, &st));
    CHECK(near(st.density, 0.5) && near(st.reciprocity, 0.0) && st.max_degree == 1);
    list_graph_destroy(c3);

    ListGraph *half = list_graph_create(3, true);
    CHECK(list_graph_add_edge(half, 0, 1, 1) && list_graph_add_edge(half, 1, 0, 1) &&
          list_graph_add_edge(half, 1, 2, 1));
    CHECK_OK(list_network_stats(half, &st));
    CHECK(near(st.reciprocity, 2.0 / 3.0));
    list_graph_destroy(half);

    ListGraph *star = list_gen_star(6);
    size_t hist[4];
    CHECK_OK(list_degree_histogram(star, hist, 4));
    CHECK(hist[0] == 0 && hist[1] == 5 && hist[2] == 0 && hist[3] == 1);  /* hub degree 5 clamps */
    size_t big[8];
    CHECK_OK(list_degree_histogram(star, big, 8));
    CHECK(big[1] == 5 && big[5] == 1);
    CHECK(list_degree_histogram(star, hist, 0) == GRAPH_ERR_INVALID_ARGUMENT);
    list_graph_destroy(star);

    /* Jaccard / common neighbours against a brute-force set computation. */
    for (int seed = 1; seed <= 20; ++seed) {
        ListGraph *g = random_graph(12, 0.3, seed % 2, 3, (unsigned)seed, 0);
        for (int u = 0; u < 12; ++u) {
            for (int v = 0; v < 12; ++v) {
                size_t inter = 0, uni = 0;
                for (int w = 0; w < 12; ++w) {
                    const bool a = list_graph_has_edge(g, u, w, NULL);
                    const bool b = list_graph_has_edge(g, v, w, NULL);
                    inter += a && b;
                    uni += a || b;
                }
                int out[12];
                size_t c;
                CHECK_OK(list_common_neighbors(g, u, v, out, &c));
                CHECK(c == inter);
                CHECK_OK(list_jaccard(g, u, v, &j));
                CHECK(near(j, uni ? (double)inter / (double)uni : 0.0));
            }
        }
        list_graph_destroy(g);
    }
}

/* ---- argument checks ---- */

static void test_bad_arguments(void)
{
    double d[4];
    CHECK(list_degree_centrality(NULL, d) == GRAPH_ERR_INVALID_ARGUMENT);
    CHECK(list_betweenness_centrality(NULL, d) == GRAPH_ERR_INVALID_ARGUMENT);
    CHECK(list_closeness_centrality(NULL, d) == GRAPH_ERR_INVALID_ARGUMENT);
    ListGraph *g = list_gen_path(3, false);
    CHECK(list_degree_centrality(g, NULL) == GRAPH_ERR_INVALID_ARGUMENT);
    ListNetworkStats st;
    CHECK(list_network_stats(NULL, &st) == GRAPH_ERR_INVALID_ARGUMENT);
    list_graph_destroy(g);
    CHECK(strcmp(graph_status_string(GRAPH_OK), "ok") == 0);
    CHECK(strlen(graph_status_string(GRAPH_ERR_NEGATIVE_CYCLE)) > 3);

    ListGraph *empty = list_graph_create(0, false);        /* nothing to do, nothing to break */
    double none[1];
    long long t;
    CHECK_OK(list_degree_centrality(empty, none));
    CHECK_OK(list_betweenness_centrality(empty, none));
    CHECK_OK(list_triangle_count(empty, &t));
    CHECK(t == 0);
    list_graph_destroy(empty);
}

int main(void)
{
    test_visit();
    test_degree_closeness_known();
    test_centrality_random();
    test_triangles_clustering();
    test_bipartite();
    test_eccentricity();
    test_transpose();
    test_similarity_and_stats();
    test_bad_arguments();
    printf("All list analysis tests passed.\n");
    return 0;
}
