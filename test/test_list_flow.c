/*
 * test_list_flow.c — Prim, max flow / min cut, bipartite matching, Hungarian.
 * Every result is checked against brute force on small random instances.
 */

#include "test_util_list.h"

/* ---- Prim ---- */

static void test_prim(void)
{
    for (int seed = 1; seed <= 80; ++seed) {
        ListGraph *g = random_graph(14, 0.12 + 0.02 * (seed % 8), false, 20, (unsigned)seed, (size_t)(seed % 3));
        ListEdge pe[16], ke[16];
        size_t pc = 0, kc = 0;
        const double prim = list_mst_prim(g, pe, &pc);
        const double krus = list_mst_kruskal(g, ke, &kc);
        if (isinf(krus)) {
            CHECK(isinf(prim) && pc == 0);
        } else {
            CHECK(near(prim, krus) && pc == kc);
            CHECK(pc + 1 == list_graph_vertex_count(g) || list_graph_vertex_count(g) <= 1);
            /* The edges are real, sum to the total, and join everything without a cycle. */
            double sum = 0.0;
            ListGraph *tree = list_graph_create(list_graph_id_bound(g), false);
            for (size_t v = 0; v < list_graph_id_bound(g); ++v) {
                if (!list_graph_is_alive(g, (int)v)) {
                    CHECK(list_graph_remove_vertex(tree, (int)v));
                }
            }
            for (size_t i = 0; i < pc; ++i) {
                CHECK(pe[i].weight == wt(g, pe[i].u, pe[i].v) && pe[i].u != pe[i].v);
                CHECK(list_graph_add_edge(tree, pe[i].u, pe[i].v, pe[i].weight));
                sum += pe[i].weight;
            }
            CHECK(near(sum, prim));
            CHECK(list_is_connected(tree) && list_is_forest(tree));
            list_graph_destroy(tree);
        }
        list_graph_destroy(g);
    }
    ListGraph *d = list_gen_path(3, true), *one = list_graph_create(1, false);
    ListEdge e[3];
    size_t c = 99;
    CHECK(isinf(list_mst_prim(d, e, &c)));                       /* directed */
    CHECK(list_mst_prim(one, e, &c) == 0.0 && c == 0);           /* single vertex */
    CHECK(isinf(list_mst_prim(NULL, e, &c)));
    list_graph_destroy(d);
    list_graph_destroy(one);
}

/* ---- max flow / min cut ---- */

/** Capacity of the cut (S, V \ S) for the subset mask S. */
static double cut_capacity(const ListGraph *g, unsigned mask)
{
    const size_t n = list_graph_id_bound(g);
    const bool directed = list_graph_is_directed(g);
    double c = 0.0;
    for (size_t u = 0; u < n; ++u) {
        for (size_t v = 0; v < n; ++v) {
            if (u == v || !list_graph_has_edge(g, (int)u, (int)v, NULL)) {
                continue;
            }
            const bool in_u = (mask >> u) & 1u, in_v = (mask >> v) & 1u;
            if (in_u && !in_v) {
                c += wt(g, (int)u, (int)v);
            } else if (!directed && !in_u && in_v) {
                /* Undirected edge {u,v} is seen from both ends; count it once above. */
            }
        }
    }
    return c;
}

static double brute_min_cut(const ListGraph *g, int s, int t)
{
    const size_t n = list_graph_id_bound(g);
    double best = INFINITY;
    for (unsigned mask = 0; mask < (1u << n); ++mask) {
        if (!((mask >> s) & 1u) || ((mask >> t) & 1u)) {
            continue;
        }
        const double c = cut_capacity(g, mask);
        if (c < best) {
            best = c;
        }
    }
    return best;
}

static void test_max_flow(void)
{
    /* The textbook network: maximum flow 23. */
    const double clrs[] = { 0,1,16, 0,2,13, 1,2,10, 2,1,4, 1,3,12, 3,2,9,
                            2,4,14, 4,3,7, 3,5,20, 4,5,4 };
    ListGraph *g = build(6, true, clrs, 10);
    double value;
    bool side[6];
    CHECK_OK(list_max_flow(g, 0, 5, LIST_FLOW_DINIC, &value, side));
    CHECK(near(value, 23.0) && side[0] && !side[5]);
    CHECK_OK(list_max_flow(g, 0, 5, LIST_FLOW_EDMONDS_KARP, &value, NULL));
    CHECK(near(value, 23.0));
    ListEdge cut[10];
    size_t nc;
    CHECK_OK(list_min_cut(g, 0, 5, &value, cut, &nc));
    double sum = 0;
    for (size_t i = 0; i < nc; ++i) {
        sum += cut[i].weight;
    }
    CHECK(near(value, 23.0) && near(sum, 23.0));
    list_graph_destroy(g);

    for (int seed = 1; seed <= 120; ++seed) {
        for (int directed = 0; directed < 2; ++directed) {
            ListGraph *r = random_graph(9, 0.3, directed, 9, (unsigned)seed, 0);
            int s = 0, t = 8;
            double dinic = -1, ek = -1;
            bool sd[9];
            CHECK_OK(list_max_flow(r, s, t, LIST_FLOW_DINIC, &dinic, sd));
            CHECK_OK(list_max_flow(r, s, t, LIST_FLOW_EDMONDS_KARP, &ek, NULL));
            const double want = brute_min_cut(r, s, t);
            CHECK(near(dinic, want) && near(ek, want));
            /* source_side is itself a minimum cut. */
            unsigned mask = 0;
            for (int v = 0; v < 9; ++v) {
                mask |= (unsigned)sd[v] << v;
            }
            CHECK(sd[s] && !sd[t] && near(cut_capacity(r, mask), want));
            ListEdge ce[100];
            size_t cn;
            double mv;
            CHECK_OK(list_min_cut(r, s, t, &mv, ce, &cn));
            double total = 0;
            for (size_t i = 0; i < cn; ++i) {
                CHECK(sd[ce[i].u] && !sd[ce[i].v]);
                total += ce[i].weight;
            }
            CHECK(near(mv, want) && near(total, want));
            list_graph_destroy(r);
        }
    }
    ListGraph *p = list_gen_path(3, true);
    CHECK(list_max_flow(p, 0, 0, LIST_FLOW_DINIC, &value, NULL) == GRAPH_ERR_INVALID_ARGUMENT);
    CHECK(list_max_flow(p, 0, 9, LIST_FLOW_DINIC, &value, NULL) == GRAPH_ERR_INVALID_ARGUMENT);
    CHECK(list_graph_add_edge(p, 0, 2, -1.0));
    CHECK(list_max_flow(p, 0, 2, LIST_FLOW_DINIC, &value, NULL) == GRAPH_ERR_INVALID_ARGUMENT);
    list_graph_destroy(p);

    ListGraph *disc = list_graph_create(4, true);        /* no path: zero flow */
    CHECK(list_graph_add_edge(disc, 0, 1, 5.0));
    CHECK_OK(list_max_flow(disc, 0, 3, LIST_FLOW_DINIC, &value, side));
    CHECK(value == 0.0 && side[0] && side[1] && !side[3]);
    list_graph_destroy(disc);
}

/* ---- bipartite matching ---- */

static int brute_matching(const ListGraph *g, const int *left, size_t nl, size_t idx, unsigned used)
{
    if (idx == nl) {
        return 0;
    }
    int best = brute_matching(g, left, nl, idx + 1, used);          /* leave this vertex free */
    for (int v = 0; v < (int)list_graph_id_bound(g); ++v) {
        if (!((used >> v) & 1u) && list_graph_is_alive(g, v) &&
            (list_graph_has_edge(g, left[idx], v, NULL) || list_graph_has_edge(g, v, left[idx], NULL))) {
            const int c = 1 + brute_matching(g, left, nl, idx + 1, used | (1u << v));
            if (c > best) {
                best = c;
            }
        }
    }
    return best;
}

static void test_matching(void)
{
    for (int seed = 1; seed <= 100; ++seed) {
        t_seed((unsigned)seed + 500u);
        const int nl = 4 + (int)t_rand(4), nr = 4 + (int)t_rand(4);
        const bool directed = seed % 2;
        ListGraph *g = list_graph_create((size_t)(nl + nr), directed);
        for (int u = 0; u < nl; ++u) {
            for (int v = nl; v < nl + nr; ++v) {
                if (t_rand(100) < 35) {
                    const bool flip = directed && t_rand(2);          /* arcs either way round */
                    CHECK(flip ? list_graph_add_edge(g, v, u, 1.0) : list_graph_add_edge(g, u, v, 1.0));
                }
            }
        }
        if (seed % 5 == 0) {
            (void)list_graph_remove_vertex(g, (int)t_rand((unsigned)(nl + nr)));
        }
        int match[16], left[16];
        size_t size = 99, k = 0;
        CHECK_OK(list_bipartite_matching(g, match, &size));
        for (int u = 0; u < nl; ++u) {
            if (list_graph_is_alive(g, u)) {
                left[k++] = u;
            }
        }
        CHECK((int)size == brute_matching(g, left, k, 0, 0));
        size_t pairs = 0;
        for (int v = 0; v < nl + nr; ++v) {
            if (match[v] < 0) {
                continue;
            }
            CHECK(match[match[v]] == v && match[v] != v);
            CHECK(list_graph_has_edge(g, v, match[v], NULL) || list_graph_has_edge(g, match[v], v, NULL));
            ++pairs;
        }
        CHECK(pairs == 2 * size);
        list_graph_destroy(g);
    }
    ListGraph *tri = list_gen_cycle(3, false);
    int m[3];
    size_t sz;
    CHECK(list_bipartite_matching(tri, m, &sz) == GRAPH_ERR_UNSUPPORTED);
    list_graph_destroy(tri);

    ListGraph *k33 = list_graph_create(6, false);                /* perfect matching exists */
    for (int u = 0; u < 3; ++u) {
        for (int v = 3; v < 6; ++v) {
            CHECK(list_graph_add_edge(k33, u, v, 1.0));
        }
    }
    int m6[6];
    CHECK_OK(list_bipartite_matching(k33, m6, &sz));
    CHECK(sz == 3);
    list_graph_destroy(k33);

    ListGraph *none = list_graph_create(4, false);
    CHECK_OK(list_bipartite_matching(none, m6, &sz));
    CHECK(sz == 0);
    list_graph_destroy(none);
}

/* ---- Hungarian ---- */

static double brute_assign(const double *cost, size_t rows, size_t cols, size_t r, unsigned used)
{
    if (r == rows) {
        return 0.0;
    }
    double best = INFINITY;
    for (size_t c = 0; c < cols; ++c) {
        if (!((used >> c) & 1u)) {
            const double t = cost[r * cols + c] + brute_assign(cost, rows, cols, r + 1, used | (1u << c));
            if (t < best) {
                best = t;
            }
        }
    }
    return best;
}

static void test_hungarian(void)
{
    const double known[] = { 4, 1, 3,   2, 0, 5,   3, 2, 2 };       /* optimum 5: (0,1)(1,0)(2,2) */
    int a[3];
    double total;
    CHECK_OK(list_hungarian(known, 3, 3, a, &total));
    CHECK(near(total, 5.0) && a[0] == 1 && a[1] == 0 && a[2] == 2);

    for (int seed = 1; seed <= 150; ++seed) {
        t_seed((unsigned)seed + 77u);
        const size_t rows = 1 + t_rand(6), cols = rows + t_rand(3);
        double cost[64];
        for (size_t i = 0; i < rows * cols; ++i) {
            cost[i] = (double)((int)t_rand(41) - 15);                  /* negatives allowed */
        }
        int assign[8];
        double got;
        CHECK_OK(list_hungarian(cost, rows, cols, assign, &got));
        CHECK(near(got, brute_assign(cost, rows, cols, 0, 0)));
        double sum = 0;
        unsigned used = 0;
        for (size_t r = 0; r < rows; ++r) {
            CHECK(assign[r] >= 0 && (size_t)assign[r] < cols && !((used >> assign[r]) & 1u));
            used |= 1u << assign[r];
            sum += cost[r * cols + (size_t)assign[r]];
        }
        CHECK(near(sum, got));
    }
    CHECK(list_hungarian(known, 4, 3, a, &total) == GRAPH_ERR_INVALID_ARGUMENT);   /* rows > cols */
    CHECK(list_hungarian(NULL, 3, 3, a, &total) == GRAPH_ERR_INVALID_ARGUMENT);
    double bad[] = { 1, NAN, 2, 3 };
    CHECK(list_hungarian(bad, 2, 2, a, &total) == GRAPH_ERR_INVALID_ARGUMENT);
    double single[] = { 7 };
    CHECK_OK(list_hungarian(single, 1, 1, a, &total));
    CHECK(total == 7.0 && a[0] == 0);
}

int main(void)
{
    test_prim();
    test_max_flow();
    test_matching();
    test_hungarian();
    printf("All list flow tests passed.\n");
    return 0;
}
