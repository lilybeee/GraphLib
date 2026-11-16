/*
 * test_list_advanced.c — known-answer and cross-check tests for list_advanced.
 *
 * Where a result has a closed form (cycle counts, traversal orders, small
 * schedules) it is checked against that.  Where the library already has an
 * independent algorithm (Dijkstra, Floyd–Warshall, A*) the advanced version is
 * compared with it on random graphs, so the two implementations do not share
 * bugs.
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "list_advanced.h"
#include "test_util_list.h"

#define N_(a) (sizeof(a) / sizeof((a)[0]))
#define BUILD(n, dir, arr) build((n), (dir), (arr), N_(arr) / 3)

static double zero_h(int v, int goal, void *ctx)
{
    (void)v; (void)goal; (void)ctx;
    return 0.0;
}

/* Checks that path is a real route in g from s to t and that its cost matches. */
static void check_path(const ListGraph *g, const ListAdvPath *p, int s, int t)
{
    CHECK(p->length >= 1 && p->vertices[0] == s && p->vertices[p->length - 1] == t);
    double cost = 0.0;
    for (size_t i = 0; i + 1 < p->length; ++i) {
        double w = 0.0;
        CHECK(list_graph_has_edge(g, p->vertices[i], p->vertices[i + 1], &w));
        cost += w;
    }
    CHECK(near(cost, p->cost));
}

static void test_longest_path(void)
{
    const double e[] = {0, 1, 2, 0, 2, 1, 1, 3, 3, 2, 3, 5};
    ListGraph *g = BUILD(5, true, e);
    double dist[5];
    int parent[5];
    CHECK_OK(list_adv_longest_path(g, 0, dist, parent));
    CHECK(near(dist[3], 6.0) && parent[3] == 2 && near(dist[0], 0.0));
    CHECK(isinf(dist[4]) && dist[4] < 0 && parent[4] == -1);
    CHECK(list_adv_longest_path(g, 9, dist, parent) == GRAPH_ERR_INVALID_VERTEX);
    CHECK(list_graph_add_edge(g, 3, 0, 1.0));
    CHECK(list_adv_longest_path(g, 0, dist, parent) == GRAPH_ERR_CYCLE);
    list_graph_destroy(g);

    ListGraph *u = list_graph_create(3, false);
    CHECK(list_adv_longest_path(u, 0, dist, parent) == GRAPH_ERR_UNSUPPORTED);
    list_graph_destroy(u);
}

static void test_pagerank(void)
{
    const double ring[] = {0, 1, 1, 1, 2, 1, 2, 0, 1};
    ListGraph *g = BUILD(3, true, ring);
    double r[3];
    size_t it = 0;
    CHECK_OK(list_adv_pagerank(g, 0.85, 1e-12, 200, false, r, &it));
    CHECK(near(r[0], 1.0 / 3) && near(r[1], 1.0 / 3) && near(r[2], 1.0 / 3));
    CHECK(list_adv_pagerank(g, 1.0, 1e-9, 10, false, r, NULL) == GRAPH_ERR_INVALID_ARGUMENT);
    CHECK(list_adv_pagerank(g, 0.85, 1e-9, 0, false, r, NULL) == GRAPH_ERR_INVALID_ARGUMENT);
    list_graph_destroy(g);

    const double in[] = {1, 0, 1, 2, 0, 1, 3, 0, 1};
    g = BUILD(4, true, in);
    double q[4];
    CHECK_OK(list_adv_pagerank(g, 0.85, 1e-12, 500, false, q, NULL));
    CHECK(q[0] > q[1] && near(q[1], q[2]) && near(q[1], q[3]));
    CHECK(near(q[0] + q[1] + q[2] + q[3], 1.0));
    list_graph_destroy(g);
}

static void test_reverse_delete(void)
{
    const double e[] = {0, 1, 1, 1, 2, 2, 0, 2, 3, 2, 3, 4};
    ListGraph *g = BUILD(5, false, e);   /* vertex 4 is isolated */
    ListEdge edges[4];
    size_t count = 0, comps = 0;
    double total = 0;
    CHECK_OK(list_adv_reverse_delete(g, edges, &count, &total, &comps));
    CHECK(count == 3 && comps == 2 && near(total, 7.0));
    CHECK(near(edges[0].weight, 1.0) && near(edges[1].weight, 2.0) && near(edges[2].weight, 4.0));

    /* Same weight as the library's own MST on random connected-ish graphs. */
    for (unsigned seed = 1; seed <= 15; ++seed) {
        ListGraph *r = random_graph(12, 0.4, false, 9, seed, 0);
        ListEdge mst[12];
        size_t mc = 0, rc = 0;
        double rd_w = 0;
        CHECK_OK(list_adv_reverse_delete(r, mst, &mc, &rd_w, &rc));
        double sum = 0;
        for (size_t i = 0; i < mc; ++i) {
            sum += mst[i].weight;
        }
        CHECK(near(sum, rd_w) && mc + rc == 12);
        list_graph_destroy(r);
    }
    list_graph_destroy(g);

    ListGraph *d = list_graph_create(3, true);
    CHECK(list_adv_reverse_delete(d, edges, &count, &total, &comps) == GRAPH_ERR_UNSUPPORTED);
    list_graph_destroy(d);
}

static void test_rocha_thatte(void)
{
    const double e[] = {0, 1, 1, 1, 2, 1, 2, 0, 1, 1, 3, 1, 3, 3, 1};
    ListGraph *g = BUILD(4, true, e);
    ListAdvCycles cy;
    CHECK_OK(list_adv_rocha_thatte(g, 0, 0, &cy));
    CHECK(cy.count == 2 && !cy.truncated);
    size_t lens[2] = {cy.offsets[1] - cy.offsets[0], cy.offsets[2] - cy.offsets[1]};
    CHECK((lens[0] == 3 && lens[1] == 1) || (lens[0] == 1 && lens[1] == 3));
    list_adv_cycles_free(&cy);

    CHECK_OK(list_adv_rocha_thatte(g, 1, 0, &cy));      /* only the self-loop */
    CHECK(cy.count == 1 && cy.members[0] == 3);
    list_adv_cycles_free(&cy);

    CHECK_OK(list_adv_rocha_thatte(g, 0, 1, &cy));      /* stop after one */
    CHECK(cy.count == 1 && cy.truncated);
    list_adv_cycles_free(&cy);

    ListGraph *u = list_graph_create(3, false);
    CHECK(list_adv_rocha_thatte(u, 0, 0, &cy) == GRAPH_ERR_UNSUPPORTED);
    list_graph_destroy(u);
    list_graph_destroy(g);
}

static void test_tree_traversal(void)
{
    const double e[] = {0, 1, 1, 0, 2, 1, 1, 3, 1, 1, 4, 1};
    ListGraph *g = BUILD(5, false, e);
    int out[5];
    size_t cnt = 0;
    const int pre[]   = {0, 1, 3, 4, 2};
    const int in[]    = {3, 1, 4, 0, 2};
    const int post[]  = {3, 4, 1, 2, 0};
    const int level[] = {0, 1, 2, 3, 4};
    CHECK_OK(list_adv_tree_traversal(g, 0, LIST_ADV_TREE_PREORDER, out, &cnt));
    CHECK(cnt == 5 && memcmp(out, pre, sizeof(pre)) == 0);
    CHECK_OK(list_adv_tree_traversal(g, 0, LIST_ADV_TREE_INORDER, out, &cnt));
    CHECK(memcmp(out, in, sizeof(in)) == 0);
    CHECK_OK(list_adv_tree_traversal(g, 0, LIST_ADV_TREE_POSTORDER, out, &cnt));
    CHECK(memcmp(out, post, sizeof(post)) == 0);
    CHECK_OK(list_adv_tree_traversal(g, 0, LIST_ADV_TREE_LEVELORDER, out, &cnt));
    CHECK(memcmp(out, level, sizeof(level)) == 0);

    CHECK(list_graph_add_edge(g, 3, 4, 1.0));            /* now a cycle, not a tree */
    CHECK(list_adv_tree_traversal(g, 0, LIST_ADV_TREE_PREORDER, out, &cnt) == GRAPH_ERR_UNSUPPORTED);
    list_graph_destroy(g);
}

static void test_zero_weight_cycle(void)
{
    ListAdvCycleKind kind;
    int cyc[3];
    size_t len = 0;
    const double zero[] = {0, 1, 2, 1, 2, -3, 2, 0, 1};
    ListGraph *g = BUILD(3, true, zero);
    CHECK_OK(list_adv_zero_weight_cycle(g, 1e-9, &kind, cyc, &len));
    CHECK(kind == LIST_ADV_CYCLE_ZERO && len == 3);
    list_graph_destroy(g);

    const double pos[] = {0, 1, 1, 1, 2, 1, 2, 0, 1};
    g = BUILD(3, true, pos);
    CHECK_OK(list_adv_zero_weight_cycle(g, 1e-9, &kind, cyc, &len));
    CHECK(kind == LIST_ADV_CYCLE_NONE);
    list_graph_destroy(g);

    const double neg[] = {0, 1, 1, 1, 2, 1, 2, 0, -3};
    g = BUILD(3, true, neg);
    CHECK_OK(list_adv_zero_weight_cycle(g, 1e-9, &kind, cyc, &len));
    CHECK(kind == LIST_ADV_CYCLE_NEGATIVE);
    list_graph_destroy(g);
}

static void test_subgraph_isomorphism(void)
{
    const double k4[] = {0, 1, 1, 0, 2, 1, 0, 3, 1, 1, 2, 1, 1, 3, 1, 2, 3, 1};
    const double tri[] = {0, 1, 1, 1, 2, 1, 0, 2, 1};
    const double path3[] = {0, 1, 1, 1, 2, 1};
    const double star[] = {0, 1, 1, 0, 2, 1, 0, 3, 1};
    ListGraph *host = BUILD(4, false, k4);
    ListGraph *t = BUILD(3, false, tri);
    ListGraph *p = BUILD(3, false, path3);
    ListGraph *s = BUILD(4, false, star);
    bool found = false;
    int map[4];

    CHECK_OK(list_adv_subgraph_isomorphism(t, host, false, false, map, &found));
    CHECK(found);
    for (int a = 0; a < 3; ++a) {
        for (int b = a + 1; b < 3; ++b) {
            CHECK(list_graph_has_edge(host, map[a], map[b], NULL));
        }
    }
    CHECK_OK(list_adv_subgraph_isomorphism(p, t, false, false, NULL, &found));
    CHECK(found);                                       /* path fits in a triangle */
    CHECK_OK(list_adv_subgraph_isomorphism(p, t, true, false, NULL, &found));
    CHECK(!found);                                      /* but not as an induced subgraph */
    CHECK_OK(list_adv_subgraph_isomorphism(t, s, false, false, NULL, &found));
    CHECK(!found);                                      /* a star has no triangle */
    list_graph_destroy(host);
    list_graph_destroy(t);
    list_graph_destroy(p);
    list_graph_destroy(s);
}

static void test_coffman_graham(void)
{
    int label[6], slot[6];
    size_t span = 0;

    ListGraph *g = list_graph_create(6, true);          /* six independent tasks */
    CHECK_OK(list_adv_coffman_graham(g, 2, label, slot, &span));
    CHECK(span == 3);
    CHECK_OK(list_adv_coffman_graham(g, 3, label, slot, &span));
    CHECK(span == 2);
    CHECK(list_adv_coffman_graham(g, 0, label, slot, &span) == GRAPH_ERR_INVALID_ARGUMENT);
    list_graph_destroy(g);

    const double two_chains[] = {0, 1, 1, 1, 2, 1, 3, 4, 1, 4, 5, 1};
    g = BUILD(6, true, two_chains);                     /* 0->1->2 and 3->4->5 */
    CHECK_OK(list_adv_coffman_graham(g, 2, label, slot, &span));
    CHECK(span == 3);
    for (int v = 0; v < 6; ++v) {
        CHECK(slot[v] >= 0 && slot[v] < 3);
    }
    CHECK(slot[0] < slot[1] && slot[1] < slot[2] && slot[3] < slot[4] && slot[4] < slot[5]);
    CHECK_OK(list_adv_coffman_graham(g, 1, label, slot, &span));
    CHECK(span == 6);
    CHECK(list_graph_add_edge(g, 2, 0, 1.0));
    CHECK(list_adv_coffman_graham(g, 2, label, slot, &span) == GRAPH_ERR_CYCLE);
    list_graph_destroy(g);
}

static void test_crossings(void)
{
    const double k4[] = {0, 1, 1, 0, 2, 1, 0, 3, 1, 1, 2, 1, 1, 3, 1, 2, 3, 1};
    ListGraph *g = BUILD(4, false, k4);
    long long c = -1;
    int order[4];
    CHECK_OK(list_adv_crossing_number_heuristic(g, 0, &c, order));
    CHECK(c == 1);                                      /* every circular order of K4 has 1 */
    list_graph_destroy(g);

    const double c4[] = {0, 1, 1, 1, 2, 1, 2, 3, 1, 3, 0, 1};
    g = BUILD(4, false, c4);
    const int good[] = {0, 1, 2, 3};
    const int bad[]  = {0, 2, 1, 3};
    const int dup[]  = {0, 1, 1, 3};
    CHECK_OK(list_adv_crossing_count(g, good, 4, &c));
    CHECK(c == 0);
    CHECK_OK(list_adv_crossing_count(g, bad, 4, &c));
    CHECK(c == 1);
    CHECK(list_adv_crossing_count(g, dup, 4, &c) == GRAPH_ERR_INVALID_VERTEX);
    CHECK(list_adv_crossing_count(g, good, 3, &c) == GRAPH_ERR_INVALID_ARGUMENT);
    CHECK_OK(list_adv_crossing_number_heuristic(g, 0, &c, order));
    CHECK(c == 0);                                      /* the heuristic finds the planar order */
    list_graph_destroy(g);

    /* K5: C(5,4) = 5 crossings in every circular order. */
    g = list_graph_create(5, false);
    for (int a = 0; a < 5; ++a) {
        for (int b = a + 1; b < 5; ++b) {
            CHECK(list_graph_add_edge(g, a, b, 1.0));
        }
    }
    int o5[5];
    CHECK_OK(list_adv_crossing_number_heuristic(g, 0, &c, o5));
    CHECK(c == 5);
    long long again = 0;
    CHECK_OK(list_adv_crossing_count(g, o5, 5, &again));
    CHECK(again == c);                                  /* the two functions agree */
    list_graph_destroy(g);
}

static void test_color_refinement(void)
{
    const double p3[] = {0, 1, 1, 1, 2, 1};
    ListGraph *g = BUILD(3, false, p3);
    int col[3];
    size_t nc = 0, rounds = 0;
    CHECK_OK(list_adv_color_refinement(g, col, &nc, &rounds));
    CHECK(nc == 2 && col[0] == col[2] && col[0] != col[1]);
    list_graph_destroy(g);

    const double c4[] = {0, 1, 1, 1, 2, 1, 2, 3, 1, 3, 0, 1};
    g = BUILD(4, false, c4);
    int c4col[4];
    CHECK_OK(list_adv_color_refinement(g, c4col, &nc, NULL));
    CHECK(nc == 1);
    list_graph_destroy(g);
}

static void test_color_coding(void)
{
    const double p5[] = {0, 1, 1, 1, 2, 1, 2, 3, 1, 3, 4, 1};
    ListGraph *g = BUILD(5, false, p5);
    bool found = false;
    int path[4];
    CHECK_OK(list_adv_color_coding_path(g, 4, 0, 12345, &found, path));
    CHECK(found);
    for (int i = 0; i < 3; ++i) {
        CHECK(list_graph_has_edge(g, path[i], path[i + 1], NULL));
    }
    for (int i = 0; i < 4; ++i) {
        for (int j = i + 1; j < 4; ++j) {
            CHECK(path[i] != path[j]);
        }
    }
    CHECK_OK(list_adv_color_coding_path(g, 6, 50, 7, &found, NULL));
    CHECK(!found);                                      /* more vertices than the graph has */
    CHECK(list_adv_color_coding_path(g, 0, 1, 1, &found, NULL) == GRAPH_ERR_INVALID_ARGUMENT);
    list_graph_destroy(g);
}

static void test_centrality(void)
{
    const double star[] = {0, 1, 1, 0, 2, 1, 0, 3, 1};
    ListGraph *g = BUILD(4, false, star);
    double h[4], ev[4], deg[4];
    CHECK_OK(list_adv_harmonic_centrality(g, h));
    CHECK(near(h[0], 1.0) && near(h[1], 2.0 / 3));
    size_t it = 0;
    CHECK_OK(list_adv_eigenvector_centrality(g, 0, 0.0, ev, &it));
    CHECK(ev[0] > ev[1] && near(ev[1], ev[2]));
    CHECK(near(ev[0] * ev[0] + ev[1] * ev[1] + ev[2] * ev[2] + ev[3] * ev[3], 1.0));
    CHECK_OK(list_adv_centrality(g, LIST_ADV_CENTRALITY_DEGREE, deg));
    CHECK(deg[0] > deg[1]);
    CHECK_OK(list_adv_centrality(g, LIST_ADV_CENTRALITY_HARMONIC, ev));
    CHECK(near(ev[0], h[0]));
    list_graph_destroy(g);
}

static void test_shortest_path_family(void)
{
    for (unsigned seed = 1; seed <= 20; ++seed) {
        ListGraph *g = random_graph(14, 0.25, seed % 2 == 0, 9, seed, seed % 4);
        const size_t n = list_graph_id_bound(g);
        double dj[14];
        int par[14];
        int s = 0;
        while (s < (int)n && !list_graph_is_alive(g, s)) {
            ++s;
        }
        if (s == (int)n) {
            list_graph_destroy(g);
            continue;
        }
        CHECK(list_dijkstra(g, s, dj, par));

        double fw[14 * 14], jo[14 * 14];
        CHECK_OK(list_adv_floyd_warshall(g, fw, NULL));
        CHECK_OK(list_adv_johnson(g, jo, NULL));
        for (size_t i = 0; i < n * n; ++i) {
            CHECK(near(fw[i], jo[i]));
        }
        for (size_t t = 0; t < n; ++t) {
            CHECK(near(fw[(size_t)s * n + t], dj[t]));
            if (!list_graph_is_alive(g, (int)t) || isinf(dj[t]) || (int)t == s) {
                continue;
            }
            ListAdvPath a, d, b;
            size_t exp = 0;
            CHECK_OK(list_adv_astar(g, s, (int)t, zero_h, NULL, &a, &exp));
            CHECK(near(a.cost, dj[t]));
            check_path(g, &a, s, (int)t);
            CHECK_OK(list_adv_ida_star(g, s, (int)t, zero_h, NULL, 0, &d));
            CHECK(near(d.cost, dj[t]));
            check_path(g, &d, s, (int)t);
            CHECK_OK(list_adv_iddfs(g, s, (int)t, SIZE_MAX, &b));
            check_path(g, &b, s, (int)t);
            ListAdvPath k[3];
            size_t found = 0;
            CHECK_OK(list_adv_k_shortest_paths(g, s, (int)t, 3, k, &found));
            CHECK(found >= 1 && near(k[0].cost, dj[t]));
            for (size_t i = 1; i < found; ++i) {
                CHECK(k[i].cost >= k[i - 1].cost - 1e-9);
            }
            for (size_t i = 0; i < found; ++i) {
                list_adv_path_free(&k[i]);
            }
            list_adv_path_free(&a);
            list_adv_path_free(&d);
            list_adv_path_free(&b);
        }
        list_graph_destroy(g);
    }
}

static void test_negative_weights(void)
{
    const double cyc[] = {0, 1, 1, 1, 0, -2};
    ListGraph *g = BUILD(2, true, cyc);
    double d[4];
    CHECK(list_adv_floyd_warshall(g, d, NULL) == GRAPH_ERR_NEGATIVE_CYCLE);
    CHECK(list_adv_johnson(g, d, NULL) == GRAPH_ERR_NEGATIVE_CYCLE);
    list_graph_destroy(g);

    const double dag[] = {0, 1, 4, 0, 2, 5, 2, 1, -3, 1, 3, 2};
    g = BUILD(4, true, dag);
    double fw[16], jo[16];
    int pred[16];
    CHECK_OK(list_adv_floyd_warshall(g, fw, pred));
    CHECK_OK(list_adv_johnson(g, jo, NULL));
    for (int i = 0; i < 16; ++i) {
        CHECK(near(fw[i], jo[i]));
    }
    CHECK(near(fw[0 * 4 + 1], 2.0) && near(fw[0 * 4 + 3], 4.0));
    ListAdvPath p;
    CHECK_OK(list_adv_apsp_path(fw, pred, 4, 0, 3, &p));
    CHECK(p.length == 4 && p.vertices[1] == 2 && p.vertices[2] == 1);
    list_adv_path_free(&p);
    list_graph_destroy(g);
}

static void test_astar_unreachable_and_limits(void)
{
    const double e[] = {0, 1, 1};
    ListGraph *g = BUILD(3, true, e);
    ListAdvPath p;
    CHECK(list_adv_astar(g, 0, 2, NULL, NULL, &p, NULL) == GRAPH_ERR_NOT_FOUND);
    CHECK(list_adv_ida_star(g, 0, 2, NULL, NULL, 0, &p) == GRAPH_ERR_NOT_FOUND);
    CHECK(list_adv_iddfs(g, 0, 2, SIZE_MAX, &p) == GRAPH_ERR_NOT_FOUND);
    CHECK(list_adv_astar(g, 0, 7, NULL, NULL, &p, NULL) == GRAPH_ERR_INVALID_VERTEX);
    list_graph_destroy(g);

    /* iddfs minimises edges, not weight: the direct heavy edge wins. */
    const double w[] = {0, 1, 1, 1, 2, 1, 2, 3, 1, 0, 3, 50};
    g = BUILD(4, false, w);
    CHECK_OK(list_adv_iddfs(g, 0, 3, SIZE_MAX, &p));
    CHECK(p.length == 2 && near(p.cost, 50.0));
    list_adv_path_free(&p);
    CHECK(list_adv_iddfs(g, 0, 3, 0, &p) == GRAPH_ERR_NOT_FOUND);
    list_graph_destroy(g);
}

static void test_flood(void)
{
    const double e[] = {0, 1, 1, 1, 2, 1, 2, 3, 1};
    ListGraph *g = BUILD(4, false, e);
    int round[4], parent[4];
    size_t msgs = 0, reached = 0;
    CHECK_OK(list_adv_flood(g, 0, SIZE_MAX, round, parent, &msgs, &reached));
    CHECK(round[0] == 0 && round[1] == 1 && round[2] == 2 && round[3] == 3);
    CHECK(reached == 4 && parent[0] == -1 && parent[3] == 2);
    CHECK_OK(list_adv_flood(g, 0, 1, round, NULL, &msgs, &reached));
    CHECK(reached == 2 && round[2] == -1);
    CHECK_OK(list_adv_flood(g, 0, 0, round, NULL, NULL, &reached));
    CHECK(reached == 1);
    list_graph_destroy(g);
}

static void test_layouts(void)
{
    const double e[] = {0, 1, 1, 1, 2, 1, 2, 3, 1, 3, 0, 1};
    ListGraph *g = BUILD(5, false, e);                  /* vertex 4 isolated */
    double x[5], y[5], x2[5], y2[5];
    CHECK_OK(list_adv_force_layout(g, 100, x, y));
    CHECK_OK(list_adv_force_layout(g, 100, x2, y2));
    for (int i = 0; i < 5; ++i) {
        CHECK(x[i] >= 0.0 && x[i] <= 1.0 && y[i] >= 0.0 && y[i] <= 1.0);
        CHECK(x[i] == x2[i] && y[i] == y2[i]);          /* deterministic */
    }
    CHECK(fabs(x[0] - x[1]) + fabs(y[0] - y[1]) > 1e-6);

    CHECK_OK(list_adv_circle_layout(g, 2.0, x, y));
    CHECK(near(x[0], 2.0) && fabs(y[0]) < 1e-12);
    for (int i = 0; i < 5; ++i) {
        CHECK(near(hypot(x[i], y[i]), 2.0));
    }
    CHECK(list_adv_circle_layout(g, -1.0, x, y) == GRAPH_ERR_INVALID_ARGUMENT);

    CHECK(list_graph_remove_vertex(g, 4));
    CHECK_OK(list_adv_circle_layout(g, 1.0, x, y));
    CHECK(x[4] == 0.0 && y[4] == 0.0 && near(y[1], 1.0));
    list_graph_destroy(g);
}

static void test_gnn(void)
{
    const double e[] = {0, 1, 1};
    ListGraph *g = BUILD(2, false, e);
    const double feat[] = {1.0, 2.0};
    const double w[]    = {1.0};
    ListAdvGnnLayer layer = {1, 1, w, NULL, LIST_ADV_ACT_NONE};
    double out[2];
    CHECK_OK(list_adv_gnn_forward(g, feat, 1, &layer, 1, out));
    CHECK(near(out[0], 1.5) && near(out[1], 1.5));      /* both average (1 + 2) / 2 */

    const double neg[] = {-1.0};
    const double bias[] = {0.25};
    ListAdvGnnLayer relu = {1, 1, neg, bias, LIST_ADV_ACT_RELU};
    CHECK_OK(list_adv_gnn_forward(g, feat, 1, &relu, 1, out));
    CHECK(near(out[0], 0.0) && near(out[1], 0.0));
    layer.in_dim = 2;
    CHECK(list_adv_gnn_forward(g, feat, 1, &layer, 1, out) == GRAPH_ERR_INVALID_ARGUMENT);
    list_graph_destroy(g);
}

static void test_isomorphism(void)
{
    const double c5[]  = {0, 1, 1, 1, 2, 1, 2, 3, 1, 3, 4, 1, 4, 0, 1};
    const double c5b[] = {0, 2, 1, 2, 4, 1, 4, 1, 1, 1, 3, 1, 3, 0, 1};
    ListGraph *a = BUILD(5, false, c5);
    ListGraph *b = BUILD(5, false, c5b);
    bool iso = false;
    int map[5];
    CHECK_OK(list_adv_graph_isomorphism(a, b, false, &iso, map));
    CHECK(iso);
    for (int u = 0; u < 5; ++u) {
        for (int v = 0; v < 5; ++v) {
            CHECK(list_graph_has_edge(a, u, v, NULL) == list_graph_has_edge(b, map[u], map[v], NULL));
        }
    }
    list_graph_destroy(a);
    list_graph_destroy(b);

    /* Both are 2-regular on 6 vertices, so colour refinement alone cannot tell them apart. */
    const double c6[]  = {0, 1, 1, 1, 2, 1, 2, 3, 1, 3, 4, 1, 4, 5, 1, 5, 0, 1};
    const double tt[]  = {0, 1, 1, 1, 2, 1, 2, 0, 1, 3, 4, 1, 4, 5, 1, 5, 3, 1};
    a = BUILD(6, false, c6);
    b = BUILD(6, false, tt);
    CHECK_OK(list_adv_graph_isomorphism(a, b, false, &iso, NULL));
    CHECK(!iso);
    list_graph_destroy(a);
    list_graph_destroy(b);

    /* Weights matter only when asked for. */
    const double w1[] = {0, 1, 1, 1, 2, 2};
    const double w2[] = {0, 1, 2, 1, 2, 2};
    a = BUILD(3, false, w1);
    b = BUILD(3, false, w2);
    CHECK_OK(list_adv_graph_isomorphism(a, b, false, &iso, NULL));
    CHECK(iso);
    CHECK_OK(list_adv_graph_isomorphism(a, b, true, &iso, NULL));
    CHECK(!iso);
    list_graph_destroy(a);
    list_graph_destroy(b);
}

static void test_jump_point_search(void)
{
    ListGraph *g = list_graph_create(25, false);        /* 5 x 5 grid, edges unused */
    ListAdvPath p;
    CHECK_OK(list_adv_jump_point_search(g, 5, 0, 24, &p));
    CHECK(near(p.cost, 4.0 * sqrt(2.0)) && p.vertices[0] == 0 && p.vertices[p.length - 1] == 24);
    list_adv_path_free(&p);

    CHECK(list_graph_remove_vertex(g, 7));              /* wall at x = 2, y = 1..3 */
    CHECK(list_graph_remove_vertex(g, 12));
    CHECK(list_graph_remove_vertex(g, 17));
    CHECK_OK(list_adv_jump_point_search(g, 5, 0, 24, &p));
    CHECK(near(p.cost, 4.0 + 2.0 * sqrt(2.0)));
    for (size_t i = 0; i < p.length; ++i) {
        CHECK(p.vertices[i] != 7 && p.vertices[i] != 12 && p.vertices[i] != 17);
    }
    list_adv_path_free(&p);
    CHECK(list_adv_jump_point_search(g, 5, 7, 24, &p) == GRAPH_ERR_INVALID_VERTEX);
    CHECK(list_adv_jump_point_search(g, 0, 0, 24, &p) == GRAPH_ERR_INVALID_ARGUMENT);
    list_graph_destroy(g);
}

static void test_k_shortest(void)
{
    const double e[] = {0, 1, 1, 1, 3, 1, 0, 2, 1, 2, 3, 2, 0, 3, 5};
    ListGraph *g = BUILD(4, true, e);
    ListAdvPath p[5];
    size_t found = 0;
    CHECK_OK(list_adv_k_shortest_paths(g, 0, 3, 5, p, &found));
    CHECK(found == 3 && near(p[0].cost, 2.0) && near(p[1].cost, 3.0) && near(p[2].cost, 5.0));
    for (size_t i = 0; i < found; ++i) {
        check_path(g, &p[i], 0, 3);
        list_adv_path_free(&p[i]);
    }
    CHECK(list_adv_k_shortest_paths(g, 3, 0, 2, p, &found) == GRAPH_ERR_NOT_FOUND && found == 0);
    list_graph_destroy(g);
}

static void test_widest_path(void)
{
    const double e[] = {0, 1, 5, 1, 2, 3, 0, 2, 2};
    ListGraph *g = BUILD(4, false, e);
    double width[4];
    int parent[4];
    CHECK_OK(list_adv_widest_path(g, 0, width, parent));
    CHECK(isinf(width[0]) && width[0] > 0);
    CHECK(near(width[1], 5.0) && near(width[2], 3.0) && parent[2] == 1);
    CHECK(isinf(width[3]) && width[3] < 0 && parent[3] == -1);
    ListAdvPath p;
    CHECK_OK(list_adv_path_from_parents(parent, 4, 0, 2, width[2], &p));
    CHECK(p.length == 3 && p.vertices[1] == 1);
    list_adv_path_free(&p);
    list_graph_destroy(g);
}

static void test_lex_bfs(void)
{
    const double p3[] = {0, 1, 1, 1, 2, 1};
    ListGraph *g = BUILD(3, false, p3);
    int order[3];
    CHECK_OK(list_adv_lex_bfs(g, 0, order));
    CHECK(order[0] == 0 && order[1] == 1 && order[2] == 2);
    list_graph_destroy(g);

    for (unsigned seed = 1; seed <= 10; ++seed) {
        g = random_graph(12, 0.3, false, 5, seed, 2);
        int ord[12];
        const size_t vc = list_graph_vertex_count(g);
        CHECK_OK(list_adv_lex_bfs(g, -1, ord));
        bool seen[12] = {false};
        for (size_t i = 0; i < vc; ++i) {
            CHECK(ord[i] >= 0 && ord[i] < 12 && list_graph_is_alive(g, ord[i]) && !seen[ord[i]]);
            seen[ord[i]] = true;
        }
        list_graph_destroy(g);
    }
}

static void test_dynamic_connectivity(void)
{
    ListGraph *g = list_graph_create(4, false);
    CHECK(list_graph_add_edge(g, 0, 1, 1.0));
    ListAdvDynConn *dc = list_adv_dynconn_create(g);
    CHECK(dc != NULL);
    bool c = false;
    size_t n = 0;
    CHECK_OK(list_adv_dynconn_connected(dc, 0, 1, &c));
    CHECK(c);
    CHECK_OK(list_adv_dynconn_connected(dc, 0, 2, &c));
    CHECK(!c);
    CHECK_OK(list_adv_dynconn_components(dc, &n));
    CHECK(n == 3);

    CHECK(list_graph_add_edge(g, 1, 2, 1.0));           /* the tracker notices on its own */
    CHECK_OK(list_adv_dynconn_connected(dc, 0, 2, &c));
    CHECK(c);
    CHECK_OK(list_adv_dynconn_components(dc, &n));
    CHECK(n == 2);

    CHECK(list_graph_remove_edge(g, 0, 1));
    CHECK_OK(list_adv_dynconn_connected(dc, 0, 2, &c));
    CHECK(!c);

    CHECK(list_graph_remove_vertex(g, 3));
    CHECK_OK(list_adv_dynconn_components(dc, &n));
    CHECK(n == 2);
    CHECK(list_adv_dynconn_connected(dc, 0, 3, &c) == GRAPH_ERR_INVALID_VERTEX);
    CHECK(list_adv_dynconn_connected(NULL, 0, 1, &c) == GRAPH_ERR_INVALID_ARGUMENT);
    CHECK(list_adv_dynconn_create(NULL) == NULL);
    list_adv_dynconn_destroy(dc);
    list_adv_dynconn_destroy(NULL);
    list_graph_destroy(g);
}

int main(void)
{
    test_longest_path();
    test_pagerank();
    test_reverse_delete();
    test_rocha_thatte();
    test_tree_traversal();
    test_zero_weight_cycle();
    test_subgraph_isomorphism();
    test_coffman_graham();
    test_crossings();
    test_color_refinement();
    test_color_coding();
    test_centrality();
    test_shortest_path_family();
    test_negative_weights();
    test_astar_unreachable_and_limits();
    test_flood();
    test_layouts();
    test_gnn();
    test_isomorphism();
    test_jump_point_search();
    test_k_shortest();
    test_widest_path();
    test_lex_bfs();
    test_dynamic_connectivity();
    puts("All list advanced-algorithm tests passed.");
    return 0;
}
