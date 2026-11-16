/*
 * test_graph_advanced.c — the unified wrappers for the advanced algorithms
 * (graph_advanced.h).  Every wrapper is run on the same graph in the list and
 * in the matrix representation and checked against a known answer, so a
 * wrapper that dispatches to the wrong backend, or a backend that disagrees
 * with the other, is caught.  Only values that are unique are compared
 * (distances, counts, colours); parent arrays may differ in tie-breaking.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "graph.h"
#include "graph_advanced.h"
#include "graph_algorithms.h"

#define CHECK(cond)                                                                  \
    do {                                                                             \
        if (!(cond)) {                                                               \
            fprintf(stderr, "%s:%d: CHECK failed: %s\n", __FILE__, __LINE__, #cond); \
            exit(1);                                                                 \
        }                                                                            \
    } while (0)

#define OK(call) CHECK((call) == GRAPH_OK)

static bool near(double a, double b)
{
    return fabs(a - b) < 1e-9 || (isinf(a) && isinf(b) && (a > 0) == (b > 0));
}

static const GraphRep REPS[2] = {GRAPH_LIST, GRAPH_MATRIX};

/* A graph with fixed edges in a given representation (adaptive off, so the
 * representation stays what the test asked for). */
static Graph *make(GraphRep rep, size_t n, bool directed, const double (*edges)[3], size_t m)
{
    Graph *g = graph_create_sized(n, directed);
    CHECK(g != NULL);
    graph_set_adaptive(g, false);
    for (size_t i = 0; i < m; i++) {
        CHECK(graph_add_edge(g, (int)edges[i][0], (int)edges[i][1], edges[i][2]));
    }
    CHECK(graph_convert(g, rep));
    CHECK(graph_representation(g) == rep);
    return g;
}

static void test_dag(void)
{
    /* 0->1 (2), 1->2 (3), 0->2 (1), 2->3 (4) */
    const double e[][3] = {{0, 1, 2}, {1, 2, 3}, {0, 2, 1}, {2, 3, 4}};
    for (int r = 0; r < 2; r++) {
        Graph *g = make(REPS[r], 4, true, e, 4);
        double dist[4];
        int    parent[4];
        OK(graph_adv_longest_path(g, 0, dist, parent));
        CHECK(near(dist[0], 0) && near(dist[1], 2) && near(dist[2], 5) && near(dist[3], 9));
        CHECK(parent[3] == 2 && parent[2] == 1 && parent[1] == 0 && parent[0] == -1);

        /* undirected input is rejected, bad source reported */
        CHECK(graph_adv_longest_path(g, 9, dist, parent) == GRAPH_ERR_INVALID_VERTEX);

        /* Coffman-Graham on 2 processors: 4 tasks in a chain-ish DAG */
        int label[4], slot[4];
        size_t makespan = 0;
        OK(graph_adv_coffman_graham(g, 2, label, slot, &makespan));
        CHECK(makespan == 4);   /* 0 < 1 < 2 < 3 is a chain after reduction */
        graph_destroy(g);
    }
    CHECK(graph_adv_longest_path(NULL, 0, NULL, NULL) == GRAPH_ERR_INVALID_ARGUMENT);
}

static void test_cycles(void)
{
    /* 0->1->2->0 (weights 1, 2, -3 = a zero cycle) and 3->3 self-loop */
    const double e[][3] = {{0, 1, 1}, {1, 2, 2}, {2, 0, -3}, {3, 3, 1}};
    const double pos[][3] = {{0, 1, 1}, {1, 2, 2}, {2, 0, 3}};
    const double neg[][3] = {{0, 1, 1}, {1, 2, 2}, {2, 0, -5}};
    for (int r = 0; r < 2; r++) {
        Graph *g = make(REPS[r], 4, true, e, 4);
        GraphAdvCycles cyc = {0};
        OK(graph_adv_rocha_thatte(g, 0, 0, &cyc));
        CHECK(cyc.count == 2 && !cyc.truncated);   /* the triangle and the self-loop */
        size_t lens[2] = {cyc.offsets[1] - cyc.offsets[0], cyc.offsets[2] - cyc.offsets[1]};
        CHECK((lens[0] == 3 && lens[1] == 1) || (lens[0] == 1 && lens[1] == 3));
        graph_adv_cycles_free(&cyc);
        CHECK(cyc.members == NULL && cyc.offsets == NULL && cyc.count == 0);
        graph_adv_cycles_free(NULL);

        GraphAdvCycles one = {0};
        OK(graph_adv_rocha_thatte(g, 0, 1, &one));
        CHECK(one.count == 1 && one.truncated);
        graph_adv_cycles_free(&one);

        GraphAdvCycleKind kind;
        int    cycle[4];
        size_t len = 0;
        /* the self-loop has weight 1 and the triangle 0, so a zero cycle exists */
        OK(graph_adv_zero_weight_cycle(g, 1e-9, &kind, cycle, &len));
        CHECK(kind == GRAPH_ADV_CYCLE_ZERO && len == 3);
        graph_destroy(g);

        g = make(REPS[r], 3, true, pos, 3);
        OK(graph_adv_zero_weight_cycle(g, 1e-9, &kind, cycle, &len));
        CHECK(kind == GRAPH_ADV_CYCLE_NONE);
        graph_destroy(g);

        g = make(REPS[r], 3, true, neg, 3);
        OK(graph_adv_zero_weight_cycle(g, 1e-9, &kind, cycle, &len));
        CHECK(kind == GRAPH_ADV_CYCLE_NEGATIVE);
        graph_destroy(g);

        /* undirected graphs are unsupported */
        g = make(REPS[r], 3, false, pos, 3);
        CHECK(graph_adv_rocha_thatte(g, 0, 0, &one) == GRAPH_ERR_UNSUPPORTED);
        graph_destroy(g);
    }
}

static void test_trees_and_mst(void)
{
    /* tree: 0-1, 0-2, 1-3, 1-4 */
    const double t[][3] = {{0, 1, 1}, {0, 2, 1}, {1, 3, 1}, {1, 4, 1}};
    /* triangle with weights 1, 2, 3 plus an isolated vertex 3 */
    const double tri[][3] = {{0, 1, 1}, {1, 2, 2}, {0, 2, 3}};
    for (int r = 0; r < 2; r++) {
        Graph *g = make(REPS[r], 5, false, t, 4);
        int    out[5];
        size_t count = 0;
        OK(graph_adv_tree_traversal(g, 0, GRAPH_ADV_TREE_PREORDER, out, &count));
        CHECK(count == 5);
        const int pre[5] = {0, 1, 3, 4, 2};
        CHECK(memcmp(out, pre, sizeof pre) == 0);
        OK(graph_adv_tree_traversal(g, 0, GRAPH_ADV_TREE_POSTORDER, out, &count));
        const int post[5] = {3, 4, 1, 2, 0};
        CHECK(memcmp(out, post, sizeof post) == 0);
        OK(graph_adv_tree_traversal(g, 0, GRAPH_ADV_TREE_LEVELORDER, out, &count));
        const int lvl[5] = {0, 1, 2, 3, 4};
        CHECK(memcmp(out, lvl, sizeof lvl) == 0);
        graph_destroy(g);

        g = make(REPS[r], 4, false, tri, 3);
        GraphEdge edges[4];
        size_t    ecount = 0, comps = 0;
        double    total = 0;
        OK(graph_adv_reverse_delete(g, edges, &ecount, &total, &comps));
        CHECK(ecount == 2 && near(total, 3.0) && comps == 2);   /* triangle + isolated vertex */
        /* Kruskal refuses a disconnected graph; reverse-delete returns a forest */
        GraphEdge ke[4];
        size_t    kc = 0;
        CHECK(isinf(graph_mst_kruskal(g, ke, &kc)));
        CHECK(graph_remove_vertex(g, 3));
        OK(graph_adv_reverse_delete(g, edges, &ecount, &total, &comps));
        CHECK(comps == 1 && near(total, graph_mst_kruskal(g, ke, &kc)));   /* same tree weight */
        graph_destroy(g);
    }
}

static void test_structure(void)
{
    /* path 0-1-2-3 plus a triangle 4-5-6 */
    const double e[][3] = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {4, 5, 1}, {5, 6, 1}, {4, 6, 1}};
    for (int r = 0; r < 2; r++) {
        Graph *g = make(REPS[r], 7, false, e, 6);
        int    color[7];
        size_t ncol = 0, rounds = 0;
        OK(graph_adv_color_refinement(g, color, &ncol, &rounds));
        /* path ends / path middles / triangle vertices */
        CHECK(ncol == 3);
        CHECK(color[0] == color[3] && color[1] == color[2] && color[4] == color[5] &&
              color[5] == color[6] && color[0] != color[1] && color[1] != color[4]);

        double h[7];
        OK(graph_adv_harmonic_centrality(g, h));
        CHECK(near(h[4], h[5]) && h[4] > 0 && h[0] < h[1]);

        bool found = false;
        int  path[3];
        OK(graph_adv_color_coding_path(g, 3, 0, 42, &found, path));
        CHECK(found);
        OK(graph_adv_color_coding_path(g, 5, 0, 42, &found, NULL));
        CHECK(!found);   /* no simple path on 5 vertices exists */

        int    order[7] = {0, 1, 2, 3, 4, 5, 6};
        long long crossings = -1;
        OK(graph_adv_crossing_count(g, order, 7, &crossings));
        CHECK(crossings == 0);
        long long heur = -1;
        OK(graph_adv_crossing_number_heuristic(g, 0, &heur, order));
        CHECK(heur == 0);

        int lex[7];
        OK(graph_adv_lex_bfs(g, 0, lex));
        bool seen[7] = {false};
        for (int i = 0; i < 7; i++) {
            CHECK(lex[i] >= 0 && lex[i] < 7 && !seen[lex[i]]);
            seen[lex[i]] = true;
        }
        graph_destroy(g);
    }

    /* K4 drawn on a circle in order 0,1,2,3 has exactly one crossing (the diagonals) */
    const double k4[][3] = {{0, 1, 1}, {0, 2, 1}, {0, 3, 1}, {1, 2, 1}, {1, 3, 1}, {2, 3, 1}};
    for (int r = 0; r < 2; r++) {
        Graph *g = make(REPS[r], 4, false, k4, 6);
        int    order[4] = {0, 1, 2, 3};
        long long c = 0;
        OK(graph_adv_crossing_count(g, order, 4, &c));
        CHECK(c == 1);
        graph_destroy(g);
    }
}

static double zero_h(int v, int goal, void *ctx)
{
    (void)v; (void)goal; (void)ctx;
    return 0.0;
}

static void test_paths(void)
{
    /* 0-1 (5), 1-2 (3), 0-2 (2), 2-3 (4): widest 0->3 is min(2,4) = 2 via 2, or min(5,3,4)=3 via 1 */
    const double e[][3] = {{0, 1, 5}, {1, 2, 3}, {0, 2, 2}, {2, 3, 4}};
    for (int r = 0; r < 2; r++) {
        Graph *g = make(REPS[r], 4, false, e, 4);
        double width[4];
        int    parent[4];
        OK(graph_adv_widest_path(g, 0, width, parent));
        CHECK(near(width[1], 5) && near(width[2], 3) && near(width[3], 3));

        GraphPath p = {0};
        OK(graph_adv_iddfs(g, 0, 3, SIZE_MAX, &p));
        CHECK(p.length == 3 && p.vertices[0] == 0 && p.vertices[2] == 3);   /* fewest edges: 0-2-3 */
        graph_path_free(&p);

        OK(graph_adv_ida_star(g, 0, 3, zero_h, NULL, 0, &p));
        CHECK(near(p.cost, 6.0));   /* 0-2-3 = 2 + 4 */
        graph_path_free(&p);

        /* the same path through the helper, from a parent array */
        int par[4] = {-1, 0, 0, 2};
        OK(graph_adv_path_from_parents(par, 4, 0, 3, 6.0, &p));
        CHECK(p.length == 3 && p.vertices[1] == 2 && near(p.cost, 6.0));
        graph_path_free(&p);

        CHECK(graph_adv_ida_star(g, 0, 3, NULL, NULL, 0, &p) == GRAPH_OK);
        graph_path_free(&p);
        graph_destroy(g);
    }

    /* Jump point search: 4x3 open grid with the middle column blocked at row 1 */
    for (int r = 0; r < 2; r++) {
        Graph *g = graph_create_sized(12, false);
        graph_set_adaptive(g, false);
        CHECK(graph_remove_vertex(g, 5));    /* (x=1, y=1) */
        CHECK(graph_convert(g, REPS[r]));
        GraphPath p = {0};
        OK(graph_adv_jump_point_search(g, 4, 0, 11, &p));
        CHECK(p.vertices[0] == 0 && p.vertices[p.length - 1] == 11);
        CHECK(near(p.cost, 2 * sqrt(2.0) + 1.0));   /* diagonal, diagonal, one step right */
        graph_path_free(&p);
        graph_destroy(g);
    }
}

static void test_dynamics(void)
{
    /* flood on a path of 4: reaches everyone, in rounds 0..3 */
    const double e[][3] = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}};
    for (int r = 0; r < 2; r++) {
        Graph *g = make(REPS[r], 4, false, e, 3);
        int    round[4], parent[4];
        size_t messages = 0, reached = 0;
        OK(graph_adv_flood(g, 0, 10, round, parent, &messages, &reached));
        CHECK(reached == 4 && round[3] == 3 && round[0] == 0);

        double x[4], y[4];
        OK(graph_adv_circle_layout(g, 1.0, x, y));
        for (int i = 0; i < 4; i++) {
            CHECK(near(x[i] * x[i] + y[i] * y[i], 1.0));
        }
        OK(graph_adv_force_layout(g, 50, x, y));
        for (int i = 0; i < 4; i++) {
            CHECK(x[i] >= 0 && x[i] <= 1 && y[i] >= 0 && y[i] <= 1);
        }

        /* one GCN layer with identity weights and no activation: output rows are
         * the averaged features, so equal features stay equal */
        const double feat[4] = {1, 1, 1, 1};
        const double w[1]    = {1.0};
        GraphAdvGnnLayer layer = {1, 1, w, NULL, GRAPH_ADV_ACT_NONE};
        double out[4];
        OK(graph_adv_gnn_forward(g, feat, 1, &layer, 1, out));
        for (int i = 0; i < 4; i++) {
            CHECK(out[i] > 0.5 && out[i] < 1.5);
        }
        graph_destroy(g);
    }
}

static void test_subgraph_isomorphism(void)
{
    const double k4[][3]   = {{0, 1, 1}, {0, 2, 1}, {0, 3, 1}, {1, 2, 1}, {1, 3, 1}, {2, 3, 1}};
    const double tri[][3]  = {{0, 1, 1}, {1, 2, 1}, {0, 2, 1}};
    const double path[][3] = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}};
    /* every pairing of representations, since the target is converted when they differ */
    for (int a = 0; a < 2; a++) {
        for (int b = 0; b < 2; b++) {
            Graph *pattern = make(REPS[a], 3, false, tri, 3);
            Graph *target  = make(REPS[b], 4, false, k4, 6);
            bool found = false;
            int  mapping[3];
            OK(graph_adv_subgraph_isomorphism(pattern, target, true, false, mapping, &found));
            CHECK(found);
            CHECK(graph_representation(target) == REPS[b]);   /* the target is left as it was */
            for (int i = 0; i < 3; i++) {
                for (int j = i + 1; j < 3; j++) {
                    CHECK(graph_has_edge(target, mapping[i], mapping[j], NULL));
                }
            }
            graph_destroy(target);

            target = make(REPS[b], 4, false, path, 3);
            OK(graph_adv_subgraph_isomorphism(pattern, target, false, false, NULL, &found));
            CHECK(!found);   /* a path has no triangle */
            graph_destroy(target);
            graph_destroy(pattern);
        }
    }
    CHECK(graph_adv_subgraph_isomorphism(NULL, NULL, false, false, NULL, NULL) ==
          GRAPH_ERR_INVALID_ARGUMENT);
}

static void test_dynconn(void)
{
    /* two components {0,1} and {2,3}; the tracker must follow edits and also a
     * representation change, which frees the handle it was built on */
    for (int r = 0; r < 2; r++) {
        Graph *g = graph_create_sized(4, false);
        CHECK(graph_add_edge(g, 0, 1, 1) && graph_add_edge(g, 2, 3, 1));
        CHECK(graph_convert(g, REPS[r]));
        GraphAdvDynConn *dc = graph_adv_dynconn_create(g);
        CHECK(dc != NULL);
        bool   c = true;
        size_t n = 0;
        OK(graph_adv_dynconn_connected(dc, 0, 1, &c));
        CHECK(c);
        OK(graph_adv_dynconn_connected(dc, 1, 2, &c));
        CHECK(!c);
        OK(graph_adv_dynconn_components(dc, &n));
        CHECK(n == 2);

        CHECK(graph_add_edge(g, 1, 2, 1));            /* edit */
        OK(graph_adv_dynconn_connected(dc, 0, 3, &c));
        CHECK(c);
        OK(graph_adv_dynconn_components(dc, &n));
        CHECK(n == 1);

        CHECK(graph_convert(g, REPS[1 - r]));         /* old handle freed here */
        OK(graph_adv_dynconn_connected(dc, 0, 3, &c));
        CHECK(c);
        CHECK(graph_remove_edge(g, 1, 2));
        OK(graph_adv_dynconn_connected(dc, 0, 3, &c));
        CHECK(!c);
        OK(graph_adv_dynconn_components(dc, &n));
        CHECK(n == 2);

        CHECK(graph_adv_dynconn_connected(dc, 0, 99, &c) == GRAPH_ERR_INVALID_VERTEX);
        graph_adv_dynconn_destroy(dc);
        graph_destroy(g);
    }
    CHECK(graph_adv_dynconn_create(NULL) == NULL);
    CHECK(graph_adv_dynconn_connected(NULL, 0, 1, NULL) == GRAPH_ERR_INVALID_ARGUMENT);
    graph_adv_dynconn_destroy(NULL);
}


static void test_johnson_eigenvector_isomorphism(void)
{
    /* Johnson with a negative edge: 0->1 (4), 0->2 (5), 2->1 (-3), 1->3 (2) */
    const double e[][3] = {{0, 1, 4}, {0, 2, 5}, {2, 1, -3}, {1, 3, 2}};
    for (int r = 0; r < 2; r++) {
        Graph *g = make(REPS[r], 4, true, e, 4);
        double dist[16], old_dist[16];
        int    pred[16], next[16];
        OK(graph_adv_johnson(g, dist, pred));
        OK(graph_johnson(g, old_dist, next));
        for (int i = 0; i < 16; i++) {
            CHECK(near(dist[i], old_dist[i]));   /* same distances as the existing API */
        }
        CHECK(near(dist[0 * 4 + 1], 2) && near(dist[0 * 4 + 3], 4));
        GraphPath p = {0};
        OK(graph_adv_apsp_path(dist, pred, 4, 0, 3, &p));
        CHECK(p.length == 4 && p.vertices[0] == 0 && p.vertices[1] == 2 &&
              p.vertices[2] == 1 && p.vertices[3] == 3 && near(p.cost, 4));
        graph_path_free(&p);
        CHECK(graph_adv_apsp_path(dist, pred, 4, 3, 0, &p) == GRAPH_ERR_NOT_FOUND);
        graph_destroy(g);

        /* a negative cycle is reported */
        const double nc[][3] = {{0, 1, 1}, {1, 2, -3}, {2, 0, 1}};
        g = make(REPS[r], 3, true, nc, 3);
        CHECK(graph_adv_johnson(g, dist, pred) == GRAPH_ERR_NEGATIVE_CYCLE);
        graph_destroy(g);
    }

    /* eigenvector centrality on a star: the hub scores highest, leaves tie; and on the
     * bipartite path 0-1-2 (where plain power iteration oscillates) it converges */
    const double star[][3] = {{0, 1, 1}, {0, 2, 1}, {0, 3, 1}, {0, 4, 1}};
    const double p3[][3]   = {{0, 1, 1}, {1, 2, 1}};
    for (int r = 0; r < 2; r++) {
        Graph *g = make(REPS[r], 5, false, star, 4);
        double c[5];
        size_t iters = 0;
        OK(graph_adv_eigenvector_centrality(g, 1000, 1e-12, c, &iters));
        CHECK(c[0] > c[1] && near(c[1], c[2]) && near(c[3], c[4]));
        graph_destroy(g);
        g = make(REPS[r], 3, false, p3, 2);
        OK(graph_adv_eigenvector_centrality(g, 1000, 1e-12, c, &iters));
        CHECK(c[1] > c[0] && near(c[0], c[2]));
        graph_destroy(g);
    }

    /* isomorphism: a 5-cycle against a relabelled 5-cycle, a path, and different weights */
    const double c5[][3]  = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 4, 1}, {4, 0, 1}};
    const double c5r[][3] = {{3, 0, 2}, {0, 4, 2}, {4, 1, 2}, {1, 2, 2}, {2, 3, 2}};   /* weights 2 */
    const double p5[][3]  = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 4, 1}};
    for (int a = 0; a < 2; a++) {
        for (int b = 0; b < 2; b++) {
            Graph *x = make(REPS[a], 5, false, c5, 5);
            Graph *y = make(REPS[b], 5, false, c5r, 5);
            bool iso = false;
            int  map[5];
            OK(graph_adv_graph_isomorphism(x, y, false, &iso, map));
            CHECK(iso);
            for (int i = 0; i < 5; i++) {                       /* the mapping is an isomorphism */
                CHECK(graph_has_edge(y, map[i], map[(i + 1) % 5], NULL));
            }
            CHECK(graph_representation(y) == REPS[b]);          /* y is left in its representation */
            OK(graph_adv_graph_isomorphism(x, y, true, &iso, NULL));
            CHECK(!iso);                                        /* weights 1 vs 2 */
            bool old = true;
            OK(graph_is_isomorphic(x, y, &old, NULL));
            CHECK(old);                                         /* the existing check ignores weights */
            graph_destroy(y);
            y = make(REPS[b], 5, false, p5, 4);
            OK(graph_adv_graph_isomorphism(x, y, false, &iso, NULL));
            CHECK(!iso);
            graph_destroy(y);
            graph_destroy(x);
        }
    }
    CHECK(graph_adv_graph_isomorphism(NULL, NULL, false, NULL, NULL) == GRAPH_ERR_INVALID_ARGUMENT);
}

/* The wrappers feed the observer like every other algorithm call. */
static void test_adaptive_counts(void)
{
    Graph *g = graph_create_sized(4, true);
    CHECK(graph_add_edge(g, 0, 1, 1) && graph_add_edge(g, 1, 2, 1));
    double dist[4];
    int    parent[4];
    OK(graph_adv_longest_path(g, 0, dist, parent));
    CHECK(near(dist[2], 2));
    graph_destroy(g);
}

int main(void)
{
    test_dag();
    test_cycles();
    test_trees_and_mst();
    test_structure();
    test_paths();
    test_dynamics();
    test_subgraph_isomorphism();
    test_dynconn();
    test_johnson_eigenvector_isomorphism();
    test_adaptive_counts();
    printf("All graph advanced-wrapper tests passed.\n");
    return 0;
}
