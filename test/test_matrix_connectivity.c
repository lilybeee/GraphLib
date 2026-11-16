/*
 * test_matrix_connectivity.c — SCC, cut vertices, bridges, blocks, cycle
 * extraction, closure / reduction and the structural predicates.
 */

#include "test_util_matrix.h"

static void test_predicates(void)
{
    MatrixGraph *path = matrix_gen_path(5, false), *cyc = matrix_gen_cycle(5, false);
    MatrixGraph *dag = matrix_gen_path(5, true), *dcyc = matrix_gen_cycle(5, true);
    MatrixGraph *empty = matrix_graph_create(0, false), *iso = matrix_graph_create(3, false);

    CHECK(matrix_is_connected(path) && matrix_is_connected(cyc) && matrix_is_connected(empty));
    CHECK(!matrix_is_connected(iso));
    CHECK(matrix_is_tree(path) && !matrix_is_tree(cyc) && !matrix_is_tree(iso) && !matrix_is_tree(empty));
    CHECK(matrix_is_forest(path) && matrix_is_forest(iso) && !matrix_is_forest(cyc));
    CHECK(!matrix_is_tree(dag) && !matrix_is_forest(dag));            /* undirected notions */
    CHECK(matrix_is_dag(dag) && !matrix_is_dag(dcyc) && !matrix_is_dag(path));
    CHECK(matrix_is_strongly_connected(dcyc) && !matrix_is_strongly_connected(dag));
    CHECK(matrix_is_connected(dag));                                   /* weakly */
    CHECK(matrix_is_strongly_connected(path) && !matrix_is_strongly_connected(iso));
    CHECK(!matrix_is_connected(NULL) && !matrix_is_dag(NULL));

    MatrixGraph *loop = matrix_graph_create(2, false);
    CHECK(matrix_graph_add_edge(loop, 0, 1, 1) && matrix_graph_add_edge(loop, 0, 0, 1));
    CHECK(!matrix_is_tree(loop) && !matrix_is_forest(loop));          /* a loop is a cycle */

    matrix_graph_destroy(path); matrix_graph_destroy(cyc); matrix_graph_destroy(dag);
    matrix_graph_destroy(dcyc); matrix_graph_destroy(empty); matrix_graph_destroy(iso);
    matrix_graph_destroy(loop);
}

/* ---- strongly connected components ---- */

static void test_scc(void)
{
    for (int seed = 1; seed <= 60; ++seed) {
        MatrixGraph *g = random_graph(16, 0.04 + 0.01 * (seed % 8), true, 3, (unsigned)seed, (size_t)(seed % 3));
        const size_t n = matrix_graph_id_bound(g);
        int comp[16], D[256];
        size_t count = 0;
        hops(g, D);
        CHECK_OK(matrix_strongly_connected_components(g, comp, &count));
        size_t distinct = 0;
        bool seen[16] = { false };
        for (size_t u = 0; u < n; ++u) {
            if (!matrix_graph_is_alive(g, (int)u)) {
                CHECK(comp[u] == -1);
                continue;
            }
            CHECK(comp[u] >= 0 && (size_t)comp[u] < count);
            if (!seen[comp[u]]) {
                seen[comp[u]] = true;
                ++distinct;
            }
            for (size_t v = 0; v < n; ++v) {
                if (!matrix_graph_is_alive(g, (int)v)) {
                    continue;
                }
                const bool mutual = D[u * n + v] < HOP_INF && D[v * n + u] < HOP_INF;
                CHECK((comp[u] == comp[v]) == mutual);
                /* Reverse topological numbering: arcs never go to a larger component. */
                if (matrix_graph_has_edge(g, (int)u, (int)v, NULL)) {
                    CHECK(comp[u] >= comp[v]);
                }
            }
        }
        CHECK(distinct == count);
        matrix_graph_destroy(g);
    }
    MatrixGraph *u = random_graph(12, 0.15, false, 2, 4, 1);          /* undirected: ordinary components */
    int comp[12];
    size_t count;
    CHECK_OK(matrix_strongly_connected_components(u, comp, &count));
    CHECK(count == naive_components(u));
    matrix_graph_destroy(u);
    CHECK(matrix_strongly_connected_components(NULL, comp, &count) == GRAPH_ERR_INVALID_ARGUMENT);
}

/* ---- articulation points, bridges, blocks ---- */

static MatrixGraph *drop_vertex(const MatrixGraph *g, size_t v)
{
    bool drop[64] = { false };
    int map[64];
    drop[v] = true;
    return survivors(g, drop, map);
}

static MatrixGraph *drop_edge(const MatrixGraph *g, int a, int b)
{
    MatrixGraph *r = matrix_graph_copy(g);
    CHECK(r != NULL && matrix_graph_remove_edge(r, a, b));
    return r;
}

static void test_cut_structure(void)
{
    for (int seed = 1; seed <= 80; ++seed) {
        MatrixGraph *g = random_graph(14, 0.08 + 0.01 * (seed % 10), false, 3, (unsigned)seed, (size_t)(seed % 3));
        const size_t n = matrix_graph_id_bound(g);
        const size_t base = naive_components(g);
        bool cut[14];
        size_t cuts = 0;
        CHECK_OK(matrix_articulation_points(g, cut, &cuts));
        size_t expect_cuts = 0;
        for (size_t v = 0; v < n; ++v) {
            bool expect = false;
            if (matrix_graph_is_alive(g, (int)v)) {
                MatrixGraph *h = drop_vertex(g, v);
                expect = naive_components(h) > base;
                matrix_graph_destroy(h);
            }
            CHECK(cut[v] == expect);
            expect_cuts += expect;
        }
        CHECK(cuts == expect_cuts);

        MatrixEdge br[200];
        size_t nb = 0;
        CHECK_OK(matrix_bridges(g, br, &nb));
        size_t expect_bridges = 0;
        for (size_t u = 0; u < n; ++u) {
            for (size_t v = u + 1; v < n; ++v) {
                if (!matrix_graph_has_edge(g, (int)u, (int)v, NULL)) {
                    continue;
                }
                MatrixGraph *h = drop_edge(g, (int)u, (int)v);
                const bool is_bridge = naive_components(h) > base;
                matrix_graph_destroy(h);
                expect_bridges += is_bridge;
                bool listed = false;
                for (size_t i = 0; i < nb; ++i) {
                    listed = listed || (br[i].u == (int)u && br[i].v == (int)v);
                }
                CHECK(listed == is_bridge);
            }
        }
        CHECK(nb == expect_bridges);
        for (size_t i = 0; i < nb; ++i) {
            CHECK(br[i].u < br[i].v && br[i].weight == wt(g, br[i].u, br[i].v));
        }
        matrix_graph_destroy(g);
    }

    MatrixGraph *d = matrix_gen_path(3, true);
    bool c[3];
    MatrixEdge e[3];
    size_t k;
    CHECK(matrix_articulation_points(d, c, &k) == GRAPH_ERR_UNSUPPORTED);
    CHECK(matrix_bridges(d, e, &k) == GRAPH_ERR_UNSUPPORTED);
    MatrixGroups gr;
    CHECK(matrix_biconnected_components(d, &gr) == GRAPH_ERR_UNSUPPORTED);
    matrix_graph_destroy(d);
}

static int group_size(const MatrixGroups *g, size_t i)
{
    return (int)(g->offsets[i + 1] - g->offsets[i]);
}

static void test_blocks(void)
{
    /* Two triangles sharing vertex 2: two blocks of size 3. */
    const double bowtie[] = { 0,1,1, 1,2,1, 0,2,1, 2,3,1, 3,4,1, 2,4,1 };
    MatrixGraph *b = build(5, false, bowtie, 6);
    MatrixGroups gr;
    CHECK_OK(matrix_biconnected_components(b, &gr));
    CHECK(gr.count == 2 && group_size(&gr, 0) == 3 && group_size(&gr, 1) == 3);
    matrix_groups_free(&gr);
    matrix_graph_destroy(b);

    MatrixGraph *path = matrix_gen_path(5, false);          /* four blocks of size 2 */
    CHECK_OK(matrix_biconnected_components(path, &gr));
    CHECK(gr.count == 4);
    for (size_t i = 0; i < gr.count; ++i) {
        CHECK(group_size(&gr, i) == 2);
    }
    matrix_groups_free(&gr);
    matrix_graph_destroy(path);

    MatrixGraph *cyc = matrix_gen_cycle(6, false);
    CHECK_OK(matrix_biconnected_components(cyc, &gr));
    CHECK(gr.count == 1 && group_size(&gr, 0) == 6);
    matrix_groups_free(&gr);
    matrix_graph_destroy(cyc);

    MatrixGraph *none = matrix_graph_create(3, false);       /* isolated vertices have no block */
    CHECK_OK(matrix_biconnected_components(none, &gr));
    CHECK(gr.count == 0);
    matrix_groups_free(&gr);
    matrix_graph_destroy(none);

    for (int seed = 1; seed <= 60; ++seed) {
        MatrixGraph *g = random_graph(13, 0.1 + 0.015 * (seed % 8), false, 3, (unsigned)seed, (size_t)(seed % 3));
        const size_t n = matrix_graph_id_bound(g);
        CHECK_OK(matrix_biconnected_components(g, &gr));
        /* Every edge lies in exactly one block. */
        for (size_t u = 0; u < n; ++u) {
            for (size_t v = u + 1; v < n; ++v) {
                if (!matrix_graph_has_edge(g, (int)u, (int)v, NULL)) {
                    continue;
                }
                int holders = 0;
                for (size_t i = 0; i < gr.count; ++i) {
                    bool hu = false, hv = false;
                    for (size_t k = gr.offsets[i]; k < gr.offsets[i + 1]; ++k) {
                        hu = hu || gr.members[k] == (int)u;
                        hv = hv || gr.members[k] == (int)v;
                    }
                    holders += hu && hv;
                }
                CHECK(holders == 1);
            }
        }
        /* Each block is 2-connected (or a single edge) and blocks overlap in <= 1 vertex. */
        for (size_t i = 0; i < gr.count; ++i) {
            const size_t sz = (size_t)group_size(&gr, i);
            CHECK(sz >= 2);
            MatrixGraph *sub = NULL;
            CHECK_OK(matrix_induced_subgraph(g, gr.members + gr.offsets[i], sz, &sub));
            CHECK(matrix_is_connected(sub));
            if (sz >= 3) {
                bool cut[16];
                size_t k;
                CHECK_OK(matrix_articulation_points(sub, cut, &k));
                CHECK(k == 0);
            }
            matrix_graph_destroy(sub);
            for (size_t j = i + 1; j < gr.count; ++j) {
                int shared = 0;
                for (size_t a = gr.offsets[i]; a < gr.offsets[i + 1]; ++a) {
                    for (size_t c = gr.offsets[j]; c < gr.offsets[j + 1]; ++c) {
                        shared += gr.members[a] == gr.members[c];
                    }
                }
                CHECK(shared <= 1);
            }
        }
        matrix_groups_free(&gr);
        matrix_graph_destroy(g);
    }
}

/* ---- cycle extraction ---- */

static void check_cycle(const MatrixGraph *g, const int *cyc, size_t len)
{
    CHECK(len >= 1);
    if (!matrix_graph_is_directed(g) && len > 1) {
        CHECK(len >= 3);
    }
    for (size_t i = 0; i < len; ++i) {
        CHECK(matrix_graph_has_edge(g, cyc[i], cyc[(i + 1) % len], NULL));
        for (size_t j = i + 1; j < len; ++j) {
            CHECK(cyc[i] != cyc[j]);
        }
    }
}

static void test_find_cycle(void)
{
    for (int seed = 1; seed <= 100; ++seed) {
        for (int directed = 0; directed < 2; ++directed) {
            MatrixGraph *g = random_graph(12, 0.05 + 0.01 * (seed % 9), directed, 3, (unsigned)seed,
                                          (size_t)(seed % 3));
            int cyc[12];
            const size_t len = matrix_find_cycle(g, cyc);
            CHECK((len > 0) == matrix_has_cycle(g));
            if (len > 0) {
                check_cycle(g, cyc, len);
            }
            matrix_graph_destroy(g);
        }
    }
    int cyc[8];
    MatrixGraph *loop = matrix_graph_create(4, false);
    CHECK(matrix_graph_add_edge(loop, 2, 2, 1.0));
    CHECK(matrix_find_cycle(loop, cyc) == 1 && cyc[0] == 2);
    matrix_graph_destroy(loop);
    MatrixGraph *dloop = matrix_graph_create(3, true);
    CHECK(matrix_graph_add_edge(dloop, 1, 1, 1.0));
    CHECK(matrix_find_cycle(dloop, cyc) == 1 && cyc[0] == 1);
    matrix_graph_destroy(dloop);
    MatrixGraph *tree = matrix_gen_random_tree(8, 3);
    CHECK(matrix_find_cycle(tree, cyc) == 0);
    matrix_graph_destroy(tree);
    CHECK(matrix_find_cycle(NULL, cyc) == 0);
    /* A 2-vertex undirected "cycle" would just reuse the one edge: not a cycle. */
    MatrixGraph *edge = matrix_gen_path(2, false);
    CHECK(matrix_find_cycle(edge, cyc) == 0);
    matrix_graph_destroy(edge);
}

/* ---- closure and reduction ---- */

static void test_closure_reduction(void)
{
    for (int seed = 1; seed <= 40; ++seed) {
        for (int directed = 0; directed < 2; ++directed) {
            MatrixGraph *g = random_graph(11, 0.12, directed, 4, (unsigned)seed, (size_t)(seed % 3));
            MatrixGraph *c = NULL;
            const size_t n = matrix_graph_id_bound(g);
            int D[121];
            hops(g, D);
            CHECK_OK(matrix_transitive_closure(g, &c));
            CHECK(matrix_graph_is_directed(c) == (bool)directed);
            for (size_t u = 0; u < n; ++u) {
                CHECK(matrix_graph_is_alive(c, (int)u) == matrix_graph_is_alive(g, (int)u));
                for (size_t v = 0; v < n; ++v) {
                    const bool want = u != v && D[u * n + v] < HOP_INF;
                    CHECK(matrix_graph_has_edge(c, (int)u, (int)v, NULL) == want);
                }
            }
            matrix_graph_destroy(c);
            matrix_graph_destroy(g);
        }
    }
    for (int seed = 1; seed <= 40; ++seed) {
        MatrixGraph *dag = matrix_gen_random_dag(10, 0.35, (unsigned)seed, 1.0, 5.0), *red = NULL;
        const size_t n = matrix_graph_id_bound(dag);
        int Dg[100], Dr[100];
        hops(dag, Dg);
        CHECK_OK(matrix_transitive_reduction(dag, &red));
        hops(red, Dr);
        for (size_t i = 0; i < n * n; ++i) {
            CHECK((Dg[i] < HOP_INF) == (Dr[i] < HOP_INF));   /* same reachability */
        }
        for (size_t u = 0; u < n; ++u) {
            for (size_t v = 0; v < n; ++v) {
                if (!matrix_graph_has_edge(red, (int)u, (int)v, NULL)) {
                    continue;
                }
                CHECK(wt(red, (int)u, (int)v) == wt(dag, (int)u, (int)v));
                /* Minimal: dropping the edge must lose the reachability u -> v. */
                MatrixGraph *h = drop_edge(red, (int)u, (int)v);
                int Dh[100];
                hops(h, Dh);
                CHECK(Dh[u * n + v] >= HOP_INF);
                matrix_graph_destroy(h);
            }
        }
        matrix_graph_destroy(red);
        matrix_graph_destroy(dag);
    }
    MatrixGraph *cyc = matrix_gen_cycle(4, true), *out = NULL;
    CHECK(matrix_transitive_reduction(cyc, &out) == GRAPH_ERR_UNSUPPORTED);
    matrix_graph_destroy(cyc);
    MatrixGraph *und = matrix_gen_path(4, false);
    CHECK(matrix_transitive_reduction(und, &out) == GRAPH_ERR_UNSUPPORTED);
    matrix_graph_destroy(und);
}

int main(void)
{
    test_predicates();
    test_scc();
    test_cut_structure();
    test_blocks();
    test_find_cycle();
    test_closure_reduction();
    printf("All matrix connectivity tests passed.\n");
    return 0;
}
