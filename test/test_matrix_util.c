/*
 * test_matrix_util.c — the internal containers in src/matrix_util.h: min-heap,
 * union-find, bitsets, hop-distance BFS and the symmetrised adjacency (CSR).
 */

#include "matrix_util.h"
#include "test_util_matrix.h"

static void test_heap(void)
{
    MuHeap h;
    CHECK(mu_heap_init(&h, 2));                        /* tiny start forces growth */
    double expect[500];
    for (int i = 0; i < 500; ++i) {
        expect[i] = (double)t_rand(1000) - 400.0;     /* includes negatives and duplicates */
        CHECK(mu_heap_push(&h, expect[i], i));
    }
    CHECK(h.size == 500);
    for (int i = 0; i < 500; ++i) {                    /* selection sort the oracle */
        int m = i;
        for (int j = i + 1; j < 500; ++j) {
            if (expect[j] < expect[m]) m = j;
        }
        const double tmp = expect[i]; expect[i] = expect[m]; expect[m] = tmp;
    }
    for (int i = 0; i < 500; ++i) {
        double key; int val;
        mu_heap_pop(&h, &key, &val);
        CHECK(key == expect[i]);
    }
    CHECK(h.size == 0);
    mu_heap_free(&h);
}

static void test_union_find(void)
{
    MuUF uf;
    CHECK(mu_uf_init(&uf, 60));
    int label[60];
    for (int i = 0; i < 60; ++i) label[i] = i;          /* naive labelling oracle */
    for (int step = 0; step < 200; ++step) {
        const int a = (int)t_rand(60), b = (int)t_rand(60);
        const bool joined = mu_uf_union(&uf, a, b);
        CHECK(joined == (label[a] != label[b]));
        if (label[a] != label[b]) {
            const int from = label[b];
            for (int i = 0; i < 60; ++i) if (label[i] == from) label[i] = label[a];
        }
        for (int i = 0; i < 60; i += 7) {
            for (int j = 0; j < 60; j += 5) {
                CHECK((mu_uf_find(&uf, i) == mu_uf_find(&uf, j)) == (label[i] == label[j]));
            }
        }
    }
    mu_uf_free(&uf);
}

static void test_bitset(void)
{
    uint64_t b[mu_words(200)];
    memset(b, 0, sizeof(b));
    bool ref[200] = { false };
    size_t count = 0;
    for (int step = 0; step < 1000; ++step) {
        const size_t i = t_rand(200);
        if (t_rand(2)) { mu_bit_set(b, i); ref[i] = true; } else { mu_bit_clear(b, i); ref[i] = false; }
    }
    for (size_t i = 0; i < 200; ++i) {
        CHECK(mu_bit_test(b, i) == ref[i]);
        count += ref[i];
    }
    CHECK(mu_bit_count(b, mu_words(200)) == count);
    CHECK(mu_words(0) == 0 && mu_words(1) == 1 && mu_words(64) == 1 && mu_words(65) == 2);
    CHECK(mu_mul(3, 4) == 12 && mu_mul(SIZE_MAX, 2) == 0 && mu_mul(0, 5) == 0);
}

static void test_hops(void)
{
    for (int seed = 1; seed <= 20; ++seed) {
        MatrixGraph *g = random_graph(15, 0.15, seed % 2, 3, (unsigned)seed, (size_t)(seed % 3));
        int D[225], dist[15], queue[15];
        hops(g, D);
        for (int s = 0; s < 15; ++s) {
            if (!matrix_graph_is_alive(g, s)) continue;
            const size_t reached = mu_bfs_hops(g, s, dist, queue);
            size_t want = 0;
            for (int v = 0; v < 15; ++v) {
                const bool r = D[s * 15 + v] < HOP_INF;
                CHECK(dist[v] == (r ? D[s * 15 + v] : -1));
                want += r;
            }
            CHECK(reached == want && queue[0] == s);
        }
        matrix_graph_destroy(g);
    }
}

/** The symmetrised adjacency against a straightforward dense construction. */
static void test_sym(void)
{
    for (int seed = 1; seed <= 60; ++seed) {
        for (int directed = 0; directed < 2; ++directed) {
            MatrixGraph *g = random_graph(12, 0.25, directed, 9, (unsigned)seed, (size_t)(seed % 3));
            if (seed % 5 == 0 && matrix_graph_is_alive(g, 2)) {
                (void)matrix_graph_add_edge(g, 2, 2, 4.0);        /* loops are dropped */
            }
            MuSym s;
            CHECK(mu_sym_build(g, &s));
            CHECK(s.n == 12);
            for (int u = 0; u < 12; ++u) {
                int prev = -1;
                size_t entries = 0;
                for (size_t k = s.off[u]; k < s.off[u + 1]; ++k) {
                    const int v = s.to[k];
                    CHECK(v > prev && v != u);                      /* sorted, unique, no loops */
                    prev = v;
                    ++entries;
                    const double want = directed ? wt(g, u, v) + wt(g, v, u) : wt(g, u, v);
                    CHECK(want > 0.0 && s.w[k] == want);
                    /* Symmetric: the mirrored entry exists with the same weight. */
                    bool mirrored = false;
                    for (size_t j = s.off[v]; j < s.off[v + 1]; ++j) {
                        mirrored = mirrored || (s.to[j] == u && s.w[j] == s.w[k]);
                    }
                    CHECK(mirrored);
                }
                CHECK(entries == mu_sym_degree(&s, u));
                size_t expect = 0;
                for (int v = 0; v < 12; ++v) {
                    expect += v != u && (matrix_graph_has_edge(g, u, v, NULL) || matrix_graph_has_edge(g, v, u, NULL));
                }
                CHECK(entries == expect);
            }
            mu_sym_free(&s);
            matrix_graph_destroy(g);
        }
    }
    MatrixGraph *empty = matrix_graph_create(0, true);
    MuSym s;
    CHECK(mu_sym_build(empty, &s) && s.n == 0 && s.off[0] == 0);
    mu_sym_free(&s);
    matrix_graph_destroy(empty);
}

/* mu_sym_build_ex: weights are optional, and rows longer than the insertion-sort limit (hubs) still come out
 * sorted with every weight attached to its own neighbour. */
static void test_sym_options(void)
{
    for (int directed = 0; directed < 2; ++directed) {
        const int n = 400;
        MatrixGraph *g = matrix_graph_create((size_t)n, directed != 0);
        CHECK(g != NULL);
        for (int v = 1; v < n; ++v) {                        /* hub 0: ~399 neighbours, added in a scrambled order */
            const int a = (v * 7919) % (n - 1) + 1;
            CHECK(matrix_graph_add_edge(g, directed && (v & 1) ? a : 0, directed && (v & 1) ? 0 : a, 1.0 + (double)(a % 9)) ||
                  matrix_graph_has_edge(g, 0, a, NULL) || matrix_graph_has_edge(g, a, 0, NULL));
        }
        for (int v = 1; v + 1 < n; v += 3) {
            (void)matrix_graph_add_edge(g, v, v + 1, 2.5);
        }
        MuSym with, without;
        CHECK(mu_sym_build_ex(g, &with, true));
        CHECK(mu_sym_build_ex(g, &without, false));
        CHECK(with.w != NULL && without.w == NULL);
        CHECK(with.n == without.n && with.off[with.n] == without.off[without.n]);
        CHECK(mu_sym_degree(&with, 0) > 32);                 /* the long-row path was exercised */
        for (size_t k = 0; k < with.off[with.n]; ++k) {
            CHECK(with.to[k] == without.to[k]);
        }
        for (int u = 0; u < n; ++u) {
            CHECK(with.off[u] == without.off[u]);
            for (size_t k = with.off[u]; k < with.off[u + 1]; ++k) {
                CHECK(k == with.off[u] || with.to[k - 1] < with.to[k]);
                const int v = with.to[k];
                const double want = directed ? wt(g, u, v) + wt(g, v, u) : wt(g, u, v);
                CHECK(want > 0.0 && with.w[k] == want);      /* the weight stayed with its neighbour */
            }
        }
        mu_sym_free(&with);
        mu_sym_free(&without);
        matrix_graph_destroy(g);
    }
}

int main(void)
{
    t_seed(11);
    test_heap();
    test_union_find();
    test_bitset();
    test_hops();
    test_sym();
    test_sym_options();
    printf("All matrix utility tests passed.\n");
    return 0;
}
