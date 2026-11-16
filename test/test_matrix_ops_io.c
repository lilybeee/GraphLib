/*
 * test_matrix_ops_io.c — copy, subgraph, union / intersection / complement,
 * isomorphism, generators, and every input / output format.
 */

#include "test_util_matrix.h"

/* ---- copy, subgraph, set operations ---- */

static void test_copy_subgraph(void)
{
    for (int seed = 1; seed <= 40; ++seed) {
        for (int directed = 0; directed < 2; ++directed) {
            MatrixGraph *g = random_graph(12, 0.25, directed, 9, (unsigned)seed, (size_t)(seed % 4));
            MatrixGraph *c = matrix_graph_copy(g);
            CHECK(c != NULL && c != g && same_graph(g, c));
            CHECK(matrix_graph_add_vertex(c) >= 0);                  /* the copy is independent */
            CHECK(!same_graph(g, c));
            matrix_graph_destroy(c);

            /* Induced subgraph on a shuffled subset of alive vertices. */
            int pick[12], k = 0;
            for (int v = 11; v >= 0; --v) {
                if (matrix_graph_is_alive(g, v) && (v + seed) % 3 != 0) {
                    pick[k++] = v;
                }
            }
            MatrixGraph *sub = NULL;
            CHECK_OK(matrix_induced_subgraph(g, pick, (size_t)k, &sub));
            CHECK(matrix_graph_vertex_count(sub) == (size_t)k && matrix_graph_id_bound(sub) == (size_t)k);
            size_t edges = 0;
            for (int i = 0; i < k; ++i) {
                for (int j = 0; j < k; ++j) {
                    const bool want = matrix_graph_has_edge(g, pick[i], pick[j], NULL);
                    CHECK(matrix_graph_has_edge(sub, i, j, NULL) == want);
                    CHECK(wt(sub, i, j) == wt(g, pick[i], pick[j]));
                    edges += want && (directed || i <= j);
                }
            }
            CHECK(matrix_graph_edge_count(sub) == edges);
            matrix_graph_destroy(sub);
            matrix_graph_destroy(g);
        }
    }
    MatrixGraph *g = matrix_gen_path(5, false), *out = NULL;
    const int dup[] = { 1, 1 }, dead[] = { 9 }, ok[] = { 2, 0 };
    CHECK(matrix_induced_subgraph(g, dup, 2, &out) == GRAPH_ERR_INVALID_ARGUMENT);
    CHECK(matrix_induced_subgraph(g, dead, 1, &out) == GRAPH_ERR_INVALID_ARGUMENT);
    CHECK_OK(matrix_induced_subgraph(g, ok, 2, &out));
    CHECK(matrix_graph_edge_count(out) == 0);
    matrix_graph_destroy(out);
    CHECK_OK(matrix_induced_subgraph(g, NULL, 0, &out));
    CHECK(matrix_graph_vertex_count(out) == 0);
    matrix_graph_destroy(out);
    CHECK(matrix_graph_copy(NULL) == NULL);
    matrix_graph_destroy(g);
}

static void test_set_operations(void)
{
    for (int seed = 1; seed <= 40; ++seed) {
        for (int directed = 0; directed < 2; ++directed) {
            MatrixGraph *a = random_graph(10, 0.25, directed, 9, (unsigned)seed, (size_t)(seed % 3));
            MatrixGraph *b = random_graph(12, 0.25, directed, 9, (unsigned)seed + 500u, (size_t)(seed % 2));
            MatrixGraph *u = NULL, *x = NULL, *comp = NULL;
            CHECK_OK(matrix_graph_union(a, b, &u));
            CHECK_OK(matrix_graph_intersection(a, b, &x));
            CHECK(matrix_graph_id_bound(u) == 12 && matrix_graph_id_bound(x) == 12);
            for (int s = 0; s < 12; ++s) {
                const bool in_a = s < 10 && matrix_graph_is_alive(a, s), in_b = matrix_graph_is_alive(b, s);
                CHECK(matrix_graph_is_alive(u, s) == (in_a || in_b));
                CHECK(matrix_graph_is_alive(x, s) == (in_a && in_b));
                for (int t = 0; t < 12; ++t) {
                    const bool ea = s < 10 && t < 10 && matrix_graph_has_edge(a, s, t, NULL);
                    const bool eb = matrix_graph_has_edge(b, s, t, NULL);
                    CHECK(matrix_graph_has_edge(u, s, t, NULL) == (ea || eb));
                    CHECK(matrix_graph_has_edge(x, s, t, NULL) == (ea && eb));
                    if (ea) {
                        CHECK(wt(u, s, t) == wt(a, s, t));            /* a's weight wins */
                    } else if (eb) {
                        CHECK(wt(u, s, t) == wt(b, s, t));
                    }
                    if (ea && eb) {
                        CHECK(wt(x, s, t) == wt(a, s, t));
                    }
                }
            }
            CHECK_OK(matrix_graph_complement(a, &comp));
            for (int s = 0; s < 10; ++s) {
                CHECK(matrix_graph_is_alive(comp, s) == matrix_graph_is_alive(a, s));
                for (int t = 0; t < 10; ++t) {
                    const bool both = matrix_graph_is_alive(a, s) && matrix_graph_is_alive(a, t) && s != t;
                    CHECK(matrix_graph_has_edge(comp, s, t, NULL) == (both && !matrix_graph_has_edge(a, s, t, NULL)));
                }
            }
            MatrixGraph *back = NULL;                                /* complement twice: same structure */
            CHECK_OK(matrix_graph_complement(comp, &back));
            for (int s = 0; s < 10; ++s) {
                for (int t = 0; t < 10; ++t) {
                    CHECK(matrix_graph_has_edge(back, s, t, NULL) ==
                          (s != t && matrix_graph_has_edge(a, s, t, NULL)));
                }
            }
            matrix_graph_destroy(a); matrix_graph_destroy(b); matrix_graph_destroy(u);
            matrix_graph_destroy(x); matrix_graph_destroy(comp); matrix_graph_destroy(back);
        }
    }
    MatrixGraph *d = matrix_gen_path(3, true), *un = matrix_gen_path(3, false), *out = NULL;
    CHECK(matrix_graph_union(d, un, &out) == GRAPH_ERR_INVALID_ARGUMENT);
    CHECK(matrix_graph_intersection(d, un, &out) == GRAPH_ERR_INVALID_ARGUMENT);
    matrix_graph_destroy(d);
    matrix_graph_destroy(un);
}

/* ---- isomorphism ---- */

static MatrixGraph *relabel(const MatrixGraph *g, const int *perm)
{
    const size_t n = matrix_graph_id_bound(g);
    MatrixGraph *r = matrix_graph_create(n, matrix_graph_is_directed(g));
    for (size_t u = 0; u < n; ++u) {
        for (size_t v = 0; v < n; ++v) {
            if (matrix_graph_has_edge(g, (int)u, (int)v, NULL) &&
                !matrix_graph_has_edge(r, perm[u], perm[v], NULL)) {
                CHECK(matrix_graph_add_edge(r, perm[u], perm[v], 1.0));
            }
        }
    }
    return r;
}

/** Brute force: try every permutation of the vertices (n <= 8). */
static bool brute_iso_rec(const MatrixGraph *a, const MatrixGraph *b, int n, int *perm, unsigned used, int pos)
{
    if (pos == n) {
        for (int u = 0; u < n; ++u) {
            for (int v = 0; v < n; ++v) {
                if (matrix_graph_has_edge(a, u, v, NULL) != matrix_graph_has_edge(b, perm[u], perm[v], NULL)) {
                    return false;
                }
            }
        }
        return true;
    }
    for (int w = 0; w < n; ++w) {
        if (!((used >> w) & 1u)) {
            perm[pos] = w;
            if (brute_iso_rec(a, b, n, perm, used | (1u << w), pos + 1)) {
                return true;
            }
        }
    }
    return false;
}

static bool brute_iso(const MatrixGraph *a, const MatrixGraph *b)
{
    int perm[8];
    return brute_iso_rec(a, b, (int)matrix_graph_id_bound(a), perm, 0, 0);
}

static void test_isomorphism(void)
{
    /* Random pairs of small graphs against the permutation oracle. */
    int agree_true = 0, agree_false = 0;
    for (int seed = 1; seed <= 300; ++seed) {
        for (int directed = 0; directed < 2; ++directed) {
            MatrixGraph *a = random_graph(7, 0.35, directed, 1, (unsigned)seed, 0);
            MatrixGraph *b = random_graph(7, 0.35, directed, 1, (unsigned)seed + 9000u, 0);
            bool iso = false;
            CHECK_OK(matrix_is_isomorphic(a, b, &iso, NULL));
            const bool want = brute_iso(a, b);
            CHECK(iso == want);
            agree_true += want;
            agree_false += !want;
            matrix_graph_destroy(a);
            matrix_graph_destroy(b);
        }
    }
    CHECK(agree_true > 0 && agree_false > 0);            /* the sample really contains both answers */
    for (int seed = 1; seed <= 60; ++seed) {
        for (int directed = 0; directed < 2; ++directed) {
            MatrixGraph *a = random_graph(9, 0.3, directed, 1, (unsigned)seed, 0);
            int perm[9];
            for (int i = 0; i < 9; ++i) {
                perm[i] = i;
            }
            for (int i = 8; i > 0; --i) {
                const int j = (int)t_rand((unsigned)i + 1);
                const int tmp = perm[i]; perm[i] = perm[j]; perm[j] = tmp;
            }
            MatrixGraph *b = relabel(a, perm);
            bool iso = false;
            int map[9];
            CHECK_OK(matrix_is_isomorphic(a, b, &iso, map));
            CHECK(iso);
            for (int u = 0; u < 9; ++u) {                           /* the mapping is a real isomorphism */
                for (int v = 0; v < 9; ++v) {
                    CHECK(matrix_graph_has_edge(a, u, v, NULL) == matrix_graph_has_edge(b, map[u], map[v], NULL));
                }
            }
            matrix_graph_destroy(a);
            matrix_graph_destroy(b);
        }
    }
    /* Moving one edge keeps the counts; the answer must still match the oracle. */
    for (int seed = 1; seed <= 100; ++seed) {
        for (int directed = 0; directed < 2; ++directed) {
            MatrixGraph *a = random_graph(7, 0.4, directed, 1, (unsigned)seed + 300u, 0);
            MatrixGraph *b = matrix_graph_copy(a);
            int eu = -1, ev = -1;
            for (int u = 0; u < 7 && eu < 0; ++u) {
                for (int v = 0; v < 7 && eu < 0; ++v) {
                    if (u != v && matrix_graph_has_edge(b, u, v, NULL)) { eu = u; ev = v; }
                }
            }
            bool moved = false;
            if (eu >= 0) {
                CHECK(matrix_graph_remove_edge(b, eu, ev));
                for (int u = 0; u < 7 && !moved; ++u) {
                    for (int v = 0; v < 7 && !moved; ++v) {
                        if (u != v && !matrix_graph_has_edge(a, u, v, NULL) &&
                            !(!directed && matrix_graph_has_edge(a, v, u, NULL))) {
                            moved = matrix_graph_add_edge(b, u, v, 1.0);
                        }
                    }
                }
            }
            if (moved) {
                bool iso = true;
                CHECK_OK(matrix_is_isomorphic(a, b, &iso, NULL));
                CHECK(iso == brute_iso(a, b));
            }
            matrix_graph_destroy(a);
            matrix_graph_destroy(b);
        }
    }
    /* Same degree sequences, different structure. */
    MatrixGraph *c6 = matrix_gen_cycle(6, false), *two = matrix_graph_create(6, false);
    for (int i = 0; i < 3; ++i) {
        CHECK(matrix_graph_add_edge(two, i, (i + 1) % 3, 1) && matrix_graph_add_edge(two, 3 + i, 3 + (i + 1) % 3, 1));
    }
    bool iso;
    CHECK_OK(matrix_is_isomorphic(c6, two, &iso, NULL));
    CHECK(!iso);
    matrix_graph_destroy(c6);
    matrix_graph_destroy(two);

    MatrixGraph *k33 = matrix_graph_create(6, false), *prism = matrix_graph_create(6, false);
    for (int u = 0; u < 3; ++u) {
        for (int v = 3; v < 6; ++v) {
            CHECK(matrix_graph_add_edge(k33, u, v, 1));
        }
        CHECK(matrix_graph_add_edge(prism, u, (u + 1) % 3, 1) && matrix_graph_add_edge(prism, 3 + u, 3 + (u + 1) % 3, 1) &&
              matrix_graph_add_edge(prism, u, u + 3, 1));
    }
    CHECK_OK(matrix_is_isomorphic(k33, prism, &iso, NULL));          /* both 3-regular on 6 vertices */
    CHECK(!iso);
    matrix_graph_destroy(k33);
    matrix_graph_destroy(prism);

    MatrixGraph *d = matrix_gen_path(3, true), *u = matrix_gen_path(3, false);
    CHECK_OK(matrix_is_isomorphic(d, u, &iso, NULL));
    CHECK(!iso);                                                      /* directedness differs */
    matrix_graph_destroy(d);
    matrix_graph_destroy(u);

    MatrixGraph *l1 = matrix_graph_create(2, false), *l2 = matrix_graph_create(2, false);
    CHECK(matrix_graph_add_edge(l1, 0, 0, 1) && matrix_graph_add_edge(l1, 0, 1, 1));
    CHECK(matrix_graph_add_edge(l2, 1, 1, 1) && matrix_graph_add_edge(l2, 0, 1, 1));
    int map[2];
    CHECK_OK(matrix_is_isomorphic(l1, l2, &iso, map));               /* loops must map to loops */
    CHECK(iso && map[0] == 1 && map[1] == 0);
    matrix_graph_destroy(l1);
    matrix_graph_destroy(l2);
    CHECK(matrix_is_isomorphic(NULL, NULL, &iso, NULL) == GRAPH_ERR_INVALID_ARGUMENT);
}

/* ---- generators ---- */

static void test_generators(void)
{
    MatrixGraph *a = matrix_gen_erdos_renyi(60, 0.2, 7, false, 1.0, 5.0);
    MatrixGraph *b = matrix_gen_erdos_renyi(60, 0.2, 7, false, 1.0, 5.0);
    MatrixGraph *c = matrix_gen_erdos_renyi(60, 0.2, 8, false, 1.0, 5.0);
    CHECK(a && b && c && same_graph(a, b) && !same_graph(a, c));    /* seeded, reproducible */
    const double pairs = 60.0 * 59.0 / 2.0;
    CHECK(fabs((double)matrix_graph_edge_count(a) / pairs - 0.2) < 0.06);
    for (int u = 0; u < 60; ++u) {
        for (int v = 0; v < 60; ++v) {
            double w;
            if (matrix_graph_has_edge(a, u, v, &w)) {
                CHECK(w >= 1.0 && w <= 5.0 && u != v);
            }
        }
    }
    matrix_graph_destroy(a); matrix_graph_destroy(b); matrix_graph_destroy(c);

    MatrixGraph *empty = matrix_gen_erdos_renyi(10, 0.0, 1, true, 1, 1);
    MatrixGraph *full  = matrix_gen_erdos_renyi(10, 1.0, 1, true, 1, 1);
    CHECK(matrix_graph_edge_count(empty) == 0 && matrix_graph_edge_count(full) == 90);
    matrix_graph_destroy(empty); matrix_graph_destroy(full);
    CHECK(matrix_gen_erdos_renyi(5, 1.5, 1, false, 1, 1) == NULL);
    CHECK(matrix_gen_erdos_renyi(5, 0.5, 1, false, 3, 1) == NULL);
    CHECK(matrix_gen_erdos_renyi(5, 0.5, 1, false, 1, INFINITY) == NULL);

    MatrixGraph *ba = matrix_gen_barabasi_albert(100, 3, 42);
    CHECK(ba != NULL && matrix_graph_vertex_count(ba) == 100 && matrix_is_connected(ba));
    CHECK(matrix_graph_edge_count(ba) == 6 + 96 * 3);                 /* K4 seed + 3 per newcomer */
    size_t minimum = 1000;
    for (int v = 0; v < 100; ++v) {
        if (matrix_graph_degree(ba, v) < minimum) {
            minimum = matrix_graph_degree(ba, v);
        }
    }
    CHECK(minimum >= 3);
    MatrixGraph *ba2 = matrix_gen_barabasi_albert(100, 3, 42);
    CHECK(same_graph(ba, ba2));
    matrix_graph_destroy(ba); matrix_graph_destroy(ba2);
    CHECK(matrix_gen_barabasi_albert(10, 0, 1) == NULL);
    MatrixGraph *tiny = matrix_gen_barabasi_albert(3, 5, 1);          /* n <= m: a complete graph */
    CHECK(tiny != NULL && matrix_graph_edge_count(tiny) == 3);
    matrix_graph_destroy(tiny);

    MatrixGraph *grid = matrix_gen_grid(4, 6);
    CHECK(matrix_graph_vertex_count(grid) == 24 && matrix_graph_edge_count(grid) == 4 * 5 + 3 * 6);
    CHECK(matrix_graph_degree(grid, 0) == 2 && matrix_graph_degree(grid, 7) == 4);
    matrix_graph_destroy(grid);

    MatrixGraph *p = matrix_gen_path(10, false), *cy = matrix_gen_cycle(10, false);
    MatrixGraph *dcy = matrix_gen_cycle(10, true), *st = matrix_gen_star(10), *kn = matrix_gen_complete(7, false);
    MatrixGraph *kd = matrix_gen_complete(7, true);
    CHECK(matrix_graph_edge_count(p) == 9 && matrix_graph_edge_count(cy) == 10);
    CHECK(matrix_graph_edge_count(dcy) == 10 && matrix_is_strongly_connected(dcy));
    CHECK(matrix_graph_edge_count(st) == 9 && matrix_graph_degree(st, 0) == 9);
    CHECK(matrix_graph_edge_count(kn) == 21 && matrix_graph_edge_count(kd) == 42);
    matrix_graph_destroy(p); matrix_graph_destroy(cy); matrix_graph_destroy(dcy);
    matrix_graph_destroy(st); matrix_graph_destroy(kn); matrix_graph_destroy(kd);

    MatrixGraph *c1 = matrix_gen_cycle(1, true), *c2 = matrix_gen_cycle(2, false), *z = matrix_gen_path(0, false);
    CHECK(matrix_graph_edge_count(c1) == 1 && matrix_graph_edge_count(c2) == 1 && matrix_graph_vertex_count(z) == 0);
    matrix_graph_destroy(c1); matrix_graph_destroy(c2); matrix_graph_destroy(z);

    for (int seed = 1; seed <= 10; ++seed) {
        MatrixGraph *dag = matrix_gen_random_dag(40, 0.2, (unsigned)seed, 1, 3);
        MatrixGraph *tree = matrix_gen_random_tree(40, (unsigned)seed);
        CHECK(matrix_is_dag(dag));
        CHECK(matrix_is_tree(tree));
        matrix_graph_destroy(dag);
        matrix_graph_destroy(tree);
    }
    CHECK(matrix_gen_random_dag(5, -0.1, 1, 1, 2) == NULL);
}

/* ---- text formats ---- */

static FILE *scratch(void)
{
    FILE *f = tmpfile();
    CHECK(f != NULL);
    return f;
}

static FILE *from_text(const char *text)
{
    FILE *f = scratch();
    fputs(text, f);
    rewind(f);
    return f;
}

static void test_dot(void)
{
    MatrixGraph *g = matrix_graph_create(3, true);
    CHECK(matrix_graph_add_edge(g, 0, 1, 2.5) && matrix_graph_add_edge(g, 1, 2, 1.0));
    CHECK(matrix_graph_remove_vertex(g, 2));
    const char *labels[] = { "Paris", "Ro\"me", NULL };
    FILE *f = scratch();
    CHECK_OK(matrix_write_dot(g, f, "demo", labels, true));
    rewind(f);
    char text[512] = { 0 };
    CHECK(fread(text, 1, sizeof(text) - 1, f) > 0);
    fclose(f);
    CHECK(strstr(text, "digraph \"demo\" {") != NULL);
    CHECK(strstr(text, "0 [label=\"Paris\"];") != NULL);
    CHECK(strstr(text, "1 [label=\"Ro\\\"me\"];") != NULL);          /* quote escaped */
    CHECK(strstr(text, "0 -> 1 [label=\"2.5\"];") != NULL);
    CHECK(strstr(text, " 2") == NULL);                              /* the dead vertex is omitted */
    matrix_graph_destroy(g);

    MatrixGraph *u = matrix_gen_path(3, false);
    f = scratch();
    CHECK_OK(matrix_write_dot(u, f, NULL, NULL, false));
    rewind(f);
    memset(text, 0, sizeof(text));
    CHECK(fread(text, 1, sizeof(text) - 1, f) > 0);
    fclose(f);
    CHECK(strstr(text, "graph \"G\" {") != NULL && strstr(text, "0 -- 1;") != NULL);
    CHECK(strstr(text, "1 -- 0") == NULL);                          /* undirected edge listed once */
    matrix_graph_destroy(u);
    CHECK(matrix_write_dot(NULL, stdout, "x", NULL, false) == GRAPH_ERR_INVALID_ARGUMENT);
}

static void test_edge_list(void)
{
    for (int seed = 1; seed <= 20; ++seed) {
        for (int directed = 0; directed < 2; ++directed) {
            MatrixGraph *g = matrix_gen_erdos_renyi(15, 0.2, (unsigned)seed, directed, 0.1, 99.9), *back = NULL;
            FILE *f = scratch();
            CHECK_OK(matrix_write_edge_list(g, f));
            rewind(f);
            CHECK_OK(matrix_read_edge_list(f, !directed, &back));      /* header overrides the argument */
            fclose(f);
            CHECK(same_graph(g, back));                                /* weights survive exactly */
            matrix_graph_destroy(g);
            matrix_graph_destroy(back);
        }
    }
    FILE *f = from_text("# a comment\n\n0 1 2.5\n  1   2\n2 0 0.25 \n# directed\n3 3\n");
    MatrixGraph *g = NULL;
    CHECK_OK(matrix_read_edge_list(f, false, &g));
    fclose(f);
    CHECK(matrix_graph_is_directed(g) && matrix_graph_vertex_count(g) == 4 && matrix_graph_edge_count(g) == 4);
    CHECK(wt(g, 0, 1) == 2.5 && wt(g, 1, 2) == 1.0 && wt(g, 2, 0) == 0.25 && matrix_graph_has_edge(g, 3, 3, NULL));
    matrix_graph_destroy(g);

    f = from_text("0 1 5\n1 0 9\n0 1 7\n");                          /* duplicates: first one wins */
    CHECK_OK(matrix_read_edge_list(f, false, &g));
    fclose(f);
    CHECK(matrix_graph_edge_count(g) == 1 && wt(g, 0, 1) == 5.0);
    matrix_graph_destroy(g);

    f = from_text("# vertices 8\n0 1\n");                            /* isolated vertices via the header */
    CHECK_OK(matrix_read_edge_list(f, false, &g));
    fclose(f);
    CHECK(matrix_graph_vertex_count(g) == 8);
    matrix_graph_destroy(g);

    f = from_text("");
    CHECK_OK(matrix_read_edge_list(f, false, &g));
    fclose(f);
    CHECK(matrix_graph_vertex_count(g) == 0);
    matrix_graph_destroy(g);

    const char *bad[] = { "0 x\n", "-1 2\n", "0 1 abc\n", "0\n", "0 1 nan\n", "0 1 2 3\n", "99999999999 1\n" };
    for (size_t i = 0; i < sizeof(bad) / sizeof(bad[0]); ++i) {
        f = from_text(bad[i]);
        g = NULL;
        CHECK(matrix_read_edge_list(f, false, &g) == GRAPH_ERR_IO);
        CHECK(g == NULL);
        fclose(f);
    }
    char longline[5000];                                              /* lines longer than any fixed buffer */
    memset(longline, ' ', sizeof(longline));
    strcpy(longline + 4000, "0 1 3\n");
    f = from_text(longline);
    CHECK_OK(matrix_read_edge_list(f, false, &g));
    fclose(f);
    CHECK(wt(g, 0, 1) == 3.0);
    matrix_graph_destroy(g);
}

static void test_named_edge_list(void)
{
    FILE *f = from_text("# undirected\nalice bob 2.5\nbob carol\ncarol alice 4\nalice dave\n");
    MatrixGraph *g = NULL;
    char **names = NULL;
    size_t count = 0;
    CHECK_OK(matrix_read_named_edge_list(f, true, &g, &names, &count));
    fclose(f);
    CHECK(count == 4 && !matrix_graph_is_directed(g) && matrix_graph_edge_count(g) == 4);
    CHECK(strcmp(names[0], "alice") == 0 && strcmp(names[1], "bob") == 0 &&
          strcmp(names[2], "carol") == 0 && strcmp(names[3], "dave") == 0);
    CHECK(wt(g, 0, 1) == 2.5 && wt(g, 1, 2) == 1.0 && wt(g, 2, 0) == 4.0);
    matrix_graph_destroy(g);
    matrix_names_free(names, count);

    f = from_text("a b\nb c d e\n");                                  /* too many fields */
    CHECK(matrix_read_named_edge_list(f, false, &g, &names, &count) == GRAPH_ERR_IO);
    fclose(f);
    f = from_text("a\n");
    CHECK(matrix_read_named_edge_list(f, false, &g, &names, &count) == GRAPH_ERR_IO);
    fclose(f);
    f = from_text("a b x\n");
    CHECK(matrix_read_named_edge_list(f, false, &g, &names, &count) == GRAPH_ERR_IO);
    fclose(f);

    /* Many names exercise the hash table's growth. */
    FILE *big = scratch();
    for (int i = 0; i < 2000; ++i) {
        fprintf(big, "node%d node%d\n", i, (i + 1) % 2000);
    }
    rewind(big);
    CHECK_OK(matrix_read_named_edge_list(big, false, &g, &names, &count));
    fclose(big);
    CHECK(count == 2000 && matrix_graph_edge_count(g) == 2000 && strcmp(names[1999], "node1999") == 0);
    matrix_graph_destroy(g);
    matrix_names_free(names, count);
    matrix_names_free(NULL, 0);
}

static void test_binary(void)
{
    for (int seed = 1; seed <= 30; ++seed) {
        for (int directed = 0; directed < 2; ++directed) {
            MatrixGraph *g = random_graph(20, 0.2, directed, 1000, (unsigned)seed, (size_t)(seed % 6));
            CHECK(matrix_graph_add_edge(g, 19, 19, 0.125) || !matrix_graph_is_alive(g, 19) || true);
            FILE *f = scratch();
            MatrixGraph *back = NULL;
            CHECK_OK(matrix_save_binary(g, f));
            rewind(f);
            CHECK_OK(matrix_load_binary(f, &back));
            fclose(f);
            CHECK(same_graph(g, back));                               /* dead IDs and all */
            matrix_graph_destroy(g);
            matrix_graph_destroy(back);
        }
    }
    /* A graph with irrational-looking weights must round-trip bit for bit. */
    MatrixGraph *g = matrix_graph_create(3, false), *back = NULL;
    CHECK(matrix_graph_add_edge(g, 0, 1, 0.1) && matrix_graph_add_edge(g, 1, 2, 1.0 / 3.0));
    FILE *f = scratch();
    CHECK_OK(matrix_save_binary(g, f));
    long size = ftell(f);
    rewind(f);
    CHECK_OK(matrix_load_binary(f, &back));
    fclose(f);
    CHECK(wt(back, 0, 1) == 0.1 && wt(back, 1, 2) == 1.0 / 3.0);
    matrix_graph_destroy(back);

    /* Truncated at every possible length, and with a corrupted magic number: always an error. */
    f = scratch();
    CHECK_OK(matrix_save_binary(g, f));
    unsigned char image[256];
    rewind(f);
    const size_t got = fread(image, 1, sizeof(image), f);
    fclose(f);
    CHECK((long)got == size);
    for (size_t cut = 0; cut < got; ++cut) {
        f = scratch();
        fwrite(image, 1, cut, f);
        rewind(f);
        back = NULL;
        CHECK(matrix_load_binary(f, &back) == GRAPH_ERR_IO && back == NULL);
        fclose(f);
    }
    image[0] ^= 0xFF;
    f = scratch();
    fwrite(image, 1, got, f);
    rewind(f);
    CHECK(matrix_load_binary(f, &back) == GRAPH_ERR_IO);
    fclose(f);
    matrix_graph_destroy(g);
}

static size_t count_substr(const char *hay, const char *needle)
{
    size_t c = 0;
    for (const char *p = hay; (p = strstr(p, needle)) != NULL; p += strlen(needle)) {
        ++c;
    }
    return c;
}

static void test_export_formats(void)
{
    MatrixGraph *g = matrix_graph_create(4, true);
    CHECK(matrix_graph_add_edge(g, 0, 1, 1.5) && matrix_graph_add_edge(g, 1, 2, 2.0) &&
          matrix_graph_add_edge(g, 2, 2, 3.0));
    CHECK(matrix_graph_remove_vertex(g, 3));
    char text[2048];
    FILE *f = scratch();

    CHECK_OK(matrix_write_graphml(g, f));
    rewind(f);
    memset(text, 0, sizeof(text));
    CHECK(fread(text, 1, sizeof(text) - 1, f) > 0);
    CHECK(strstr(text, "edgedefault=\"directed\"") && count_substr(text, "<node ") == 3 &&
          count_substr(text, "<edge ") == 3 && strstr(text, "source=\"n0\" target=\"n1\"") &&
          strstr(text, "<data key=\"w\">1.5</data>"));
    fclose(f);

    f = scratch();
    CHECK_OK(matrix_write_gml(g, f));
    rewind(f);
    memset(text, 0, sizeof(text));
    CHECK(fread(text, 1, sizeof(text) - 1, f) > 0);
    CHECK(strstr(text, "directed 1") && count_substr(text, "node [") == 3 && count_substr(text, "edge [") == 3 &&
          strstr(text, "source 0 target 1 weight 1.5"));
    fclose(f);

    f = scratch();
    CHECK_OK(matrix_write_json(g, f));
    rewind(f);
    memset(text, 0, sizeof(text));
    CHECK(fread(text, 1, sizeof(text) - 1, f) > 0);
    CHECK(strstr(text, "\"directed\": true") && strstr(text, "\"vertices\": [0, 1, 2]") &&
          count_substr(text, "\"source\"") == 3 && strstr(text, "\"weight\": 1.5"));
    fclose(f);
    matrix_graph_destroy(g);

    MatrixGraph *u = matrix_gen_path(3, false), *none = matrix_graph_create(2, false);
    f = scratch();
    CHECK_OK(matrix_write_json(u, f));
    rewind(f);
    memset(text, 0, sizeof(text));
    CHECK(fread(text, 1, sizeof(text) - 1, f) > 0);
    CHECK(strstr(text, "\"directed\": false") && count_substr(text, "\"source\"") == 2);
    fclose(f);
    f = scratch();
    CHECK_OK(matrix_write_json(none, f));
    rewind(f);
    memset(text, 0, sizeof(text));
    CHECK(fread(text, 1, sizeof(text) - 1, f) > 0 && strstr(text, "\"edges\": []"));
    fclose(f);
    matrix_graph_destroy(u);
    matrix_graph_destroy(none);
    CHECK(matrix_write_graphml(NULL, stdout) == GRAPH_ERR_INVALID_ARGUMENT);
    CHECK(matrix_write_gml(NULL, stdout) == GRAPH_ERR_INVALID_ARGUMENT);
    CHECK(matrix_write_json(NULL, stdout) == GRAPH_ERR_INVALID_ARGUMENT);
}

int main(void)
{
    test_copy_subgraph();
    test_set_operations();
    test_isomorphism();
    test_generators();
    test_dot();
    test_edge_list();
    test_named_edge_list();
    test_binary();
    test_export_formats();
    printf("All matrix ops and I/O tests passed.\n");
    return 0;
}
