/*
 * test_list_ops_io.c — copy, subgraph, union / intersection / complement,
 * isomorphism, generators, and every input / output format.
 */

#include "test_util_list.h"

/* ---- copy, subgraph, set operations ---- */

static void test_copy_subgraph(void)
{
    for (int seed = 1; seed <= 40; ++seed) {
        for (int directed = 0; directed < 2; ++directed) {
            ListGraph *g = random_graph(12, 0.25, directed, 9, (unsigned)seed, (size_t)(seed % 4));
            ListGraph *c = list_graph_copy(g);
            CHECK(c != NULL && c != g && same_graph(g, c));
            CHECK(list_graph_add_vertex(c) >= 0);                  /* the copy is independent */
            CHECK(!same_graph(g, c));
            list_graph_destroy(c);

            /* Induced subgraph on a shuffled subset of alive vertices. */
            int pick[12], k = 0;
            for (int v = 11; v >= 0; --v) {
                if (list_graph_is_alive(g, v) && (v + seed) % 3 != 0) {
                    pick[k++] = v;
                }
            }
            ListGraph *sub = NULL;
            CHECK_OK(list_induced_subgraph(g, pick, (size_t)k, &sub));
            CHECK(list_graph_vertex_count(sub) == (size_t)k && list_graph_id_bound(sub) == (size_t)k);
            size_t edges = 0;
            for (int i = 0; i < k; ++i) {
                for (int j = 0; j < k; ++j) {
                    const bool want = list_graph_has_edge(g, pick[i], pick[j], NULL);
                    CHECK(list_graph_has_edge(sub, i, j, NULL) == want);
                    CHECK(wt(sub, i, j) == wt(g, pick[i], pick[j]));
                    edges += want && (directed || i <= j);
                }
            }
            CHECK(list_graph_edge_count(sub) == edges);
            list_graph_destroy(sub);
            list_graph_destroy(g);
        }
    }
    ListGraph *g = list_gen_path(5, false), *out = NULL;
    const int dup[] = { 1, 1 }, dead[] = { 9 }, ok[] = { 2, 0 };
    CHECK(list_induced_subgraph(g, dup, 2, &out) == GRAPH_ERR_INVALID_ARGUMENT);
    CHECK(list_induced_subgraph(g, dead, 1, &out) == GRAPH_ERR_INVALID_ARGUMENT);
    CHECK_OK(list_induced_subgraph(g, ok, 2, &out));
    CHECK(list_graph_edge_count(out) == 0);
    list_graph_destroy(out);
    CHECK_OK(list_induced_subgraph(g, NULL, 0, &out));
    CHECK(list_graph_vertex_count(out) == 0);
    list_graph_destroy(out);
    CHECK(list_graph_copy(NULL) == NULL);
    list_graph_destroy(g);
}

static void test_set_operations(void)
{
    for (int seed = 1; seed <= 40; ++seed) {
        for (int directed = 0; directed < 2; ++directed) {
            ListGraph *a = random_graph(10, 0.25, directed, 9, (unsigned)seed, (size_t)(seed % 3));
            ListGraph *b = random_graph(12, 0.25, directed, 9, (unsigned)seed + 500u, (size_t)(seed % 2));
            ListGraph *u = NULL, *x = NULL, *comp = NULL;
            CHECK_OK(list_graph_union(a, b, &u));
            CHECK_OK(list_graph_intersection(a, b, &x));
            CHECK(list_graph_id_bound(u) == 12 && list_graph_id_bound(x) == 12);
            for (int s = 0; s < 12; ++s) {
                const bool in_a = s < 10 && list_graph_is_alive(a, s), in_b = list_graph_is_alive(b, s);
                CHECK(list_graph_is_alive(u, s) == (in_a || in_b));
                CHECK(list_graph_is_alive(x, s) == (in_a && in_b));
                for (int t = 0; t < 12; ++t) {
                    const bool ea = s < 10 && t < 10 && list_graph_has_edge(a, s, t, NULL);
                    const bool eb = list_graph_has_edge(b, s, t, NULL);
                    CHECK(list_graph_has_edge(u, s, t, NULL) == (ea || eb));
                    CHECK(list_graph_has_edge(x, s, t, NULL) == (ea && eb));
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
            CHECK_OK(list_graph_complement(a, &comp));
            for (int s = 0; s < 10; ++s) {
                CHECK(list_graph_is_alive(comp, s) == list_graph_is_alive(a, s));
                for (int t = 0; t < 10; ++t) {
                    const bool both = list_graph_is_alive(a, s) && list_graph_is_alive(a, t) && s != t;
                    CHECK(list_graph_has_edge(comp, s, t, NULL) == (both && !list_graph_has_edge(a, s, t, NULL)));
                }
            }
            ListGraph *back = NULL;                                /* complement twice: same structure */
            CHECK_OK(list_graph_complement(comp, &back));
            for (int s = 0; s < 10; ++s) {
                for (int t = 0; t < 10; ++t) {
                    CHECK(list_graph_has_edge(back, s, t, NULL) ==
                          (s != t && list_graph_has_edge(a, s, t, NULL)));
                }
            }
            list_graph_destroy(a); list_graph_destroy(b); list_graph_destroy(u);
            list_graph_destroy(x); list_graph_destroy(comp); list_graph_destroy(back);
        }
    }
    ListGraph *d = list_gen_path(3, true), *un = list_gen_path(3, false), *out = NULL;
    CHECK(list_graph_union(d, un, &out) == GRAPH_ERR_INVALID_ARGUMENT);
    CHECK(list_graph_intersection(d, un, &out) == GRAPH_ERR_INVALID_ARGUMENT);
    list_graph_destroy(d);
    list_graph_destroy(un);
}

/* ---- isomorphism ---- */

static ListGraph *relabel(const ListGraph *g, const int *perm)
{
    const size_t n = list_graph_id_bound(g);
    ListGraph *r = list_graph_create(n, list_graph_is_directed(g));
    for (size_t u = 0; u < n; ++u) {
        for (size_t v = 0; v < n; ++v) {
            if (list_graph_has_edge(g, (int)u, (int)v, NULL) &&
                !list_graph_has_edge(r, perm[u], perm[v], NULL)) {
                CHECK(list_graph_add_edge(r, perm[u], perm[v], 1.0));
            }
        }
    }
    return r;
}

/** Brute force: try every permutation of the vertices (n <= 8). */
static bool brute_iso_rec(const ListGraph *a, const ListGraph *b, int n, int *perm, unsigned used, int pos)
{
    if (pos == n) {
        for (int u = 0; u < n; ++u) {
            for (int v = 0; v < n; ++v) {
                if (list_graph_has_edge(a, u, v, NULL) != list_graph_has_edge(b, perm[u], perm[v], NULL)) {
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

static bool brute_iso(const ListGraph *a, const ListGraph *b)
{
    int perm[8];
    return brute_iso_rec(a, b, (int)list_graph_id_bound(a), perm, 0, 0);
}

static void test_isomorphism(void)
{
    /* Random pairs of small graphs against the permutation oracle. */
    int agree_true = 0, agree_false = 0;
    for (int seed = 1; seed <= 300; ++seed) {
        for (int directed = 0; directed < 2; ++directed) {
            ListGraph *a = random_graph(7, 0.35, directed, 1, (unsigned)seed, 0);
            ListGraph *b = random_graph(7, 0.35, directed, 1, (unsigned)seed + 9000u, 0);
            bool iso = false;
            CHECK_OK(list_is_isomorphic(a, b, &iso, NULL));
            const bool want = brute_iso(a, b);
            CHECK(iso == want);
            agree_true += want;
            agree_false += !want;
            list_graph_destroy(a);
            list_graph_destroy(b);
        }
    }
    CHECK(agree_true > 0 && agree_false > 0);            /* the sample really contains both answers */
    for (int seed = 1; seed <= 60; ++seed) {
        for (int directed = 0; directed < 2; ++directed) {
            ListGraph *a = random_graph(9, 0.3, directed, 1, (unsigned)seed, 0);
            int perm[9];
            for (int i = 0; i < 9; ++i) {
                perm[i] = i;
            }
            for (int i = 8; i > 0; --i) {
                const int j = (int)t_rand((unsigned)i + 1);
                const int tmp = perm[i]; perm[i] = perm[j]; perm[j] = tmp;
            }
            ListGraph *b = relabel(a, perm);
            bool iso = false;
            int map[9];
            CHECK_OK(list_is_isomorphic(a, b, &iso, map));
            CHECK(iso);
            for (int u = 0; u < 9; ++u) {                           /* the mapping is a real isomorphism */
                for (int v = 0; v < 9; ++v) {
                    CHECK(list_graph_has_edge(a, u, v, NULL) == list_graph_has_edge(b, map[u], map[v], NULL));
                }
            }
            list_graph_destroy(a);
            list_graph_destroy(b);
        }
    }
    /* Moving one edge keeps the counts; the answer must still match the oracle. */
    for (int seed = 1; seed <= 100; ++seed) {
        for (int directed = 0; directed < 2; ++directed) {
            ListGraph *a = random_graph(7, 0.4, directed, 1, (unsigned)seed + 300u, 0);
            ListGraph *b = list_graph_copy(a);
            int eu = -1, ev = -1;
            for (int u = 0; u < 7 && eu < 0; ++u) {
                for (int v = 0; v < 7 && eu < 0; ++v) {
                    if (u != v && list_graph_has_edge(b, u, v, NULL)) { eu = u; ev = v; }
                }
            }
            bool moved = false;
            if (eu >= 0) {
                CHECK(list_graph_remove_edge(b, eu, ev));
                for (int u = 0; u < 7 && !moved; ++u) {
                    for (int v = 0; v < 7 && !moved; ++v) {
                        if (u != v && !list_graph_has_edge(a, u, v, NULL) &&
                            !(!directed && list_graph_has_edge(a, v, u, NULL))) {
                            moved = list_graph_add_edge(b, u, v, 1.0);
                        }
                    }
                }
            }
            if (moved) {
                bool iso = true;
                CHECK_OK(list_is_isomorphic(a, b, &iso, NULL));
                CHECK(iso == brute_iso(a, b));
            }
            list_graph_destroy(a);
            list_graph_destroy(b);
        }
    }
    /* Same degree sequences, different structure. */
    ListGraph *c6 = list_gen_cycle(6, false), *two = list_graph_create(6, false);
    for (int i = 0; i < 3; ++i) {
        CHECK(list_graph_add_edge(two, i, (i + 1) % 3, 1) && list_graph_add_edge(two, 3 + i, 3 + (i + 1) % 3, 1));
    }
    bool iso;
    CHECK_OK(list_is_isomorphic(c6, two, &iso, NULL));
    CHECK(!iso);
    list_graph_destroy(c6);
    list_graph_destroy(two);

    ListGraph *k33 = list_graph_create(6, false), *prism = list_graph_create(6, false);
    for (int u = 0; u < 3; ++u) {
        for (int v = 3; v < 6; ++v) {
            CHECK(list_graph_add_edge(k33, u, v, 1));
        }
        CHECK(list_graph_add_edge(prism, u, (u + 1) % 3, 1) && list_graph_add_edge(prism, 3 + u, 3 + (u + 1) % 3, 1) &&
              list_graph_add_edge(prism, u, u + 3, 1));
    }
    CHECK_OK(list_is_isomorphic(k33, prism, &iso, NULL));          /* both 3-regular on 6 vertices */
    CHECK(!iso);
    list_graph_destroy(k33);
    list_graph_destroy(prism);

    ListGraph *d = list_gen_path(3, true), *u = list_gen_path(3, false);
    CHECK_OK(list_is_isomorphic(d, u, &iso, NULL));
    CHECK(!iso);                                                      /* directedness differs */
    list_graph_destroy(d);
    list_graph_destroy(u);

    ListGraph *l1 = list_graph_create(2, false), *l2 = list_graph_create(2, false);
    CHECK(list_graph_add_edge(l1, 0, 0, 1) && list_graph_add_edge(l1, 0, 1, 1));
    CHECK(list_graph_add_edge(l2, 1, 1, 1) && list_graph_add_edge(l2, 0, 1, 1));
    int map[2];
    CHECK_OK(list_is_isomorphic(l1, l2, &iso, map));               /* loops must map to loops */
    CHECK(iso && map[0] == 1 && map[1] == 0);
    list_graph_destroy(l1);
    list_graph_destroy(l2);
    CHECK(list_is_isomorphic(NULL, NULL, &iso, NULL) == GRAPH_ERR_INVALID_ARGUMENT);
}

/* ---- generators ---- */

static void test_generators(void)
{
    ListGraph *a = list_gen_erdos_renyi(60, 0.2, 7, false, 1.0, 5.0);
    ListGraph *b = list_gen_erdos_renyi(60, 0.2, 7, false, 1.0, 5.0);
    ListGraph *c = list_gen_erdos_renyi(60, 0.2, 8, false, 1.0, 5.0);
    CHECK(a && b && c && same_graph(a, b) && !same_graph(a, c));    /* seeded, reproducible */
    const double pairs = 60.0 * 59.0 / 2.0;
    CHECK(fabs((double)list_graph_edge_count(a) / pairs - 0.2) < 0.06);
    for (int u = 0; u < 60; ++u) {
        for (int v = 0; v < 60; ++v) {
            double w;
            if (list_graph_has_edge(a, u, v, &w)) {
                CHECK(w >= 1.0 && w <= 5.0 && u != v);
            }
        }
    }
    list_graph_destroy(a); list_graph_destroy(b); list_graph_destroy(c);

    ListGraph *empty = list_gen_erdos_renyi(10, 0.0, 1, true, 1, 1);
    ListGraph *full  = list_gen_erdos_renyi(10, 1.0, 1, true, 1, 1);
    CHECK(list_graph_edge_count(empty) == 0 && list_graph_edge_count(full) == 90);
    list_graph_destroy(empty); list_graph_destroy(full);
    CHECK(list_gen_erdos_renyi(5, 1.5, 1, false, 1, 1) == NULL);
    CHECK(list_gen_erdos_renyi(5, 0.5, 1, false, 3, 1) == NULL);
    CHECK(list_gen_erdos_renyi(5, 0.5, 1, false, 1, INFINITY) == NULL);

    ListGraph *ba = list_gen_barabasi_albert(100, 3, 42);
    CHECK(ba != NULL && list_graph_vertex_count(ba) == 100 && list_is_connected(ba));
    CHECK(list_graph_edge_count(ba) == 6 + 96 * 3);                 /* K4 seed + 3 per newcomer */
    size_t minimum = 1000;
    for (int v = 0; v < 100; ++v) {
        if (list_graph_degree(ba, v) < minimum) {
            minimum = list_graph_degree(ba, v);
        }
    }
    CHECK(minimum >= 3);
    ListGraph *ba2 = list_gen_barabasi_albert(100, 3, 42);
    CHECK(same_graph(ba, ba2));
    list_graph_destroy(ba); list_graph_destroy(ba2);
    CHECK(list_gen_barabasi_albert(10, 0, 1) == NULL);
    ListGraph *tiny = list_gen_barabasi_albert(3, 5, 1);          /* n <= m: a complete graph */
    CHECK(tiny != NULL && list_graph_edge_count(tiny) == 3);
    list_graph_destroy(tiny);

    ListGraph *grid = list_gen_grid(4, 6);
    CHECK(list_graph_vertex_count(grid) == 24 && list_graph_edge_count(grid) == 4 * 5 + 3 * 6);
    CHECK(list_graph_degree(grid, 0) == 2 && list_graph_degree(grid, 7) == 4);
    list_graph_destroy(grid);

    ListGraph *p = list_gen_path(10, false), *cy = list_gen_cycle(10, false);
    ListGraph *dcy = list_gen_cycle(10, true), *st = list_gen_star(10), *kn = list_gen_complete(7, false);
    ListGraph *kd = list_gen_complete(7, true);
    CHECK(list_graph_edge_count(p) == 9 && list_graph_edge_count(cy) == 10);
    CHECK(list_graph_edge_count(dcy) == 10 && list_is_strongly_connected(dcy));
    CHECK(list_graph_edge_count(st) == 9 && list_graph_degree(st, 0) == 9);
    CHECK(list_graph_edge_count(kn) == 21 && list_graph_edge_count(kd) == 42);
    list_graph_destroy(p); list_graph_destroy(cy); list_graph_destroy(dcy);
    list_graph_destroy(st); list_graph_destroy(kn); list_graph_destroy(kd);

    ListGraph *c1 = list_gen_cycle(1, true), *c2 = list_gen_cycle(2, false), *z = list_gen_path(0, false);
    CHECK(list_graph_edge_count(c1) == 1 && list_graph_edge_count(c2) == 1 && list_graph_vertex_count(z) == 0);
    list_graph_destroy(c1); list_graph_destroy(c2); list_graph_destroy(z);

    for (int seed = 1; seed <= 10; ++seed) {
        ListGraph *dag = list_gen_random_dag(40, 0.2, (unsigned)seed, 1, 3);
        ListGraph *tree = list_gen_random_tree(40, (unsigned)seed);
        CHECK(list_is_dag(dag));
        CHECK(list_is_tree(tree));
        list_graph_destroy(dag);
        list_graph_destroy(tree);
    }
    CHECK(list_gen_random_dag(5, -0.1, 1, 1, 2) == NULL);
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
    ListGraph *g = list_graph_create(3, true);
    CHECK(list_graph_add_edge(g, 0, 1, 2.5) && list_graph_add_edge(g, 1, 2, 1.0));
    CHECK(list_graph_remove_vertex(g, 2));
    const char *labels[] = { "Paris", "Ro\"me", NULL };
    FILE *f = scratch();
    CHECK_OK(list_write_dot(g, f, "demo", labels, true));
    rewind(f);
    char text[512] = { 0 };
    CHECK(fread(text, 1, sizeof(text) - 1, f) > 0);
    fclose(f);
    CHECK(strstr(text, "digraph \"demo\" {") != NULL);
    CHECK(strstr(text, "0 [label=\"Paris\"];") != NULL);
    CHECK(strstr(text, "1 [label=\"Ro\\\"me\"];") != NULL);          /* quote escaped */
    CHECK(strstr(text, "0 -> 1 [label=\"2.5\"];") != NULL);
    CHECK(strstr(text, " 2") == NULL);                              /* the dead vertex is omitted */
    list_graph_destroy(g);

    ListGraph *u = list_gen_path(3, false);
    f = scratch();
    CHECK_OK(list_write_dot(u, f, NULL, NULL, false));
    rewind(f);
    memset(text, 0, sizeof(text));
    CHECK(fread(text, 1, sizeof(text) - 1, f) > 0);
    fclose(f);
    CHECK(strstr(text, "graph \"G\" {") != NULL && strstr(text, "0 -- 1;") != NULL);
    CHECK(strstr(text, "1 -- 0") == NULL);                          /* undirected edge listed once */
    list_graph_destroy(u);
    CHECK(list_write_dot(NULL, stdout, "x", NULL, false) == GRAPH_ERR_INVALID_ARGUMENT);
}

static void test_edge_list(void)
{
    for (int seed = 1; seed <= 20; ++seed) {
        for (int directed = 0; directed < 2; ++directed) {
            ListGraph *g = list_gen_erdos_renyi(15, 0.2, (unsigned)seed, directed, 0.1, 99.9), *back = NULL;
            FILE *f = scratch();
            CHECK_OK(list_write_edge_list(g, f));
            rewind(f);
            CHECK_OK(list_read_edge_list(f, !directed, &back));      /* header overrides the argument */
            fclose(f);
            CHECK(same_graph(g, back));                                /* weights survive exactly */
            list_graph_destroy(g);
            list_graph_destroy(back);
        }
    }
    FILE *f = from_text("# a comment\n\n0 1 2.5\n  1   2\n2 0 0.25 \n# directed\n3 3\n");
    ListGraph *g = NULL;
    CHECK_OK(list_read_edge_list(f, false, &g));
    fclose(f);
    CHECK(list_graph_is_directed(g) && list_graph_vertex_count(g) == 4 && list_graph_edge_count(g) == 4);
    CHECK(wt(g, 0, 1) == 2.5 && wt(g, 1, 2) == 1.0 && wt(g, 2, 0) == 0.25 && list_graph_has_edge(g, 3, 3, NULL));
    list_graph_destroy(g);

    f = from_text("0 1 5\n1 0 9\n0 1 7\n");                          /* duplicates: first one wins */
    CHECK_OK(list_read_edge_list(f, false, &g));
    fclose(f);
    CHECK(list_graph_edge_count(g) == 1 && wt(g, 0, 1) == 5.0);
    list_graph_destroy(g);

    f = from_text("# vertices 8\n0 1\n");                            /* isolated vertices via the header */
    CHECK_OK(list_read_edge_list(f, false, &g));
    fclose(f);
    CHECK(list_graph_vertex_count(g) == 8);
    list_graph_destroy(g);

    f = from_text("");
    CHECK_OK(list_read_edge_list(f, false, &g));
    fclose(f);
    CHECK(list_graph_vertex_count(g) == 0);
    list_graph_destroy(g);

    const char *bad[] = { "0 x\n", "-1 2\n", "0 1 abc\n", "0\n", "0 1 nan\n", "0 1 2 3\n", "99999999999 1\n" };
    for (size_t i = 0; i < sizeof(bad) / sizeof(bad[0]); ++i) {
        f = from_text(bad[i]);
        g = NULL;
        CHECK(list_read_edge_list(f, false, &g) == GRAPH_ERR_IO);
        CHECK(g == NULL);
        fclose(f);
    }
    char longline[5000];                                              /* lines longer than any fixed buffer */
    memset(longline, ' ', sizeof(longline));
    strcpy(longline + 4000, "0 1 3\n");
    f = from_text(longline);
    CHECK_OK(list_read_edge_list(f, false, &g));
    fclose(f);
    CHECK(wt(g, 0, 1) == 3.0);
    list_graph_destroy(g);
}

static void test_named_edge_list(void)
{
    FILE *f = from_text("# undirected\nalice bob 2.5\nbob carol\ncarol alice 4\nalice dave\n");
    ListGraph *g = NULL;
    char **names = NULL;
    size_t count = 0;
    CHECK_OK(list_read_named_edge_list(f, true, &g, &names, &count));
    fclose(f);
    CHECK(count == 4 && !list_graph_is_directed(g) && list_graph_edge_count(g) == 4);
    CHECK(strcmp(names[0], "alice") == 0 && strcmp(names[1], "bob") == 0 &&
          strcmp(names[2], "carol") == 0 && strcmp(names[3], "dave") == 0);
    CHECK(wt(g, 0, 1) == 2.5 && wt(g, 1, 2) == 1.0 && wt(g, 2, 0) == 4.0);
    list_graph_destroy(g);
    list_names_free(names, count);

    f = from_text("a b\nb c d e\n");                                  /* too many fields */
    CHECK(list_read_named_edge_list(f, false, &g, &names, &count) == GRAPH_ERR_IO);
    fclose(f);
    f = from_text("a\n");
    CHECK(list_read_named_edge_list(f, false, &g, &names, &count) == GRAPH_ERR_IO);
    fclose(f);
    f = from_text("a b x\n");
    CHECK(list_read_named_edge_list(f, false, &g, &names, &count) == GRAPH_ERR_IO);
    fclose(f);

    /* Many names exercise the hash table's growth. */
    FILE *big = scratch();
    for (int i = 0; i < 2000; ++i) {
        fprintf(big, "node%d node%d\n", i, (i + 1) % 2000);
    }
    rewind(big);
    CHECK_OK(list_read_named_edge_list(big, false, &g, &names, &count));
    fclose(big);
    CHECK(count == 2000 && list_graph_edge_count(g) == 2000 && strcmp(names[1999], "node1999") == 0);
    list_graph_destroy(g);
    list_names_free(names, count);
    list_names_free(NULL, 0);
}

static void test_binary(void)
{
    for (int seed = 1; seed <= 30; ++seed) {
        for (int directed = 0; directed < 2; ++directed) {
            ListGraph *g = random_graph(20, 0.2, directed, 1000, (unsigned)seed, (size_t)(seed % 6));
            CHECK(list_graph_add_edge(g, 19, 19, 0.125) || !list_graph_is_alive(g, 19) || true);
            FILE *f = scratch();
            ListGraph *back = NULL;
            CHECK_OK(list_save_binary(g, f));
            rewind(f);
            CHECK_OK(list_load_binary(f, &back));
            fclose(f);
            CHECK(same_graph(g, back));                               /* dead IDs and all */
            list_graph_destroy(g);
            list_graph_destroy(back);
        }
    }
    /* A graph with irrational-looking weights must round-trip bit for bit. */
    ListGraph *g = list_graph_create(3, false), *back = NULL;
    CHECK(list_graph_add_edge(g, 0, 1, 0.1) && list_graph_add_edge(g, 1, 2, 1.0 / 3.0));
    FILE *f = scratch();
    CHECK_OK(list_save_binary(g, f));
    long size = ftell(f);
    rewind(f);
    CHECK_OK(list_load_binary(f, &back));
    fclose(f);
    CHECK(wt(back, 0, 1) == 0.1 && wt(back, 1, 2) == 1.0 / 3.0);
    list_graph_destroy(back);

    /* Truncated at every possible length, and with a corrupted magic number: always an error. */
    f = scratch();
    CHECK_OK(list_save_binary(g, f));
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
        CHECK(list_load_binary(f, &back) == GRAPH_ERR_IO && back == NULL);
        fclose(f);
    }
    image[0] ^= 0xFF;
    f = scratch();
    fwrite(image, 1, got, f);
    rewind(f);
    CHECK(list_load_binary(f, &back) == GRAPH_ERR_IO);
    fclose(f);
    list_graph_destroy(g);
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
    ListGraph *g = list_graph_create(4, true);
    CHECK(list_graph_add_edge(g, 0, 1, 1.5) && list_graph_add_edge(g, 1, 2, 2.0) &&
          list_graph_add_edge(g, 2, 2, 3.0));
    CHECK(list_graph_remove_vertex(g, 3));
    char text[2048];
    FILE *f = scratch();

    CHECK_OK(list_write_graphml(g, f));
    rewind(f);
    memset(text, 0, sizeof(text));
    CHECK(fread(text, 1, sizeof(text) - 1, f) > 0);
    CHECK(strstr(text, "edgedefault=\"directed\"") && count_substr(text, "<node ") == 3 &&
          count_substr(text, "<edge ") == 3 && strstr(text, "source=\"n0\" target=\"n1\"") &&
          strstr(text, "<data key=\"w\">1.5</data>"));
    fclose(f);

    f = scratch();
    CHECK_OK(list_write_gml(g, f));
    rewind(f);
    memset(text, 0, sizeof(text));
    CHECK(fread(text, 1, sizeof(text) - 1, f) > 0);
    CHECK(strstr(text, "directed 1") && count_substr(text, "node [") == 3 && count_substr(text, "edge [") == 3 &&
          strstr(text, "source 0 target 1 weight 1.5"));
    fclose(f);

    f = scratch();
    CHECK_OK(list_write_json(g, f));
    rewind(f);
    memset(text, 0, sizeof(text));
    CHECK(fread(text, 1, sizeof(text) - 1, f) > 0);
    CHECK(strstr(text, "\"directed\": true") && strstr(text, "\"vertices\": [0, 1, 2]") &&
          count_substr(text, "\"source\"") == 3 && strstr(text, "\"weight\": 1.5"));
    fclose(f);
    list_graph_destroy(g);

    ListGraph *u = list_gen_path(3, false), *none = list_graph_create(2, false);
    f = scratch();
    CHECK_OK(list_write_json(u, f));
    rewind(f);
    memset(text, 0, sizeof(text));
    CHECK(fread(text, 1, sizeof(text) - 1, f) > 0);
    CHECK(strstr(text, "\"directed\": false") && count_substr(text, "\"source\"") == 2);
    fclose(f);
    f = scratch();
    CHECK_OK(list_write_json(none, f));
    rewind(f);
    memset(text, 0, sizeof(text));
    CHECK(fread(text, 1, sizeof(text) - 1, f) > 0 && strstr(text, "\"edges\": []"));
    fclose(f);
    list_graph_destroy(u);
    list_graph_destroy(none);
    CHECK(list_write_graphml(NULL, stdout) == GRAPH_ERR_INVALID_ARGUMENT);
    CHECK(list_write_gml(NULL, stdout) == GRAPH_ERR_INVALID_ARGUMENT);
    CHECK(list_write_json(NULL, stdout) == GRAPH_ERR_INVALID_ARGUMENT);
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
    printf("All list ops and I/O tests passed.\n");
    return 0;
}
