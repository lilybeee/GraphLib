/*
 * test_matrix_combinatorial.c — Euler trails, Hamiltonian search, colouring,
 * cliques, vertex cover, independent sets, DAG longest paths.
 */

#include "test_util_matrix.h"

/* ---- Euler ---- */

/** Validate a trail: consecutive edges exist and every edge is used exactly once. */
static void check_euler(const MatrixGraph *g, const int *path, size_t len, bool closed)
{
    const size_t n = matrix_graph_id_bound(g);
    const bool directed = matrix_graph_is_directed(g);
    CHECK(len == matrix_graph_edge_count(g) + 1);
    unsigned char *used = calloc(n * n + 1, 1);
    for (size_t i = 0; i + 1 < len; ++i) {
        const size_t a = (size_t)path[i], b = (size_t)path[i + 1];
        CHECK(matrix_graph_has_edge(g, (int)a, (int)b, NULL));
        const size_t k = (directed || a <= b) ? a * n + b : b * n + a;
        CHECK(!used[k]);
        used[k] = 1;
    }
    if (closed) {
        CHECK(path[0] == path[len - 1]);
    }
    free(used);
}

/** Independent existence test: degree conditions plus connectivity of the edges. */
static bool euler_possible(const MatrixGraph *g, bool closed)
{
    const size_t n = matrix_graph_id_bound(g);
    const bool directed = matrix_graph_is_directed(g);
    if (matrix_graph_edge_count(g) == 0) {
        return true;
    }
    int odd = 0, plus = 0, minus = 0;
    MatrixGraph *edges_only = NULL;
    {   /* Keep only vertices that touch an edge; they must form one component. */
        bool drop[32] = { false };
        int map[32];
        for (size_t v = 0; v < n; ++v) {
            drop[v] = matrix_graph_degree(g, (int)v) == 0;
            for (size_t u = 0; u < n && drop[v]; ++u) {
                if (directed && matrix_graph_has_edge(g, (int)u, (int)v, NULL)) {
                    drop[v] = false;
                }
            }
        }
        edges_only = survivors(g, drop, map);
    }
    const bool connected = naive_components(edges_only) == 1;
    matrix_graph_destroy(edges_only);
    if (!connected) {
        return false;
    }
    for (size_t v = 0; v < n; ++v) {
        if (!matrix_graph_is_alive(g, (int)v)) {
            continue;
        }
        if (directed) {
            int out = 0, in = 0;
            for (size_t u = 0; u < n; ++u) {
                out += matrix_graph_has_edge(g, (int)v, (int)u, NULL);
                in  += matrix_graph_has_edge(g, (int)u, (int)v, NULL);
            }
            if (out - in == 1) ++plus;
            else if (in - out == 1) ++minus;
            else if (out != in) return false;
        } else {
            int d = 0;
            for (size_t u = 0; u < n; ++u) {
                d += matrix_graph_has_edge(g, (int)v, (int)u, NULL) * (u == v ? 2 : 1);
            }
            odd += d % 2;
        }
    }
    if (directed) {
        return (plus == 0 && minus == 0) || (!closed && plus == 1 && minus == 1);
    }
    return odd == 0 || (!closed && odd == 2);
}

static void test_euler(void)
{
    int path[64];
    size_t len;
    MatrixGraph *k5 = matrix_gen_complete(5, false);
    CHECK_OK(matrix_euler_trail(k5, true, path, &len));
    check_euler(k5, path, len, true);
    matrix_graph_destroy(k5);

    MatrixGraph *k4 = matrix_gen_complete(4, false);                  /* four odd vertices */
    CHECK(matrix_euler_trail(k4, false, path, &len) == GRAPH_ERR_NOT_FOUND);
    matrix_graph_destroy(k4);

    MatrixGraph *p = matrix_gen_path(5, false);                       /* trail yes, circuit no */
    CHECK_OK(matrix_euler_trail(p, false, path, &len));
    check_euler(p, path, len, false);
    CHECK(path[0] == 0 || path[0] == 4);
    CHECK(matrix_euler_trail(p, true, path, &len) == GRAPH_ERR_NOT_FOUND);
    matrix_graph_destroy(p);

    MatrixGraph *dc = matrix_gen_cycle(5, true);
    CHECK_OK(matrix_euler_trail(dc, true, path, &len));
    check_euler(dc, path, len, true);
    matrix_graph_destroy(dc);

    MatrixGraph *dp = matrix_gen_path(4, true);                       /* must start at 0, end at 3 */
    CHECK_OK(matrix_euler_trail(dp, false, path, &len));
    check_euler(dp, path, len, false);
    CHECK(path[0] == 0 && path[3] == 3);
    matrix_graph_destroy(dp);

    MatrixGraph *loops = matrix_graph_create(2, false);               /* a loop adds degree 2 */
    CHECK(matrix_graph_add_edge(loops, 0, 0, 1) && matrix_graph_add_edge(loops, 0, 1, 1));
    CHECK_OK(matrix_euler_trail(loops, false, path, &len));
    check_euler(loops, path, len, false);
    matrix_graph_destroy(loops);

    MatrixGraph *none = matrix_graph_create(3, false);
    CHECK_OK(matrix_euler_trail(none, true, path, &len));
    CHECK(len == 0);
    matrix_graph_destroy(none);

    MatrixGraph *two = matrix_graph_create(6, false);                 /* two triangles, disconnected */
    CHECK(matrix_graph_add_edge(two, 0, 1, 1) && matrix_graph_add_edge(two, 1, 2, 1) &&
          matrix_graph_add_edge(two, 2, 0, 1) && matrix_graph_add_edge(two, 3, 4, 1) &&
          matrix_graph_add_edge(two, 4, 5, 1) && matrix_graph_add_edge(two, 5, 3, 1));
    CHECK(matrix_euler_trail(two, false, path, &len) == GRAPH_ERR_NOT_FOUND);
    matrix_graph_destroy(two);

    for (int seed = 1; seed <= 300; ++seed) {
        for (int directed = 0; directed < 2; ++directed) {
            MatrixGraph *g = random_graph(7, 0.2 + 0.03 * (seed % 6), directed, 2, (unsigned)seed, (size_t)(seed % 2));
            for (int closed = 0; closed < 2; ++closed) {
                const MatrixStatus st = matrix_euler_trail(g, closed, path, &len);
                const bool want = euler_possible(g, closed);
                CHECK((st == GRAPH_OK) == want);
                if (st == GRAPH_OK && len > 0) {
                    check_euler(g, path, len, closed);
                }
            }
            matrix_graph_destroy(g);
        }
    }
}

/* ---- Hamiltonian ---- */

static bool brute_ham(const MatrixGraph *g, const int *ids, int m, bool cycle,
                      int last, unsigned used, int depth, int first)
{
    if (depth == m) {
        return !cycle || matrix_graph_has_edge(g, ids[last], ids[first], NULL);
    }
    for (int i = 0; i < m; ++i) {
        if (!((used >> i) & 1u) && matrix_graph_has_edge(g, ids[last], ids[i], NULL) &&
            brute_ham(g, ids, m, cycle, i, used | (1u << i), depth + 1, first)) {
            return true;
        }
    }
    return false;
}

static bool ham_exists(const MatrixGraph *g, bool cycle)
{
    int ids[16], m = 0;
    for (int v = 0; v < (int)matrix_graph_id_bound(g); ++v) {
        if (matrix_graph_is_alive(g, v)) {
            ids[m++] = v;
        }
    }
    if (m == 0) {
        return false;
    }
    if (m == 1) {
        return !cycle || matrix_graph_has_edge(g, ids[0], ids[0], NULL);
    }
    if (cycle && m < (matrix_graph_is_directed(g) ? 2 : 3)) {
        return false;
    }
    for (int s = 0; s < m; ++s) {
        if (brute_ham(g, ids, m, cycle, s, 1u << s, 1, s)) {
            return true;
        }
        if (cycle) {
            break;                      /* a cycle can start anywhere */
        }
    }
    return false;
}

static void check_ham(const MatrixGraph *g, const int *path, bool cycle)
{
    const size_t m = matrix_graph_vertex_count(g);
    bool seen[32] = { false };
    for (size_t i = 0; i < m; ++i) {
        CHECK(matrix_graph_is_alive(g, path[i]) && !seen[path[i]]);
        seen[path[i]] = true;
        if (i + 1 < m) {
            CHECK(matrix_graph_has_edge(g, path[i], path[i + 1], NULL));
        }
    }
    if (cycle && m > 1) {
        CHECK(matrix_graph_has_edge(g, path[m - 1], path[0], NULL));
    }
}

static void test_hamiltonian(void)
{
    int path[32];
    MatrixGraph *c = matrix_gen_cycle(8, false), *p = matrix_gen_path(8, false);
    CHECK_OK(matrix_hamiltonian(c, true, path));
    check_ham(c, path, true);
    CHECK_OK(matrix_hamiltonian(p, false, path));
    check_ham(p, path, false);
    CHECK(matrix_hamiltonian(p, true, path) == GRAPH_ERR_NOT_FOUND);
    matrix_graph_destroy(c);
    matrix_graph_destroy(p);

    MatrixGraph *star = matrix_gen_star(5);
    CHECK(matrix_hamiltonian(star, false, path) == GRAPH_ERR_NOT_FOUND);
    matrix_graph_destroy(star);

    MatrixGraph *big = matrix_gen_cycle(25, false);
    CHECK(matrix_hamiltonian(big, false, path) == GRAPH_ERR_TOO_LARGE);
    matrix_graph_destroy(big);

    MatrixGraph *edge = matrix_gen_path(2, false), *one = matrix_graph_create(1, false);
    CHECK(matrix_hamiltonian(edge, true, path) == GRAPH_ERR_NOT_FOUND);   /* 2 vertices: no simple cycle */
    CHECK_OK(matrix_hamiltonian(edge, false, path));
    CHECK_OK(matrix_hamiltonian(one, false, path));
    CHECK(path[0] == 0);
    CHECK(matrix_hamiltonian(one, true, path) == GRAPH_ERR_NOT_FOUND);
    matrix_graph_destroy(edge);
    matrix_graph_destroy(one);

    MatrixGraph *d2 = matrix_gen_cycle(2, true);                           /* 0 -> 1 -> 0 */
    CHECK_OK(matrix_hamiltonian(d2, true, path));
    matrix_graph_destroy(d2);

    for (int seed = 1; seed <= 200; ++seed) {
        for (int directed = 0; directed < 2; ++directed) {
            MatrixGraph *g = random_graph(7, 0.3 + 0.05 * (seed % 5), directed, 2, (unsigned)seed, (size_t)(seed % 2));
            for (int cycle = 0; cycle < 2; ++cycle) {
                const MatrixStatus st = matrix_hamiltonian(g, cycle, path);
                CHECK((st == GRAPH_OK) == ham_exists(g, cycle));
                if (st == GRAPH_OK) {
                    check_ham(g, path, cycle);
                }
            }
            matrix_graph_destroy(g);
        }
    }
}

/* ---- colouring ---- */

static size_t max_degree(const MatrixGraph *g)
{
    size_t m = 0;
    for (size_t v = 0; v < matrix_graph_id_bound(g); ++v) {
        size_t d = 0;
        for (size_t u = 0; u < matrix_graph_id_bound(g); ++u) {
            d += u != v && (matrix_graph_has_edge(g, (int)v, (int)u, NULL) ||
                            matrix_graph_has_edge(g, (int)u, (int)v, NULL));
        }
        if (d > m) {
            m = d;
        }
    }
    return m;
}

static void check_coloring(const MatrixGraph *g, const int *color, size_t ncolors)
{
    const size_t n = matrix_graph_id_bound(g);
    int top = -1;
    for (size_t v = 0; v < n; ++v) {
        if (!matrix_graph_is_alive(g, (int)v)) {
            CHECK(color[v] == -1);
            continue;
        }
        CHECK(color[v] >= 0);
        top = color[v] > top ? color[v] : top;
        for (size_t u = 0; u < n; ++u) {
            if (u != v && matrix_graph_is_alive(g, (int)u) &&
                (matrix_graph_has_edge(g, (int)v, (int)u, NULL) || matrix_graph_has_edge(g, (int)u, (int)v, NULL))) {
                CHECK(color[v] != color[u]);
            }
        }
    }
    CHECK(ncolors == (size_t)(top + 1));
}

static void test_coloring(void)
{
    int col[32];
    size_t k;
    MatrixGraph *k6 = matrix_gen_complete(6, false);
    CHECK_OK(matrix_color_greedy(k6, col, &k));
    CHECK(k == 6);
    CHECK_OK(matrix_color_welsh_powell(k6, col, &k));
    CHECK(k == 6);
    CHECK_OK(matrix_color_dsatur(k6, col, &k));
    CHECK(k == 6);
    matrix_graph_destroy(k6);

    MatrixGraph *odd = matrix_gen_cycle(7, false), *even = matrix_gen_cycle(8, false);
    CHECK_OK(matrix_color_dsatur(odd, col, &k));
    CHECK(k == 3);
    CHECK_OK(matrix_color_dsatur(even, col, &k));
    CHECK(k == 2);                                  /* DSatur is exact on bipartite graphs */
    matrix_graph_destroy(odd);
    matrix_graph_destroy(even);

    MatrixGraph *none = matrix_graph_create(4, false);
    CHECK_OK(matrix_color_greedy(none, col, &k));
    CHECK(k == 1);
    matrix_graph_destroy(none);
    MatrixGraph *zero = matrix_graph_create(0, false);
    CHECK_OK(matrix_color_dsatur(zero, col, &k));
    CHECK(k == 0);
    matrix_graph_destroy(zero);

    for (int seed = 1; seed <= 100; ++seed) {
        for (int directed = 0; directed < 2; ++directed) {
            MatrixGraph *g = random_graph(14, 0.15 + 0.03 * (seed % 8), directed, 2, (unsigned)seed, (size_t)(seed % 3));
            const size_t bound = max_degree(g) + 1;
            CHECK_OK(matrix_color_greedy(g, col, &k));
            check_coloring(g, col, k);
            CHECK(k <= bound);
            CHECK_OK(matrix_color_welsh_powell(g, col, &k));
            check_coloring(g, col, k);
            CHECK(k <= bound);
            CHECK_OK(matrix_color_dsatur(g, col, &k));
            check_coloring(g, col, k);
            CHECK(k <= bound);
            matrix_graph_destroy(g);
        }
    }
    for (int seed = 1; seed <= 20; ++seed) {            /* random bipartite graphs: DSatur uses <= 2 */
        t_seed((unsigned)seed + 4000u);
        MatrixGraph *g = matrix_graph_create(12, false);
        for (int u = 0; u < 6; ++u) {
            for (int v = 6; v < 12; ++v) {
                if (t_rand(100) < 40) {
                    CHECK(matrix_graph_add_edge(g, u, v, 1.0));
                }
            }
        }
        CHECK_OK(matrix_color_dsatur(g, col, &k));
        CHECK(k <= 2);
        matrix_graph_destroy(g);
    }
}

/* ---- cliques ---- */

static void test_cliques(void)
{
    for (int seed = 1; seed <= 100; ++seed) {
        MatrixGraph *g = random_graph(12, 0.25 + 0.04 * (seed % 8), false, 2, (unsigned)seed, (size_t)(seed % 3));
        const size_t n = matrix_graph_id_bound(g);
        /* Brute force: every vertex subset that is a clique and not extendable. */
        size_t expect = 0, best = 0;
        unsigned adj[12] = { 0 };
        for (size_t u = 0; u < n; ++u) {
            for (size_t v = 0; v < n; ++v) {
                if (u != v && matrix_graph_has_edge(g, (int)u, (int)v, NULL)) {
                    adj[u] |= 1u << v;
                }
            }
        }
        unsigned alive = 0;
        for (size_t v = 0; v < n; ++v) {
            alive |= matrix_graph_is_alive(g, (int)v) ? 1u << v : 0u;
        }
        unsigned *cliques = malloc((1u << n) * sizeof(unsigned));
        for (unsigned mask = 1; mask < (1u << n); ++mask) {
            if ((mask & ~alive) != 0) {
                continue;
            }
            bool is_clique = true;
            for (size_t v = 0; v < n && is_clique; ++v) {
                if ((mask >> v) & 1u) {
                    is_clique = (mask & ~(1u << v) & ~adj[v]) == 0;
                }
            }
            if (!is_clique) {
                continue;
            }
            bool maximal = true;
            for (size_t v = 0; v < n && maximal; ++v) {
                if (((alive >> v) & 1u) && !((mask >> v) & 1u) && (mask & ~adj[v]) == 0) {
                    maximal = false;
                }
            }
            if (maximal) {
                cliques[expect++] = mask;
                best = (size_t)__builtin_popcount(mask) > best ? (size_t)__builtin_popcount(mask) : best;
            }
        }
        MatrixGroups gr;
        CHECK_OK(matrix_maximal_cliques(g, &gr));
        CHECK(gr.count == expect);
        for (size_t i = 0; i < gr.count; ++i) {             /* every reported clique is in the brute-force list */
            unsigned mask = 0;
            for (size_t k = gr.offsets[i]; k < gr.offsets[i + 1]; ++k) {
                mask |= 1u << gr.members[k];
            }
            bool found = false;
            for (size_t j = 0; j < expect; ++j) {
                found = found || cliques[j] == mask;
            }
            CHECK(found);
        }
        int clique[12];
        size_t sz;
        CHECK_OK(matrix_max_clique(g, clique, &sz));
        CHECK(sz == best);
        matrix_groups_free(&gr);
        free(cliques);
        matrix_graph_destroy(g);
    }
    MatrixGraph *d = matrix_gen_complete(3, true);
    MatrixGroups gr;
    CHECK(matrix_maximal_cliques(d, &gr) == GRAPH_ERR_UNSUPPORTED);
    matrix_graph_destroy(d);
    MatrixGraph *e = matrix_graph_create(0, false);
    int c[1];
    size_t s;
    CHECK_OK(matrix_max_clique(e, c, &s));
    CHECK(s == 0);
    matrix_graph_destroy(e);
}

/* ---- vertex cover / independent set ---- */

static void test_cover_and_independent(void)
{
    for (int seed = 1; seed <= 100; ++seed) {
        for (int directed = 0; directed < 2; ++directed) {
            MatrixGraph *g = random_graph(12, 0.15 + 0.03 * (seed % 8), directed, 2, (unsigned)seed, (size_t)(seed % 3));
            if (seed % 7 == 0 && matrix_graph_is_alive(g, 3)) {
                (void)matrix_graph_add_edge(g, 3, 3, 1.0);        /* a self-loop */
            }
            const size_t n = matrix_graph_id_bound(g);
            bool cover[12], set[12];
            size_t cs, is;
            CHECK_OK(matrix_vertex_cover_approx(g, cover, &cs));
            CHECK_OK(matrix_independent_set_greedy(g, set, &is));

            size_t cc = 0, ic = 0;
            for (size_t u = 0; u < n; ++u) {
                cc += cover[u];
                ic += set[u];
                CHECK(!cover[u] || matrix_graph_is_alive(g, (int)u));
                CHECK(!set[u] || matrix_graph_is_alive(g, (int)u));
                for (size_t v = 0; v < n; ++v) {
                    if (matrix_graph_is_alive(g, (int)u) && matrix_graph_is_alive(g, (int)v) &&
                        matrix_graph_has_edge(g, (int)u, (int)v, NULL)) {
                        CHECK(cover[u] || cover[v]);               /* a real cover */
                        CHECK(!(set[u] && set[v]));                /* a real independent set (u==v: loop) */
                    }
                }
            }
            CHECK(cc == cs && ic == is);
            /* Independent set is maximal: every outside vertex has a chosen neighbour or a loop. */
            for (size_t v = 0; v < n; ++v) {
                if (!matrix_graph_is_alive(g, (int)v) || set[v]) {
                    continue;
                }
                bool blocked = matrix_graph_has_edge(g, (int)v, (int)v, NULL);
                for (size_t u = 0; u < n && !blocked; ++u) {
                    blocked = set[u] && (matrix_graph_has_edge(g, (int)u, (int)v, NULL) ||
                                         matrix_graph_has_edge(g, (int)v, (int)u, NULL));
                }
                CHECK(blocked);
            }
            /* 2-approximation: at most twice the minimum cover (brute force). */
            size_t opt = n + 1;
            for (unsigned mask = 0; mask < (1u << n); ++mask) {
                bool ok = true;
                for (size_t u = 0; u < n && ok; ++u) {
                    for (size_t v = 0; v < n && ok; ++v) {
                        if (matrix_graph_is_alive(g, (int)u) && matrix_graph_is_alive(g, (int)v) &&
                            matrix_graph_has_edge(g, (int)u, (int)v, NULL)) {
                            ok = ((mask >> u) & 1u) || ((mask >> v) & 1u);
                        }
                    }
                }
                if (ok && (size_t)__builtin_popcount(mask) < opt) {
                    opt = (size_t)__builtin_popcount(mask);
                }
            }
            CHECK(cs <= 2 * opt);
            matrix_graph_destroy(g);
        }
    }
}

/* ---- DAG longest path ---- */

static double brute_longest_from(const MatrixGraph *g, int u)
{
    double best = 0.0;
    size_t cur = 0;
    int v;
    double w;
    while (matrix_graph_next_neighbor(g, u, &cur, &v, &w)) {
        const double c = w + brute_longest_from(g, v);
        if (c > best) {
            best = c;
        }
    }
    return best;
}

static void test_dag_paths(void)
{
    for (int seed = 1; seed <= 100; ++seed) {
        MatrixGraph *g = matrix_gen_random_dag(9, 0.35, (unsigned)seed, -3.0, 7.0);
        if (seed % 4 == 0) {
            (void)matrix_graph_remove_vertex(g, seed % 9);
        }
        const size_t n = matrix_graph_id_bound(g);
        double dist[9], cost;
        int parent[9], end, path[9];
        size_t len;
        CHECK_OK(matrix_dag_longest_path(g, dist, parent, &end));
        double best = 0.0;
        for (size_t s = 0; s < n; ++s) {
            if (matrix_graph_is_alive(g, (int)s)) {
                best = fmax(best, brute_longest_from(g, (int)s));
            }
        }
        if (matrix_graph_vertex_count(g) > 0) {
            CHECK(end >= 0 && near(dist[end], best));
        }
        for (size_t v = 0; v < n; ++v) {
            if (!matrix_graph_is_alive(g, (int)v)) {
                CHECK(isinf(dist[v]) && dist[v] < 0 && parent[v] == -1);
            } else if (parent[v] >= 0) {
                CHECK(near(dist[v], dist[parent[v]] + wt(g, parent[v], (int)v)));
            }
        }
        CHECK_OK(matrix_critical_path(g, path, &len, &cost));
        CHECK(near(cost, best));
        double sum = 0;
        for (size_t i = 0; i + 1 < len; ++i) {
            sum += wt(g, path[i], path[i + 1]);
            CHECK(matrix_graph_has_edge(g, path[i], path[i + 1], NULL));
        }
        CHECK(near(sum, best));
        matrix_graph_destroy(g);
    }
    MatrixGraph *cyc = matrix_gen_cycle(4, true), *und = matrix_gen_path(4, false);
    double d[4], c;
    int p[4], e;
    size_t l;
    CHECK(matrix_dag_longest_path(cyc, d, p, &e) == GRAPH_ERR_UNSUPPORTED);
    CHECK(matrix_dag_longest_path(und, d, p, &e) == GRAPH_ERR_UNSUPPORTED);
    CHECK(matrix_critical_path(cyc, p, &l, &c) == GRAPH_ERR_UNSUPPORTED);
    matrix_graph_destroy(cyc);
    matrix_graph_destroy(und);
    MatrixGraph *empty = matrix_graph_create(0, true);
    CHECK_OK(matrix_critical_path(empty, p, &l, &c));
    CHECK(l == 0 && c == 0.0);
    matrix_graph_destroy(empty);
}

int main(void)
{
    test_euler();
    test_hamiltonian();
    test_coloring();
    test_cliques();
    test_cover_and_independent();
    test_dag_paths();
    printf("All matrix combinatorial tests passed.\n");
    return 0;
}
