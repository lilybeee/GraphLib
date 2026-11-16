/*
 * test_matrix_network.c — PageRank, HITS, eigenvector centrality, k-cores,
 * label propagation, Louvain, modularity, assortativity.  The iterative
 * methods are compared with plain dense-matrix implementations.
 */

#include "test_util_matrix.h"

#define MAXN 24

static void dense(const MatrixGraph *g, double *A)       /* A[u*n+v] = weight or 0 */
{
    const size_t n = matrix_graph_id_bound(g);
    for (size_t u = 0; u < n; ++u) {
        for (size_t v = 0; v < n; ++v) {
            A[u * n + v] = wt(g, (int)u, (int)v);
        }
    }
}

/* ---- PageRank ---- */

static void dense_pagerank(const MatrixGraph *g, double d, double *r)
{
    const size_t n = matrix_graph_id_bound(g);
    double A[MAXN * MAXN], next[MAXN];
    dense(g, A);
    size_t alive = 0;
    for (size_t v = 0; v < n; ++v) {
        alive += matrix_graph_is_alive(g, (int)v);
    }
    for (size_t v = 0; v < n; ++v) {
        r[v] = matrix_graph_is_alive(g, (int)v) ? 1.0 / (double)alive : 0.0;
    }
    for (int it = 0; it < 2000; ++it) {
        double dangling = 0.0;
        for (size_t u = 0; u < n; ++u) {
            double out = 0;
            for (size_t v = 0; v < n; ++v) {
                out += A[u * n + v];
            }
            if (matrix_graph_is_alive(g, (int)u) && out == 0.0) {
                dangling += r[u];
            }
        }
        for (size_t v = 0; v < n; ++v) {
            next[v] = 0.0;
            if (!matrix_graph_is_alive(g, (int)v)) {
                continue;
            }
            next[v] = (1.0 - d) / (double)alive + d * dangling / (double)alive;
            for (size_t u = 0; u < n; ++u) {
                double out = 0;
                for (size_t w = 0; w < n; ++w) {
                    out += A[u * n + w];
                }
                if (out > 0.0) {
                    next[v] += d * r[u] * A[u * n + v] / out;
                }
            }
        }
        memcpy(r, next, n * sizeof(double));
    }
}

static void test_pagerank(void)
{
    double pr[MAXN], want[MAXN];
    size_t iters = 0;

    MatrixGraph *c3 = matrix_gen_cycle(3, true);                    /* symmetric: all equal */
    CHECK_OK(matrix_pagerank(c3, 0.85, 1e-12, 1000, pr, &iters));
    CHECK(near(pr[0], 1.0 / 3) && near(pr[1], 1.0 / 3) && iters > 0);
    matrix_graph_destroy(c3);

    MatrixGraph *star = matrix_gen_star(6);
    CHECK_OK(matrix_pagerank(star, 0.85, 1e-12, 1000, pr, NULL));
    CHECK(pr[0] > pr[1] && near(pr[1], pr[2]));
    double sum = 0;
    for (int i = 0; i < 6; ++i) {
        sum += pr[i];
    }
    CHECK(near(sum, 1.0));
    matrix_graph_destroy(star);

    for (int seed = 1; seed <= 50; ++seed) {
        for (int directed = 0; directed < 2; ++directed) {
            MatrixGraph *g = random_graph(12, 0.2, directed, 5, (unsigned)seed, (size_t)(seed % 3));
            CHECK_OK(matrix_pagerank(g, 0.85, 1e-13, 5000, pr, NULL));
            dense_pagerank(g, 0.85, want);
            sum = 0;
            for (size_t v = 0; v < matrix_graph_id_bound(g); ++v) {
                CHECK(fabs(pr[v] - want[v]) < 1e-9);
                sum += pr[v];
            }
            CHECK(matrix_graph_vertex_count(g) == 0 || near(sum, 1.0));
            matrix_graph_destroy(g);
        }
    }
    MatrixGraph *g = matrix_gen_path(3, true);
    CHECK(matrix_pagerank(g, 1.5, 1e-9, 10, pr, NULL) == GRAPH_ERR_INVALID_ARGUMENT);
    CHECK(matrix_pagerank(g, 0.85, 0.0, 10, pr, NULL) == GRAPH_ERR_INVALID_ARGUMENT);
    CHECK(matrix_pagerank(g, 0.85, 1e-30, 2, pr, NULL) == GRAPH_ERR_NO_CONVERGENCE);
    CHECK(matrix_graph_add_edge(g, 0, 2, -1.0));
    CHECK(matrix_pagerank(g, 0.85, 1e-9, 100, pr, NULL) == GRAPH_ERR_INVALID_ARGUMENT);
    matrix_graph_destroy(g);
}

/* ---- HITS / eigenvector ---- */

static void l2(double *x, size_t n)
{
    double s = 0;
    for (size_t i = 0; i < n; ++i) {
        s += x[i] * x[i];
    }
    s = sqrt(s);
    for (size_t i = 0; s > 0 && i < n; ++i) {
        x[i] /= s;
    }
}

static void test_hits_eigenvector(void)
{
    for (int seed = 1; seed <= 40; ++seed) {
        for (int directed = 0; directed < 2; ++directed) {
            MatrixGraph *g = random_graph(10, 0.3, directed, 1, (unsigned)seed, (size_t)(seed % 2));
            const size_t n = matrix_graph_id_bound(g);
            double A[MAXN * MAXN], hub[MAXN], auth[MAXN], h[MAXN], a[MAXN], tmp[MAXN];
            dense(g, A);
            CHECK_OK(matrix_hits(g, 1e-14, 20000, hub, auth));

            /* Reference power iteration on A A^T and A^T A. */
            for (size_t i = 0; i < n; ++i) {
                h[i] = a[i] = matrix_graph_is_alive(g, (int)i) ? 1.0 : 0.0;
            }
            l2(h, n);
            for (int it = 0; it < 20000; ++it) {
                for (size_t v = 0; v < n; ++v) {
                    tmp[v] = 0;
                    for (size_t u = 0; u < n; ++u) {
                        tmp[v] += A[u * n + v] * h[u];
                    }
                }
                l2(tmp, n);
                memcpy(a, tmp, n * sizeof(double));
                for (size_t u = 0; u < n; ++u) {
                    tmp[u] = 0;
                    for (size_t v = 0; v < n; ++v) {
                        tmp[u] += A[u * n + v] * a[v];
                    }
                }
                l2(tmp, n);
                memcpy(h, tmp, n * sizeof(double));
            }
            for (size_t i = 0; i < n; ++i) {
                CHECK(fabs(hub[i] - h[i]) < 1e-6 && fabs(auth[i] - a[i]) < 1e-6);
            }

            /* Eigenvector centrality: reference power iteration on (A^T + I).  A symmetric
             * matrix always converges; a directed one may converge only polynomially, so it
             * is allowed to report NO_CONVERGENCE but must still return a sane vector. */
            double ev[MAXN], x[MAXN];
            const MatrixStatus st = matrix_eigenvector_centrality(g, 1e-13, 50000, ev);
            CHECK(directed ? (st == GRAPH_OK || st == GRAPH_ERR_NO_CONVERGENCE) : st == GRAPH_OK);
            double norm = 0;
            for (size_t i = 0; i < n; ++i) {
                CHECK(ev[i] >= 0.0 && isfinite(ev[i]));
                norm += ev[i] * ev[i];
            }
            CHECK(matrix_graph_vertex_count(g) == 0 || near(norm, 1.0));
            if (st == GRAPH_OK) {
                for (size_t i = 0; i < n; ++i) {
                    x[i] = matrix_graph_is_alive(g, (int)i) ? 1.0 : 0.0;
                }
                l2(x, n);
                for (int it = 0; it < 80000; ++it) {
                    for (size_t v = 0; v < n; ++v) {
                        tmp[v] = x[v];
                        for (size_t u = 0; u < n; ++u) {
                            tmp[v] += (A[u * n + v] != 0.0 ? 1.0 : 0.0) * x[u];
                        }
                    }
                    l2(tmp, n);
                    memcpy(x, tmp, n * sizeof(double));
                }
                for (size_t i = 0; i < n; ++i) {
                    CHECK(fabs(ev[i] - x[i]) < 1e-5);
                }
            }
            matrix_graph_destroy(g);
        }
    }
    /* Path of three: eigenvector (1/2, 1/sqrt2, 1/2). */
    MatrixGraph *p = matrix_gen_path(3, false);
    double ev[3];
    CHECK_OK(matrix_eigenvector_centrality(p, 1e-13, 10000, ev));
    CHECK(fabs(ev[0] - 0.5) < 1e-6 && fabs(ev[1] - sqrt(0.5)) < 1e-6 && fabs(ev[2] - 0.5) < 1e-6);
    matrix_graph_destroy(p);

    /* A source node that points at everything is the best hub and not an authority. */
    MatrixGraph *d = matrix_graph_create(4, true);
    for (int v = 1; v < 4; ++v) {
        CHECK(matrix_graph_add_edge(d, 0, v, 1.0));
    }
    double hub[4], auth[4];
    CHECK_OK(matrix_hits(d, 1e-12, 1000, hub, auth));
    CHECK(hub[0] > 0.99 && auth[0] < 1e-9 && auth[1] > 0.5);
    CHECK(matrix_hits(d, 1e-30, 1, hub, auth) == GRAPH_ERR_NO_CONVERGENCE);
    matrix_graph_destroy(d);
}

/* ---- k-core ---- */

static void test_cores(void)
{
    for (int seed = 1; seed <= 100; ++seed) {
        for (int directed = 0; directed < 2; ++directed) {
            MatrixGraph *g = random_graph(14, 0.1 + 0.03 * (seed % 8), directed, 1, (unsigned)seed, (size_t)(seed % 3));
            const size_t n = matrix_graph_id_bound(g);
            int core[14], want[14];
            CHECK_OK(matrix_core_numbers(g, core));

            /* Naive peeling: repeatedly delete a minimum-degree vertex. */
            bool gone[14];
            int deg[14], k = 0;
            for (size_t v = 0; v < n; ++v) {
                gone[v] = !matrix_graph_is_alive(g, (int)v);
                want[v] = -1;
            }
            for (;;) {
                int best = -1;
                for (size_t v = 0; v < n; ++v) {
                    deg[v] = 0;
                    for (size_t u = 0; u < n; ++u) {
                        deg[v] += !gone[v] && !gone[u] && u != v &&
                                  (matrix_graph_has_edge(g, (int)v, (int)u, NULL) ||
                                   matrix_graph_has_edge(g, (int)u, (int)v, NULL));
                    }
                    if (!gone[v] && (best < 0 || deg[v] < deg[best])) {
                        best = (int)v;
                    }
                }
                if (best < 0) {
                    break;
                }
                k = deg[best] > k ? deg[best] : k;
                want[best] = k;
                gone[best] = true;
            }
            for (size_t v = 0; v < n; ++v) {
                CHECK(core[v] == want[v]);
            }
            matrix_graph_destroy(g);
        }
    }
    MatrixGraph *k5 = matrix_gen_complete(5, false);
    int core[5];
    CHECK_OK(matrix_core_numbers(k5, core));
    CHECK(core[0] == 4 && core[4] == 4);
    matrix_graph_destroy(k5);
}

/* ---- communities ---- */

/** `cliques` cliques of `size` vertices joined in a ring by single edges. */
static MatrixGraph *ring_of_cliques(int cliques, int size)
{
    MatrixGraph *g = matrix_graph_create((size_t)(cliques * size), false);
    for (int c = 0; c < cliques; ++c) {
        for (int a = 0; a < size; ++a) {
            for (int b = a + 1; b < size; ++b) {
                CHECK(matrix_graph_add_edge(g, c * size + a, c * size + b, 1.0));
            }
        }
        CHECK(matrix_graph_add_edge(g, c * size, ((c + 1) % cliques) * size + size - 1, 1.0));
    }
    return g;
}

static double brute_modularity(const MatrixGraph *g, const int *comm)
{
    const size_t n = matrix_graph_id_bound(g);
    double A[MAXN * MAXN] = { 0 }, m2 = 0, k[MAXN] = { 0 }, q = 0;
    for (size_t u = 0; u < n; ++u) {
        for (size_t v = 0; v < n; ++v) {
            if (u == v) continue;
            const double w = wt(g, (int)u, (int)v) + (matrix_graph_is_directed(g) ? wt(g, (int)v, (int)u) : 0.0);
            A[u * n + v] = matrix_graph_is_directed(g) ? w : wt(g, (int)u, (int)v);
        }
    }
    for (size_t u = 0; u < n; ++u) {
        for (size_t v = 0; v < n; ++v) {
            k[u] += A[u * n + v];
            m2   += A[u * n + v];
        }
    }
    if (m2 == 0.0) return 0.0;
    for (size_t u = 0; u < n; ++u) {
        for (size_t v = 0; v < n; ++v) {
            if (matrix_graph_is_alive(g, (int)u) && matrix_graph_is_alive(g, (int)v) && comm[u] == comm[v]) {
                q += A[u * n + v] - k[u] * k[v] / m2;
            }
        }
    }
    return q / m2;
}

static void test_communities(void)
{
    MatrixGraph *ring = ring_of_cliques(4, 5);
    int lp[20], lv[20];
    size_t nlp, nlv;
    double q;
    CHECK_OK(matrix_label_propagation(ring, lp, &nlp));
    for (int c = 0; c < 4; ++c) {                       /* each clique is entirely one label */
        for (int a = 1; a < 5; ++a) {
            CHECK(lp[c * 5] == lp[c * 5 + a]);
        }
    }
    CHECK(nlp >= 1 && nlp <= 4);
    CHECK_OK(matrix_louvain(ring, lv, &q, &nlv));
    CHECK(nlv == 4 && q > 0.6);
    for (int c = 0; c < 4; ++c) {
        for (int a = 1; a < 5; ++a) {
            CHECK(lv[c * 5] == lv[c * 5 + a]);
        }
        CHECK(lv[c * 5] != lv[((c + 1) % 4) * 5]);
    }
    double q2;
    CHECK_OK(matrix_modularity(ring, lv, &q2));
    CHECK(near(q, q2) && near(q2, brute_modularity(ring, lv)));
    int singles[20];
    for (int i = 0; i < 20; ++i) singles[i] = i;
    double qs;
    CHECK_OK(matrix_modularity(ring, singles, &qs));
    CHECK(q > qs);
    int one[20] = { 0 };
    CHECK_OK(matrix_modularity(ring, one, &qs));
    CHECK(fabs(qs) < 1e-12);                             /* everything together: Q = 0 */
    matrix_graph_destroy(ring);

    for (int seed = 1; seed <= 60; ++seed) {
        for (int directed = 0; directed < 2; ++directed) {
            MatrixGraph *g = random_graph(14, 0.2, directed, 4, (unsigned)seed, (size_t)(seed % 3));
            const size_t n = matrix_graph_id_bound(g);
            int comm[14];
            size_t nc;
            CHECK_OK(matrix_louvain(g, comm, &q, &nc));
            CHECK_OK(matrix_modularity(g, comm, &q2));
            CHECK(near(q, q2) && fabs(q2 - brute_modularity(g, comm)) < 1e-9);
            int sing[14];
            for (size_t v = 0; v < n; ++v) {
                sing[v] = (int)v;
            }
            CHECK_OK(matrix_modularity(g, sing, &qs));
            CHECK(q >= qs - 1e-12);                      /* never worse than singletons */
            size_t seen = 0;
            bool used[14] = { false };
            for (size_t v = 0; v < n; ++v) {
                if (!matrix_graph_is_alive(g, (int)v)) {
                    CHECK(comm[v] == -1);
                    continue;
                }
                CHECK(comm[v] >= 0 && (size_t)comm[v] < nc);
                if (!used[comm[v]]) {
                    used[comm[v]] = true;
                    ++seen;
                }
            }
            CHECK(seen == nc);
            /* Deterministic. */
            int again[14];
            CHECK_OK(matrix_louvain(g, again, NULL, NULL));
            CHECK(memcmp(comm, again, n * sizeof(int)) == 0);
            CHECK_OK(matrix_label_propagation(g, comm, &nc));
            CHECK_OK(matrix_label_propagation(g, again, NULL));
            CHECK(memcmp(comm, again, n * sizeof(int)) == 0);
            matrix_graph_destroy(g);
        }
    }
    MatrixGraph *none = matrix_graph_create(3, false);
    CHECK_OK(matrix_louvain(none, lv, &q, &nlv));
    CHECK(nlv == 3 && q == 0.0);
    matrix_graph_destroy(none);
    MatrixGraph *w = matrix_gen_path(3, false);
    CHECK(matrix_graph_add_edge(w, 0, 2, -1.0));
    CHECK(matrix_louvain(w, lv, NULL, NULL) == GRAPH_ERR_INVALID_ARGUMENT);
    CHECK(matrix_label_propagation(w, lv, NULL) == GRAPH_ERR_INVALID_ARGUMENT);
    int bad[3] = { -1, 0, 0 };
    CHECK(matrix_modularity(w, bad, &q) == GRAPH_ERR_INVALID_ARGUMENT);
    matrix_graph_destroy(w);
}

/* ---- assortativity ---- */

static void test_assortativity(void)
{
    double r;
    MatrixGraph *star = matrix_gen_star(6);
    CHECK_OK(matrix_degree_assortativity(star, &r));
    CHECK(near(r, -1.0));
    matrix_graph_destroy(star);

    MatrixGraph *cyc = matrix_gen_cycle(6, false);                 /* regular: undefined */
    CHECK(matrix_degree_assortativity(cyc, &r) == GRAPH_ERR_UNSUPPORTED);
    matrix_graph_destroy(cyc);
    MatrixGraph *none = matrix_graph_create(4, false);
    CHECK(matrix_degree_assortativity(none, &r) == GRAPH_ERR_UNSUPPORTED);
    matrix_graph_destroy(none);

    for (int seed = 1; seed <= 50; ++seed) {
        MatrixGraph *g = random_graph(12, 0.3, false, 1, (unsigned)seed, (size_t)(seed % 3));
        const size_t n = matrix_graph_id_bound(g);
        int deg[12] = { 0 };
        for (size_t u = 0; u < n; ++u) {
            for (size_t v = 0; v < n; ++v) {
                deg[u] += u != v && matrix_graph_has_edge(g, (int)u, (int)v, NULL);
            }
        }
        /* Pearson correlation of the degrees at the two ends, over both orientations of each edge. */
        double cnt = 0, sx = 0, sxx = 0, sxy = 0;
        for (size_t u = 0; u < n; ++u) {
            for (size_t v = 0; v < n; ++v) {
                if (u != v && matrix_graph_has_edge(g, (int)u, (int)v, NULL)) {
                    cnt += 1;
                    sx  += deg[u];
                    sxx += (double)deg[u] * deg[u];
                    sxy += (double)deg[u] * deg[v];
                }
            }
        }
        const MatrixStatus st = matrix_degree_assortativity(g, &r);
        const double var = cnt > 0 ? sxx / cnt - (sx / cnt) * (sx / cnt) : 0.0;
        if (cnt == 0 || var <= 1e-15) {
            CHECK(st == GRAPH_ERR_UNSUPPORTED);
        } else {
            CHECK(st == GRAPH_OK);
            CHECK(fabs(r - (sxy / cnt - (sx / cnt) * (sx / cnt)) / var) < 1e-9);
            CHECK(r >= -1.0 - 1e-9 && r <= 1.0 + 1e-9);
        }
        matrix_graph_destroy(g);
    }
}

int main(void)
{
    test_pagerank();
    test_hits_eigenvector();
    test_cores();
    test_communities();
    test_assortativity();
    printf("All matrix network tests passed.\n");
    return 0;
}
