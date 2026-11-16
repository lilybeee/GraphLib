/*
 * test_util.h — helpers shared by the list test suites.
 *
 * Checks use CHECK() rather than assert() so they still run if someone builds
 * with -DNDEBUG.  The oracles here are deliberately naive (dense matrices, brute
 * force) and use only the public ListGraph API, so they do not share code, or
 * bugs, with the algorithms under test.
 */

#ifndef TEST_UTIL_LIST_H
#define TEST_UTIL_LIST_H

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "list_algorithms.h"

#define CHECK(cond)                                                                  \
    do {                                                                             \
        if (!(cond)) {                                                               \
            fprintf(stderr, "%s:%d: CHECK failed: %s\n", __FILE__, __LINE__, #cond); \
            exit(1);                                                                 \
        }                                                                            \
    } while (0)

#define CHECK_OK(expr) CHECK((expr) == GRAPH_OK)

#define HOP_INF 1000000000

static inline bool near(double a, double b)
{
    if (isinf(a) || isinf(b)) {
        return isinf(a) && isinf(b) && (a > 0) == (b > 0);
    }
    return fabs(a - b) <= 1e-9 * fmax(1.0, fmax(fabs(a), fabs(b)));
}

/* A tiny deterministic generator, independent of the library's own. */
static unsigned long long t_state = 88172645463325252ULL;

static inline void t_seed(unsigned long long s)
{
    t_state = s * 2654435761ULL + 88172645463325252ULL;
}

static inline unsigned t_rand(unsigned bound)
{
    t_state ^= t_state << 13;
    t_state ^= t_state >> 7;
    t_state ^= t_state << 17;
    return (unsigned)((t_state >> 11) % bound);
}

/** Build a graph from {u, v, w} triples. */
static inline ListGraph *build(size_t n, bool directed, const double *triples, size_t m)
{
    ListGraph *g = list_graph_create(n, directed);
    CHECK(g != NULL);
    for (size_t i = 0; i < m; ++i) {
        CHECK(list_graph_add_edge(g, (int)triples[3 * i], (int)triples[3 * i + 1], triples[3 * i + 2]));
    }
    return g;
}

/** Random graph with integer weights in [1, wmax]; random vertices removed afterwards. */
static inline ListGraph *random_graph(size_t n, double p, bool directed, int wmax, unsigned seed,
                                        size_t remove_count)
{
    t_seed(seed);
    ListGraph *g = list_graph_create(n, directed);
    CHECK(g != NULL);
    for (size_t u = 0; u < n; ++u) {
        for (size_t v = directed ? 0 : u + 1; v < n; ++v) {
            if (u != v && (double)t_rand(10000) / 10000.0 < p) {
                CHECK(list_graph_add_edge(g, (int)u, (int)v, (double)(1 + t_rand((unsigned)wmax))));
            }
        }
    }
    for (size_t i = 0; i < remove_count; ++i) {
        (void)list_graph_remove_vertex(g, (int)t_rand((unsigned)n));
    }
    return g;
}

/** Same ID space, same alive set, same edges and weights. */
static inline bool same_graph(const ListGraph *a, const ListGraph *b)
{
    if (list_graph_is_directed(a) != list_graph_is_directed(b) ||
        list_graph_id_bound(a) != list_graph_id_bound(b) ||
        list_graph_vertex_count(a) != list_graph_vertex_count(b) ||
        list_graph_edge_count(a) != list_graph_edge_count(b)) {
        return false;
    }
    const size_t n = list_graph_id_bound(a);
    for (size_t u = 0; u < n; ++u) {
        if (list_graph_is_alive(a, (int)u) != list_graph_is_alive(b, (int)u)) {
            return false;
        }
        for (size_t v = 0; v < n; ++v) {
            double wa = 0, wb = 0;
            const bool ha = list_graph_has_edge(a, (int)u, (int)v, &wa);
            const bool hb = list_graph_has_edge(b, (int)u, (int)v, &wb);
            if (ha != hb || (ha && wa != wb)) {
                return false;
            }
        }
    }
    return true;
}

/**
 * A dense copy of the subgraph of g that survives when `drop[v]` is set,
 * renumbered 0..k-1 in ID order; map[old] = new or -1.  Built from scratch with
 * only add_edge, so it is an oracle for "the algorithm ignores dead vertices".
 */
static inline ListGraph *survivors(const ListGraph *g, const bool *drop, int *map)
{
    const size_t n = list_graph_id_bound(g);
    size_t k = 0;
    for (size_t v = 0; v < n; ++v) {
        map[v] = (list_graph_is_alive(g, (int)v) && (drop == NULL || !drop[v])) ? (int)k++ : -1;
    }
    ListGraph *r = list_graph_create(k, list_graph_is_directed(g));
    CHECK(r != NULL);
    for (size_t u = 0; u < n; ++u) {
        for (size_t v = 0; v < n; ++v) {
            double w;
            if (map[u] >= 0 && map[v] >= 0 && list_graph_has_edge(g, (int)u, (int)v, &w) &&
                !list_graph_has_edge(r, map[u], map[v], NULL)) {
                CHECK(list_graph_add_edge(r, map[u], map[v], w));
            }
        }
    }
    return r;
}

/** All-pairs hop distances by Floyd-Warshall over has_edge; D is n*n ints. */
static inline void hops(const ListGraph *g, int *D)
{
    const size_t n = list_graph_id_bound(g);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            D[i * n + j] = HOP_INF;
            if (list_graph_is_alive(g, (int)i) && list_graph_is_alive(g, (int)j)) {
                if (i == j) {
                    D[i * n + j] = 0;
                } else if (list_graph_has_edge(g, (int)i, (int)j, NULL)) {
                    D[i * n + j] = 1;
                }
            }
        }
    }
    for (size_t k = 0; k < n; ++k) {
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                if (D[i * n + k] < HOP_INF && D[k * n + j] < HOP_INF &&
                    D[i * n + k] + D[k * n + j] < D[i * n + j]) {
                    D[i * n + j] = D[i * n + k] + D[k * n + j];
                }
            }
        }
    }
}

/** Number of weakly connected components, by flood fill over has_edge. */
static inline size_t naive_components(const ListGraph *g)
{
    const size_t n = list_graph_id_bound(g);
    bool *seen = calloc(n + 1, sizeof(bool));
    int  *stack = malloc((n + 1) * sizeof(int));
    CHECK(seen != NULL && stack != NULL);
    size_t comps = 0;
    for (size_t s = 0; s < n; ++s) {
        if (!list_graph_is_alive(g, (int)s) || seen[s]) {
            continue;
        }
        ++comps;
        size_t sp = 0;
        stack[sp++] = (int)s;
        seen[s] = true;
        while (sp > 0) {
            const int u = stack[--sp];
            for (size_t v = 0; v < n; ++v) {
                if (!seen[v] && list_graph_is_alive(g, (int)v) &&
                    (list_graph_has_edge(g, u, (int)v, NULL) ||
                     list_graph_has_edge(g, (int)v, u, NULL))) {
                    seen[v] = true;
                    stack[sp++] = (int)v;
                }
            }
        }
    }
    free(seen);
    free(stack);
    return comps;
}

/** Edge weight or 0 if absent (for oracles that sum weights). */
static inline double wt(const ListGraph *g, int u, int v)
{
    double w = 0.0;
    return list_graph_has_edge(g, u, v, &w) ? w : 0.0;
}

#endif /* TEST_UTIL_LIST_H */
