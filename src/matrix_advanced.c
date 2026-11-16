/*
 * matrix_advanced.c — twenty-five further algorithms on the adjacency-matrix
 * backend.  See matrix_advanced.h.
 *
 * Reuse
 * -----
 *   Everything is written against the public MatrixGraph API (is_alive,
 *   id_bound, next_neighbor, has_edge, add_edge, ...) and the containers in
 *   matrix_util.h (MuHeap, MuSym, bitsets).  Where the backend already has an
 *   algorithm it is called, not rewritten: matrix_topological_sort,
 *   matrix_has_cycle, matrix_bfs, matrix_dijkstra and the degree / closeness /
 *   betweenness centralities of matrix_analysis.c.
 *
 * Conventions
 * -----------
 *   - Results are MatrixStatus values, as in matrix_analysis.c.
 *   - Caller-supplied arrays are indexed by vertex ID and sized
 *     matrix_graph_id_bound(); dead IDs receive a documented sentinel.
 *   - Dead IDs are skipped everywhere; the graph is never modified.
 *   - Complexity is stated per function.  V = alive vertices, E = edges as
 *     counted by matrix_graph_edge_count, N = id_bound.
 *
 * Compile: cc -std=c11 -Wall -Wextra -Werror -pedantic -c matrix_advanced.c
 */

#include "matrix_advanced.h"
#include "matrix_util.h"

#include <limits.h>

#define MA_PI        3.14159265358979323846
#define MA_SQRT2     1.41421356237309504880
#define MA_CC_MAX_K  20

/* Negative-cycle failures are reported with this code;
   all uses go through this macro so it can be repointed in one place. */
#define MA_ERR_NEG_CYCLE GRAPH_ERR_NEGATIVE_CYCLE

/* ------------------------------------------------------------------ */
/*  Small shared helpers                                                */
/* ------------------------------------------------------------------ */

/** malloc of n elements; never NULL for n == 0, NULL on overflow. */
static void *ma_malloc(size_t n, size_t size)
{
    const size_t bytes = mu_mul(n > 0 ? n : 1, size);
    return bytes == 0 ? NULL : malloc(bytes);
}

/** Grows buf to hold at least `need` elements.  NULL on failure (buf intact). */
static void *ma_grow(void *buf, size_t *cap, size_t need, size_t elem)
{
    if (need <= *cap) {
        return buf;
    }
    size_t c = *cap != 0 ? *cap : 64;
    while (c < need) {
        if (c > SIZE_MAX / 2) {
            return NULL;
        }
        c *= 2;
    }
    if (c > SIZE_MAX / elem) {
        return NULL;
    }
    void *p = realloc(buf, c * elem);
    if (p != NULL) {
        *cap = c;
    }
    return p;
}

static inline int ma_ctz64(uint64_t x)
{
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_ctzll(x);
#else
    int n = 0;
    while ((x & 1) == 0) {
        x >>= 1;
        ++n;
    }
    return n;
#endif
}

static int ma_cmp_int_desc(const void *l, const void *r)
{
    const int a = *(const int *)l;
    const int b = *(const int *)r;
    return (a < b) - (a > b);
}

static int ma_cmp_int_asc(const void *l, const void *r)
{
    const int a = *(const int *)l;
    const int b = *(const int *)r;
    return (a > b) - (a < b);
}

/** Comparison with context (C11 has no qsort_r): <0, 0, >0 on items a and b. */
typedef int (*MaCmp)(const void *ctx, int a, int b);

/** Stable bottom-up merge sort of idx[0..n).  Returns false if out of memory. */
static bool ma_sort_idx(int *idx, size_t n, MaCmp cmp, const void *ctx)
{
    if (n < 2) {
        return true;
    }
    int *tmp = ma_malloc(n, sizeof(*tmp));
    if (tmp == NULL) {
        return false;
    }
    for (size_t width = 1; width < n; width *= 2) {
        for (size_t lo = 0; lo < n; lo += 2 * width) {
            const size_t mid = lo + width < n ? lo + width : n;
            const size_t hi  = lo + 2 * width < n ? lo + 2 * width : n;
            size_t i = lo, j = mid, k = lo;
            while (i < mid && j < hi) {
                tmp[k++] = cmp(ctx, idx[j], idx[i]) < 0 ? idx[j++] : idx[i++];
            }
            while (i < mid) {
                tmp[k++] = idx[i++];
            }
            while (j < hi) {
                tmp[k++] = idx[j++];
            }
        }
        memcpy(idx, tmp, n * sizeof(*idx));
    }
    free(tmp);
    return true;
}

/**
 * ma_path_from_parents - build a MatrixAdvPath from a parent array.
 * NOT_FOUND if target is not linked back to source.  O(path length).
 */
static MatrixStatus ma_path_from_parents(const int *parent, size_t n, int source, int target,
                                         double cost, MatrixAdvPath *out)
{
    size_t len = 1;
    int    v   = target;
    while (v != source) {
        v = parent[v];
        if (v < 0 || (size_t)v >= n || ++len > n) {
            return GRAPH_ERR_NOT_FOUND;
        }
    }
    int *vs = malloc(len * sizeof(*vs));
    if (vs == NULL) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    v = target;
    for (size_t i = len; i-- > 0;) {
        vs[i] = v;
        v     = parent[v];
    }
    out->vertices = vs;
    out->length   = len;
    out->cost     = cost;
    return GRAPH_OK;
}

/**
 * ma_clone - rebuild g with the same IDs.  If pot is non-NULL every weight is
 * replaced by w + pot[u] - pot[v], clamped at 0 (Johnson's reweighting).
 * The graph API has no weight setter, so the copy is built edge by edge, the
 * same way matrix_transpose does.
 */
static MatrixGraph *ma_clone(const MatrixGraph *g, const double *pot)
{
    const size_t n        = matrix_graph_id_bound(g);
    const bool   directed = matrix_graph_is_directed(g);
    MatrixGraph *c        = matrix_graph_create(n, directed);
    if (c == NULL) {
        return NULL;
    }
    for (size_t u = 0; u < n; ++u) {
        if (!matrix_graph_is_alive(g, (int)u) && !matrix_graph_remove_vertex(c, (int)u)) {
            matrix_graph_destroy(c);
            return NULL;
        }
    }
    for (size_t ui = 0; ui < n; ++ui) {
        const int u = (int)ui;
        if (!matrix_graph_is_alive(g, u)) {
            continue;
        }
        size_t cur = 0;
        int    v;
        double w;
        while (matrix_graph_next_neighbor(g, u, &cur, &v, &w)) {
            if (!directed && v < u) {
                continue;               /* the mirror cell is added with (v, u) */
            }
            if (pot != NULL) {
                w += pot[u] - pot[v];
                if (w < 0.0) {
                    w = 0.0;
                }
            }
            if (!matrix_graph_add_edge(c, u, v, w)) {
                matrix_graph_destroy(c);
                return NULL;
            }
        }
    }
    return c;
}

/**
 * ma_potentials - Bellman-Ford from a virtual source joined to every alive
 * vertex by a 0-weight arc.  pot[v] <= 0 is the shortest distance, par[] the
 * tree.  Returns MA_ERR_NEG_CYCLE (with *witness set) if a negative cycle is
 * reachable.  A relaxation must improve by more than eps.  O(V * E).
 */
static MatrixStatus ma_potentials(const MatrixGraph *g, double eps, double *pot, int *par,
                                  int *witness)
{
    const size_t n    = matrix_graph_id_bound(g);
    const size_t live = matrix_graph_vertex_count(g);
    mu_fill_double(pot, n, 0.0);
    mu_fill_int(par, n, -1);
    *witness = -1;
    for (size_t round = 0; round < live; ++round) {
        bool changed = false;
        int  last    = -1;
        for (size_t ui = 0; ui < n; ++ui) {
            const int u = (int)ui;
            if (!matrix_graph_is_alive(g, u)) {
                continue;
            }
            size_t cur = 0;
            int    v;
            double w;
            while (matrix_graph_next_neighbor(g, u, &cur, &v, &w)) {
                if (pot[u] + w < pot[v] - eps) {
                    pot[v]  = pot[u] + w;
                    par[v]  = u;
                    changed = true;
                    last    = v;
                }
            }
        }
        if (!changed) {
            return GRAPH_OK;
        }
        *witness = last;
    }
    return live == 0 ? GRAPH_OK : MA_ERR_NEG_CYCLE;
}

/** Walks parents from a vertex that was relaxed in the last round into the
 *  cycle and writes it in arc order.  Returns its length (0 on failure). */
static size_t ma_extract_cycle(const int *par, size_t n, int witness, int *cycle)
{
    int v = witness;
    for (size_t i = 0; i < n; ++i) {
        if (v < 0 || par[v] < 0) {
            return 0;
        }
        v = par[v];
    }
    size_t len = 0;
    int    x   = v;
    do {
        cycle[len++] = x;
        x = par[x];
        if (x < 0 || len > n) {
            return 0;
        }
    } while (x != v);
    for (size_t i = 0; i < len / 2; ++i) {
        const int t = cycle[i];
        cycle[i] = cycle[len - 1 - i];
        cycle[len - 1 - i] = t;
    }
    return len;
}

/** Fails with GRAPH_ERR_NEGATIVE_WEIGHT if any weight reachable by iteration is < 0. */
static bool ma_has_negative_weight(const MatrixGraph *g)
{
    const size_t n = matrix_graph_id_bound(g);
    for (size_t ui = 0; ui < n; ++ui) {
        const int u = (int)ui;
        if (!matrix_graph_is_alive(g, u)) {
            continue;
        }
        size_t cur = 0;
        int    v;
        double w;
        while (matrix_graph_next_neighbor(g, u, &cur, &v, &w)) {
            if (w < 0.0) {
                return true;
            }
        }
    }
    return false;
}

/* ------------------------------------------------------------------ */
/*  Result helpers                                                      */
/* ------------------------------------------------------------------ */

void matrix_adv_path_free(MatrixAdvPath *path)
{
    if (path == NULL) {
        return;
    }
    free(path->vertices);
    memset(path, 0, sizeof(*path));
}

void matrix_adv_cycles_free(MatrixAdvCycles *cycles)
{
    if (cycles == NULL) {
        return;
    }
    free(cycles->offsets);
    free(cycles->members);
    memset(cycles, 0, sizeof(*cycles));
}

/**
 * matrix_adv_path_from_parents - public wrapper: path source -> target from a
 * parent[] array of n entries (as produced by matrix_dijkstra, matrix_bfs, the
 * widest-path routine, ...).  cost is stored in the result unchanged.
 * Time O(path length).
 */
MatrixStatus matrix_adv_path_from_parents(const int *parent, size_t n, int source, int target,
                                          double cost, MatrixAdvPath *out)
{
    if (parent == NULL || out == NULL || source < 0 || target < 0 || (size_t)source >= n ||
        (size_t)target >= n) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    memset(out, 0, sizeof(*out));
    return ma_path_from_parents(parent, n, source, target, cost, out);
}

/**
 * matrix_adv_apsp_path - path from an all-pairs result (Floyd-Warshall or
 * Johnson): dist and pred are n x n, row-major, pred[s*n+v] being the vertex
 * before v on the shortest path from s.  NOT_FOUND if unreachable.
 */
MatrixStatus matrix_adv_apsp_path(const double *dist, const int *pred, size_t n,
                                  int source, int target, MatrixAdvPath *out)
{
    if (dist == NULL || pred == NULL || out == NULL || source < 0 || target < 0 ||
        (size_t)source >= n || (size_t)target >= n) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    memset(out, 0, sizeof(*out));
    const double d = dist[(size_t)source * n + (size_t)target];
    if (!isfinite(d)) {
        return GRAPH_ERR_NOT_FOUND;
    }
    return ma_path_from_parents(pred + (size_t)source * n, n, source, target, d, out);
}

/* ================================================================== */
/*  1. Longest path (DAG)                                               */
/* ================================================================== */

/**
 * matrix_adv_longest_path - longest weighted path from source in a directed
 * acyclic graph, by dynamic programming over the topological order given by
 * matrix_topological_sort.  Weights may be negative.
 *
 * dist[v] is the longest distance (-INFINITY if unreachable or dead);
 * parent[v] the predecessor (-1 for the source, unreachable and dead IDs).
 * parent may be NULL.
 *
 * Errors: UNSUPPORTED for undirected graphs, CYCLE if the graph has a cycle
 * (the problem is NP-hard there), INVALID_VERTEX for a bad source.
 * Time O(V + E) on top of the sort; the sort is O(N * N/64 + V + E).
 */
MatrixStatus matrix_adv_longest_path(const MatrixGraph *g, int source, double *dist, int *parent)
{
    if (g == NULL || dist == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    if (!matrix_graph_is_alive(g, source)) {
        return GRAPH_ERR_INVALID_VERTEX;
    }
    if (!matrix_graph_is_directed(g)) {
        return GRAPH_ERR_UNSUPPORTED;
    }
    const size_t n    = matrix_graph_id_bound(g);
    const size_t live = matrix_graph_vertex_count(g);
    int *order = ma_malloc(n, sizeof(*order));
    if (order == NULL) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    if (!matrix_topological_sort(g, order)) {
        free(order);
        return matrix_has_cycle(g) ? GRAPH_ERR_CYCLE : GRAPH_ERR_OUT_OF_MEMORY;
    }
    mu_fill_double(dist, n, -INFINITY);
    if (parent != NULL) {
        mu_fill_int(parent, n, -1);
    }
    dist[source] = 0.0;
    for (size_t i = 0; i < live; ++i) {
        const int u = order[i];
        if (dist[u] == -INFINITY) {
            continue;
        }
        size_t cur = 0;
        int    v;
        double w;
        while (matrix_graph_next_neighbor(g, u, &cur, &v, &w)) {
            if (dist[u] + w > dist[v]) {
                dist[v] = dist[u] + w;
                if (parent != NULL) {
                    parent[v] = u;
                }
            }
        }
    }
    free(order);
    return GRAPH_OK;
}

/* ================================================================== */
/*  2. PageRank                                                         */
/* ================================================================== */

/**
 * matrix_adv_pagerank - power iteration
 *     r' = (1-d)/V + d * (sum over in-arcs r[u]*share(u,v) + dangling/V).
 * Out-arcs of u are weighted by w/strength(u) if `weighted` (weights must be
 * >= 0), else uniformly.  Vertices with no outgoing weight are dangling and
 * spread their rank evenly over all alive vertices.  Undirected edges act as
 * arcs in both directions.  Self-loops count as ordinary arcs.
 *
 * rank[v] sums to 1 over alive vertices; dead IDs get 0.  Stops when the L1
 * change drops below tol (tol == 0 runs max_iter rounds).  *iterations
 * (optional) receives the number of rounds performed.
 * Errors: INVALID_ARGUMENT for damping outside [0,1), negative tol or
 * max_iter == 0; NEGATIVE_WEIGHT.  Time O(max_iter * (V + E)) over the bit rows.
 */
MatrixStatus matrix_adv_pagerank(const MatrixGraph *g, double damping, double tol,
                                 size_t max_iter, bool weighted, double *rank,
                                 size_t *iterations)
{
    if (g == NULL || rank == NULL || !(damping >= 0.0 && damping < 1.0) || !(tol >= 0.0) ||
        max_iter == 0) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n  = matrix_graph_id_bound(g);
    const size_t vc = matrix_graph_vertex_count(g);
    if (iterations != NULL) {
        *iterations = 0;
    }
    mu_fill_double(rank, n, 0.0);
    if (vc == 0) {
        return GRAPH_OK;
    }
    double *strength = ma_malloc(n, sizeof(*strength));
    double *next     = ma_malloc(n, sizeof(*next));
    if (strength == NULL || next == NULL) {
        free(strength);
        free(next);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    for (size_t ui = 0; ui < n; ++ui) {
        const int u = (int)ui;
        strength[ui] = 0.0;
        if (!matrix_graph_is_alive(g, u)) {
            continue;
        }
        size_t cur = 0;
        int    v;
        double w;
        while (matrix_graph_next_neighbor(g, u, &cur, &v, &w)) {
            if (weighted && w < 0.0) {
                free(strength);
                free(next);
                return GRAPH_ERR_NEGATIVE_WEIGHT;
            }
            strength[ui] += weighted ? w : 1.0;
        }
        rank[ui] = 1.0 / (double)vc;
    }

    size_t it = 0;
    while (it < max_iter) {
        double dangling = 0.0;
        for (size_t ui = 0; ui < n; ++ui) {
            if (matrix_graph_is_alive(g, (int)ui) && strength[ui] <= 0.0) {
                dangling += rank[ui];
            }
        }
        const double base = ((1.0 - damping) + damping * dangling) / (double)vc;
        for (size_t ui = 0; ui < n; ++ui) {
            next[ui] = matrix_graph_is_alive(g, (int)ui) ? base : 0.0;
        }
        for (size_t ui = 0; ui < n; ++ui) {
            const int u = (int)ui;
            if (!matrix_graph_is_alive(g, u) || strength[ui] <= 0.0) {
                continue;
            }
            const double scale = damping * rank[ui] / strength[ui];
            size_t cur = 0;
            int    v;
            double w;
            while (matrix_graph_next_neighbor(g, u, &cur, &v, &w)) {
                next[v] += scale * (weighted ? w : 1.0);
            }
        }
        double diff = 0.0;
        for (size_t ui = 0; ui < n; ++ui) {
            diff += fabs(next[ui] - rank[ui]);
            rank[ui] = next[ui];
        }
        ++it;
        if (diff < tol) {
            break;
        }
    }
    if (iterations != NULL) {
        *iterations = it;
    }
    free(strength);
    free(next);
    return GRAPH_OK;
}

/* ================================================================== */
/*  3. Reverse-delete minimum spanning tree                             */
/* ================================================================== */

static int ma_cmp_edge_desc(const void *left, const void *right)
{
    const MatrixEdge *a = left;
    const MatrixEdge *b = right;
    if (a->weight > b->weight) {
        return -1;
    }
    if (a->weight < b->weight) {
        return 1;
    }
    if (a->u != b->u) {
        return a->u > b->u ? -1 : 1;
    }
    if (a->v != b->v) {
        return a->v > b->v ? -1 : 1;
    }
    return 0;
}

/**
 * matrix_adv_reverse_delete - minimum spanning forest of an undirected graph.
 * Edges are examined from heaviest to lightest and deleted from a private
 * working copy unless the deletion disconnects their endpoints (tested with
 * matrix_bfs).  The input graph is not modified.
 *
 * edges must have room for vertex_count - 1 entries; they come back sorted by
 * ascending weight (ties by endpoints).  A disconnected graph yields a
 * spanning forest rather than an error; *components (optional) receives its
 * number of trees.  *total_weight (optional) is the forest's weight.
 * Self-loops are never kept.  Negative weights are fine.
 * Errors: UNSUPPORTED for directed graphs.
 * Time O(E * (N * N/64 + E/64)) — one BFS per edge; intended for modest graphs.
 */
MatrixStatus matrix_adv_reverse_delete(const MatrixGraph *g, MatrixEdge *edges,
                                       size_t *edge_count, double *total_weight,
                                       size_t *components)
{
    if (g == NULL || edge_count == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    *edge_count = 0;
    if (total_weight != NULL) {
        *total_weight = 0.0;
    }
    if (matrix_graph_is_directed(g)) {
        return GRAPH_ERR_UNSUPPORTED;
    }
    const size_t live = matrix_graph_vertex_count(g);
    if (components != NULL) {
        *components = live;
    }
    if (live <= 1) {
        return GRAPH_OK;
    }
    if (edges == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }

    const size_t n = matrix_graph_id_bound(g);
    MatrixEdge  *cand   = ma_malloc(matrix_graph_edge_count(g), sizeof(*cand));
    int         *order  = ma_malloc(n, sizeof(*order));
    int         *parent = ma_malloc(n, sizeof(*parent));
    MatrixGraph *work   = ma_clone(g, NULL);
    MatrixStatus status = GRAPH_OK;
    size_t       count  = 0;
    size_t       kept   = 0;
    double       total  = 0.0;

    if (cand == NULL || order == NULL || parent == NULL || work == NULL) {
        status = GRAPH_ERR_OUT_OF_MEMORY;
    } else {
        for (size_t ui = 0; ui < n; ++ui) {
            const int u = (int)ui;
            size_t cur = 0;
            int    v;
            double w;
            if (!matrix_graph_is_alive(g, u)) {
                continue;
            }
            while (matrix_graph_next_neighbor(g, u, &cur, &v, &w)) {
                if (v > u && count < matrix_graph_edge_count(g)) {
                    cand[count++] = (MatrixEdge){u, v, w};
                }
            }
        }
        qsort(cand, count, sizeof(*cand), ma_cmp_edge_desc);

        for (size_t i = 0; i < count && status == GRAPH_OK; ++i) {
            const MatrixEdge e = cand[i];
            if (!matrix_graph_remove_edge(work, e.u, e.v)) {
                status = GRAPH_ERR_INVALID_ARGUMENT;
                break;
            }
            if (matrix_bfs(work, e.u, order, parent) == 0) {
                status = GRAPH_ERR_OUT_OF_MEMORY;
                break;
            }
            if (parent[e.v] == -1) {
                /* A bridge: put it back and keep it. */
                if (!matrix_graph_add_edge(work, e.u, e.v, e.weight)) {
                    status = GRAPH_ERR_OUT_OF_MEMORY;
                    break;
                }
                edges[kept++] = e;
                total += e.weight;
            }
        }
    }
    if (status == GRAPH_OK) {
        for (size_t i = 0; i < kept / 2; ++i) {      /* heaviest-first -> ascending */
            const MatrixEdge t = edges[i];
            edges[i] = edges[kept - 1 - i];
            edges[kept - 1 - i] = t;
        }
        *edge_count = kept;
        if (total_weight != NULL) {
            *total_weight = total;
        }
        if (components != NULL) {
            *components = live - kept;
        }
    }
    free(cand);
    free(order);
    free(parent);
    matrix_graph_destroy(work);
    return status;
}

/* ================================================================== */
/*  4. Rocha–Thatte cycle detection                                     */
/* ================================================================== */

static bool ma_seq_contains(const int *seq, size_t len, int v)
{
    for (size_t i = 0; i < len; ++i) {
        if (seq[i] == v) {
            return true;
        }
    }
    return false;
}

/**
 * matrix_adv_rocha_thatte - enumerate the simple cycles of a directed graph
 * with the Rocha–Thatte superstep scheme.  In superstep L every vertex holds
 * the vertex sequences of length L that end in it; it appends each out-neighbour
 * and forwards the result.  A sequence that reaches its own first vertex is a
 * cycle; one that meets another vertex it already contains is dropped.
 * Following the usual refinement, a sequence is extended only through vertices
 * with a larger ID than its first vertex, so each cycle is found exactly once,
 * starting at its smallest vertex.  Self-loops are cycles of length 1.
 *
 * max_len limits the cycle length (0 = no limit); max_cycles stops the search
 * once that many cycles are stored (0 = no limit) and sets out->truncated.
 * Results: out->offsets / out->members (see MatrixAdvCycles); free with
 * matrix_adv_cycles_free.  Errors: UNSUPPORTED for undirected graphs.
 * Time O(P * L) where P is the number of simple paths examined (exponential in
 * the worst case; the output itself can be exponential).
 */
MatrixStatus matrix_adv_rocha_thatte(const MatrixGraph *g, size_t max_len, size_t max_cycles,
                                     MatrixAdvCycles *out)
{
    if (g == NULL || out == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    memset(out, 0, sizeof(*out));
    if (!matrix_graph_is_directed(g)) {
        return GRAPH_ERR_UNSUPPORTED;
    }
    const size_t n     = matrix_graph_id_bound(g);
    const size_t live  = matrix_graph_vertex_count(g);
    const size_t limit = (max_len == 0 || max_len > live) ? live : max_len;

    int    *cur = NULL, *next = NULL;
    size_t  cur_cap = 0, next_cap = 0, cur_n = 0;
    size_t  mem_cap = 0, off_cap = 0, mem_n = 0;
    MatrixStatus status = GRAPH_OK;

    out->offsets = ma_grow(NULL, &off_cap, 1, sizeof(*out->offsets));
    if (out->offsets == NULL) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    out->offsets[0] = 0;
    cur = ma_grow(NULL, &cur_cap, live > 0 ? live : 1, sizeof(*cur));
    if (cur == NULL) {
        matrix_adv_cycles_free(out);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    for (size_t v = 0; v < n; ++v) {
        if (matrix_graph_is_alive(g, (int)v)) {
            cur[cur_n++] = (int)v;
        }
    }

    for (size_t len = 1; len <= limit && cur_n > 0 && status == GRAPH_OK; ++len) {
        size_t next_n = 0;
        for (size_t q = 0; q < cur_n && status == GRAPH_OK; ++q) {
            const int *seq = cur + q * len;
            const int  first = seq[0];
            size_t cursor = 0;
            int    v;
            while (matrix_graph_next_neighbor(g, seq[len - 1], &cursor, &v, NULL)) {
                if (v == first) {
                    int *m = ma_grow(out->members, &mem_cap, mem_n + len, sizeof(*m));
                    if (m == NULL) {
                        status = GRAPH_ERR_OUT_OF_MEMORY;
                        break;
                    }
                    out->members = m;
                    size_t *o = ma_grow(out->offsets, &off_cap, out->count + 2, sizeof(*o));
                    if (o == NULL) {
                        status = GRAPH_ERR_OUT_OF_MEMORY;
                        break;
                    }
                    out->offsets = o;
                    memcpy(out->members + mem_n, seq, len * sizeof(*seq));
                    mem_n += len;
                    out->offsets[++out->count] = mem_n;
                    if (max_cycles != 0 && out->count >= max_cycles) {
                        out->truncated = true;
                        status = GRAPH_OK;
                        q = cur_n;        /* leave both loops */
                        len = limit;
                        break;
                    }
                } else if (len < limit && v > first && !ma_seq_contains(seq, len, v)) {
                    const size_t need = (next_n + 1) * (len + 1);
                    int *p = ma_grow(next, &next_cap, need, sizeof(*p));
                    if (p == NULL) {
                        status = GRAPH_ERR_OUT_OF_MEMORY;
                        break;
                    }
                    next = p;
                    memcpy(next + next_n * (len + 1), seq, len * sizeof(*seq));
                    next[next_n * (len + 1) + len] = v;
                    ++next_n;
                }
            }
        }
        int *swap = cur;       /* the new sequences become the current ones */
        cur = next;
        next = swap;
        const size_t swap_cap = cur_cap;
        cur_cap = next_cap;
        next_cap = swap_cap;
        cur_n = next_n;
    }
    free(cur);
    free(next);
    if (status != GRAPH_OK) {
        matrix_adv_cycles_free(out);
    }
    return status;
}

/* ================================================================== */
/*  5. Tree traversal                                                   */
/* ================================================================== */

typedef struct {
    int    v;
    int    parent;
    size_t cursor;
    size_t kids;
    bool   emitted;
} MaFrame;

/**
 * matrix_adv_tree_traversal - pre-, in-, post- or level-order listing of a
 * tree.  Directed graphs are read as rooted trees (children = out-neighbours);
 * undirected graphs are rooted at `root` (children = neighbours except the
 * parent).  Children are visited in ascending ID order.  In-order is the
 * natural generalisation of the binary case: the subtree of the first child,
 * then the node, then the subtrees of the remaining children (so a node with
 * one child treats it as the left child).
 *
 * The graph must be a tree: every alive vertex reachable from root and exactly
 * V - 1 edges.  Otherwise UNSUPPORTED.  out needs vertex_count entries;
 * *count receives that number.  Iterative (no recursion).  Time O(V * N / 64).
 */
MatrixStatus matrix_adv_tree_traversal(const MatrixGraph *g, int root, MatrixAdvTreeOrder order,
                                       int *out, size_t *count)
{
    if (g == NULL || out == NULL || count == NULL || (unsigned)order > MATRIX_ADV_TREE_LEVELORDER) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    *count = 0;
    if (!matrix_graph_is_alive(g, root)) {
        return GRAPH_ERR_INVALID_VERTEX;
    }
    const size_t vc = matrix_graph_vertex_count(g);
    if (matrix_graph_edge_count(g) + 1 != vc) {
        return GRAPH_ERR_UNSUPPORTED;
    }
    /* matrix_bfs doubles as the tree check and as the level-order listing. */
    int *bfs = ma_malloc(matrix_graph_id_bound(g), sizeof(*bfs));
    if (bfs == NULL) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    const size_t reached = matrix_bfs(g, root, bfs, NULL);
    if (reached != vc) {
        free(bfs);
        return reached == 0 ? GRAPH_ERR_OUT_OF_MEMORY : GRAPH_ERR_UNSUPPORTED;
    }
    if (order == MATRIX_ADV_TREE_LEVELORDER) {
        memcpy(out, bfs, vc * sizeof(*out));
        free(bfs);
        *count = vc;
        return GRAPH_OK;
    }
    free(bfs);

    MaFrame *stack = ma_malloc(vc, sizeof(*stack));
    if (stack == NULL) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    size_t top = 0, k = 0;
    stack[0] = (MaFrame){root, -1, 0, 0, false};
    if (order == MATRIX_ADV_TREE_PREORDER) {
        out[k++] = root;
    }
    for (;;) {
        MaFrame *f = &stack[top];
        if (order == MATRIX_ADV_TREE_INORDER && f->kids >= 1 && !f->emitted) {
            out[k++]   = f->v;      /* the first child's subtree is finished */
            f->emitted = true;
        }
        int  c;
        bool found = false;
        while (matrix_graph_next_neighbor(g, f->v, &f->cursor, &c, NULL)) {
            if (c != f->parent) {
                found = true;
                break;
            }
        }
        if (found) {
            ++f->kids;
            stack[top + 1] = (MaFrame){c, f->v, 0, 0, false};
            ++top;
            if (order == MATRIX_ADV_TREE_PREORDER) {
                out[k++] = c;
            }
            continue;
        }
        if (order == MATRIX_ADV_TREE_INORDER && !f->emitted) {
            out[k++] = f->v;        /* a leaf */
        } else if (order == MATRIX_ADV_TREE_POSTORDER) {
            out[k++] = f->v;
        }
        if (top == 0) {
            break;
        }
        --top;
    }
    free(stack);
    *count = k;
    return GRAPH_OK;
}

/* ================================================================== */
/*  6. Zero-weight cycle                                                */
/* ================================================================== */

/**
 * matrix_adv_zero_weight_cycle - decide whether a directed graph has a cycle
 * of weight exactly 0 (|weight| <= eps) and return one.
 *
 * Method: Bellman-Ford potentials pi from a virtual source give every arc a
 * reduced cost w + pi(u) - pi(v) >= 0, and a cycle's reduced cost equals its
 * weight.  So when no negative cycle exists, a zero-weight cycle exists iff
 * the "tight" arcs (reduced cost 0) contain a cycle — found by an iterative
 * DFS.
 *
 * *kind: NONE (all cycles positive, or no cycle), ZERO (cycle returned), or
 * NEGATIVE (a negative cycle exists; the potentials do not exist, a negative
 * cycle is returned instead).  A zero cycle that coexists with a negative one
 * is not searched for: with real weights that is a subset-sum style problem.
 * cycle needs vertex_count entries (may be NULL); *cycle_len (may be NULL)
 * receives its length.  eps >= 0 is the tolerance for both the relaxation and
 * the tightness test.  Errors: UNSUPPORTED for undirected graphs.
 * Time O(V * E) for the potentials, O(V + E) for the DFS.
 */
MatrixStatus matrix_adv_zero_weight_cycle(const MatrixGraph *g, double eps,
                                          MatrixAdvCycleKind *kind, int *cycle,
                                          size_t *cycle_len)
{
    if (g == NULL || kind == NULL || !(eps >= 0.0)) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    *kind = MATRIX_ADV_CYCLE_NONE;
    if (cycle_len != NULL) {
        *cycle_len = 0;
    }
    if (!matrix_graph_is_directed(g)) {
        return GRAPH_ERR_UNSUPPORTED;
    }
    const size_t n = matrix_graph_id_bound(g);
    double *pot   = ma_malloc(n, sizeof(*pot));
    int    *par   = ma_malloc(n, sizeof(*par));
    int    *stack = ma_malloc(n, sizeof(*stack));
    size_t *cursor = ma_malloc(n, sizeof(*cursor));
    unsigned char *color = mu_calloc(n, sizeof(*color));
    MatrixStatus status = GRAPH_OK;

    if (pot == NULL || par == NULL || stack == NULL || cursor == NULL || color == NULL) {
        status = GRAPH_ERR_OUT_OF_MEMORY;
    } else {
        int witness;
        const MatrixStatus st = ma_potentials(g, eps, pot, par, &witness);
        if (st == MA_ERR_NEG_CYCLE) {
            *kind = MATRIX_ADV_CYCLE_NEGATIVE;
            if (cycle != NULL && witness >= 0) {
                const size_t len = ma_extract_cycle(par, n, witness, cycle);
                if (cycle_len != NULL) {
                    *cycle_len = len;
                }
            }
        } else {
            for (size_t root = 0; root < n && *kind == MATRIX_ADV_CYCLE_NONE; ++root) {
                if (color[root] != 0 || !matrix_graph_is_alive(g, (int)root)) {
                    continue;
                }
                size_t top = 0;
                stack[0]   = (int)root;
                cursor[0]  = 0;
                color[root] = 1;
                /* cursor[] is indexed by stack depth here, not by vertex. */
                while (*kind == MATRIX_ADV_CYCLE_NONE) {
                    const int u = stack[top];
                    int    v;
                    double w;
                    bool   found = false;
                    while (matrix_graph_next_neighbor(g, u, &cursor[top], &v, &w)) {
                        if (fabs(w + pot[u] - pot[v]) <= eps) {
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        color[u] = 2;
                        if (top == 0) {
                            break;
                        }
                        --top;
                        continue;
                    }
                    if (color[v] == 0) {
                        color[v] = 1;
                        ++top;
                        stack[top]  = v;
                        cursor[top] = 0;
                    } else if (color[v] == 1) {
                        size_t from = top;
                        while (stack[from] != v) {
                            --from;
                        }
                        *kind = MATRIX_ADV_CYCLE_ZERO;
                        if (cycle != NULL) {
                            memcpy(cycle, stack + from, (top - from + 1) * sizeof(*cycle));
                        }
                        if (cycle_len != NULL) {
                            *cycle_len = top - from + 1;
                        }
                    }
                }
            }
        }
    }
    free(pot);
    free(par);
    free(stack);
    free(cursor);
    free(color);
    return status;
}

/* ================================================================== */
/*  Shared backtracking matcher (subgraph + graph isomorphism)          */
/* ================================================================== */

typedef struct {
    const MatrixGraph *p;       /* pattern */
    const MatrixGraph *t;       /* target  */
    bool         induced;
    bool         weights;
    const int   *pcol;          /* optional vertex colours; must match */
    const int   *tcol;
    const size_t *pin;          /* in-degrees (directed graphs) */
    const size_t *tin;
    bool         directed;
} MaMatch;

/** In-degree of every ID (self-loops included). */
static bool ma_indegrees(const MatrixGraph *g, size_t *in)
{
    const size_t n = matrix_graph_id_bound(g);
    memset(in, 0, (n > 0 ? n : 1) * sizeof(*in));
    for (size_t ui = 0; ui < n; ++ui) {
        if (!matrix_graph_is_alive(g, (int)ui)) {
            continue;
        }
        size_t cur = 0;
        int    v;
        while (matrix_graph_next_neighbor(g, (int)ui, &cur, &v, NULL)) {
            ++in[v];
        }
    }
    return true;
}

/** Can pattern vertex pv map to target vertex tv, given ord[0..depth) mapped? */
static bool ma_feasible(const MaMatch *m, int pv, int tv, const int *ord, size_t depth,
                        const int *map)
{
    if (m->pcol != NULL && m->pcol[pv] != m->tcol[tv]) {
        return false;
    }
    if (matrix_graph_degree(m->p, pv) > matrix_graph_degree(m->t, tv)) {
        return false;
    }
    if (m->directed && m->pin[pv] > m->tin[tv]) {
        return false;
    }
    double wp, wt;
    const bool ploop = matrix_graph_has_edge(m->p, pv, pv, &wp);
    const bool tloop = matrix_graph_has_edge(m->t, tv, tv, &wt);
    if (ploop ? (!tloop || (m->weights && wp != wt)) : (m->induced && tloop)) {
        return false;
    }
    for (size_t i = 0; i < depth; ++i) {
        const int q  = ord[i];
        const int tq = map[q];
        for (int dir = 0; dir < (m->directed ? 2 : 1); ++dir) {
            const int a  = dir == 0 ? pv : q;
            const int b  = dir == 0 ? q : pv;
            const int ta = dir == 0 ? tv : tq;
            const int tb = dir == 0 ? tq : tv;
            const bool pe = matrix_graph_has_edge(m->p, a, b, &wp);
            const bool te = matrix_graph_has_edge(m->t, ta, tb, &wt);
            if (pe ? (!te || (m->weights && wp != wt)) : (m->induced && te)) {
                return false;
            }
        }
    }
    return true;
}

/**
 * ma_match - backtracking search for an injective map pattern -> target that
 * preserves edges (and non-edges if `induced`, weights if `weights`), with
 * optional colour classes.  Pattern vertices are matched in a connectivity-
 * first order so each new vertex is checked against already-placed neighbours.
 * Worst case O(Nt^Vp); iterative, so recursion depth is no issue.
 */
static MatrixStatus ma_match(const MatrixGraph *p, const MatrixGraph *t, bool induced,
                             bool weights, const int *pcol, const int *tcol, int *mapping,
                             bool *found)
{
    const size_t np = matrix_graph_id_bound(p);
    const size_t nt = matrix_graph_id_bound(t);
    const size_t vp = matrix_graph_vertex_count(p);
    const size_t vt = matrix_graph_vertex_count(t);
    *found = false;
    if (mapping != NULL) {
        mu_fill_int(mapping, np, -1);
    }
    if (vp == 0) {
        *found = true;
        return GRAPH_OK;
    }
    if (vp > vt) {
        return GRAPH_OK;
    }

    int    *map  = ma_malloc(np, sizeof(*map));
    int    *ord  = ma_malloc(vp, sizeof(*ord));
    size_t *conn = mu_calloc(np, sizeof(*conn));
    bool   *placed = mu_calloc(np, sizeof(*placed));
    bool   *used = mu_calloc(nt, sizeof(*used));
    size_t *cand = ma_malloc(vp + 1, sizeof(*cand));
    size_t *pin  = ma_malloc(np, sizeof(*pin));
    size_t *tin  = ma_malloc(nt, sizeof(*tin));
    MatrixStatus status = GRAPH_OK;

    if (map == NULL || ord == NULL || conn == NULL || placed == NULL || used == NULL ||
        cand == NULL || pin == NULL || tin == NULL) {
        status = GRAPH_ERR_OUT_OF_MEMORY;
    } else {
        mu_fill_int(map, np, -1);
        MaMatch m = {p, t, induced, weights, pcol, tcol, pin, tin,
                     matrix_graph_is_directed(p)};
        ma_indegrees(p, pin);
        ma_indegrees(t, tin);

        /* Order: most connections to placed vertices first, then highest degree. */
        for (size_t d = 0; d < vp; ++d) {
            int best = -1;
            for (size_t v = 0; v < np; ++v) {
                if (placed[v] || !matrix_graph_is_alive(p, (int)v)) {
                    continue;
                }
                if (best < 0 || conn[v] > conn[best] ||
                    (conn[v] == conn[best] &&
                     matrix_graph_degree(p, (int)v) > matrix_graph_degree(p, best))) {
                    best = (int)v;
                }
            }
            ord[d]       = best;
            placed[best] = true;
            for (size_t v = 0; v < np; ++v) {
                if (!placed[v] && matrix_graph_is_alive(p, (int)v) &&
                    (matrix_graph_has_edge(p, best, (int)v, NULL) ||
                     matrix_graph_has_edge(p, (int)v, best, NULL))) {
                    ++conn[v];
                }
            }
        }

        size_t d = 0;
        cand[0]  = 0;
        for (;;) {
            if (d == vp) {
                *found = true;
                break;
            }
            const int pv = ord[d];
            bool advanced = false;
            while (cand[d] < nt) {
                const int tv = (int)cand[d]++;
                if (!matrix_graph_is_alive(t, tv) || used[tv]) {
                    continue;
                }
                if (ma_feasible(&m, pv, tv, ord, d, map)) {
                    map[pv]  = tv;
                    used[tv] = true;
                    advanced = true;
                    break;
                }
            }
            if (advanced) {
                ++d;
                if (d < vp) {
                    cand[d] = 0;
                }
                continue;
            }
            if (d == 0) {
                break;
            }
            --d;
            used[map[ord[d]]] = false;
            map[ord[d]]       = -1;
        }
        if (*found && mapping != NULL) {
            memcpy(mapping, map, np * sizeof(*map));
        }
    }
    free(map);
    free(ord);
    free(conn);
    free(placed);
    free(used);
    free(cand);
    free(pin);
    free(tin);
    return status;
}

/* ================================================================== */
/*  7. Subgraph isomorphism                                             */
/* ================================================================== */

/**
 * matrix_adv_subgraph_isomorphism - does `target` contain a copy of `pattern`?
 * Backtracking in the style of Ullmann / VF2: pattern vertices are placed in a
 * connectivity-first order and a candidate must agree with every already
 * placed vertex (edges both ways for directed graphs, self-loops, and — if
 * `induced` — non-edges too).  Candidates are pruned by out-/in-degree.
 *
 * induced == false finds a subgraph monomorphism; true an induced subgraph.
 * match_weights additionally requires equal edge weights (exact comparison).
 * mapping (optional, pattern id_bound entries) receives pattern -> target IDs
 * (-1 for dead pattern IDs).  *found tells whether a match exists.
 * Both graphs must have the same directedness (else UNSUPPORTED).  An empty
 * pattern matches anything.
 * Time O(Nt^Vp) worst case; memory O(Np + Nt).
 */
MatrixStatus matrix_adv_subgraph_isomorphism(const MatrixGraph *pattern, const MatrixGraph *target,
                                             bool induced, bool match_weights, int *mapping,
                                             bool *found)
{
    if (pattern == NULL || target == NULL || found == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    *found = false;
    if (matrix_graph_is_directed(pattern) != matrix_graph_is_directed(target)) {
        return GRAPH_ERR_UNSUPPORTED;
    }
    return ma_match(pattern, target, induced, match_weights, NULL, NULL, mapping, found);
}

/* ================================================================== */
/*  8. Coffman–Graham scheduling                                        */
/* ================================================================== */

static int ma_cmp_pos(const void *ctx, int a, int b)
{
    const int *pos = ctx;
    return pos[a] - pos[b];
}

/** Is descending list a (length la) lexicographically smaller than b? */
static bool ma_lex_less(const int *a, size_t la, const int *b, size_t lb)
{
    const size_t m = la < lb ? la : lb;
    for (size_t i = 0; i < m; ++i) {
        if (a[i] != b[i]) {
            return a[i] < b[i];
        }
    }
    return la < lb;
}

/**
 * matrix_adv_coffman_graham - schedule the unit-time tasks of a DAG on
 * `processors` identical machines.  An arc u -> v means u must finish before
 * v starts.  Optimal for processors == 2, a (2 - 2/m)-approximation otherwise.
 *
 * Steps: (1) transitive reduction (via reachability bitsets); (2) label
 * vertices 1..V: repeatedly give the next label to the vertex all of whose
 * successors are labelled and whose decreasing list of successor labels is
 * lexicographically smallest; (3) at each time step run the (up to m) ready
 * vertices with the highest labels.
 *
 * label[v] is the Coffman–Graham label, slot[v] the 0-based time step (both -1
 * for dead IDs); *makespan is the number of time steps.
 * Errors: UNSUPPORTED if undirected, CYCLE if not a DAG, INVALID_ARGUMENT for
 * processors == 0, TOO_LARGE if the N x N/64 bitsets overflow.
 * Time O(V^2 * D + E * V / 64) with D the largest reduced out-degree;
 * memory O(N^2 / 64) words.
 */
MatrixStatus matrix_adv_coffman_graham(const MatrixGraph *g, size_t processors, int *label,
                                       int *slot, size_t *makespan)
{
    if (g == NULL || label == NULL || slot == NULL || makespan == NULL || processors == 0) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    if (!matrix_graph_is_directed(g)) {
        return GRAPH_ERR_UNSUPPORTED;
    }
    const size_t n    = matrix_graph_id_bound(g);
    const size_t live = matrix_graph_vertex_count(g);
    mu_fill_int(label, n, -1);
    mu_fill_int(slot, n, -1);
    *makespan = 0;
    if (live == 0) {
        return GRAPH_OK;
    }
    const size_t words = mu_words(n);
    const size_t cells = mu_mul(n, words);
    if (cells == 0) {
        return GRAPH_ERR_TOO_LARGE;
    }

    int      *topo  = ma_malloc(n, sizeof(*topo));
    int      *pos   = ma_malloc(n, sizeof(*pos));
    int      *succ  = ma_malloc(n, sizeof(*succ));
    int      *npred = mu_calloc(n, sizeof(*npred));
    int      *cur   = ma_malloc(n, sizeof(*cur));
    int      *best  = ma_malloc(n, sizeof(*best));
    int      *pick  = ma_malloc(n, sizeof(*pick));
    uint64_t *reach = mu_calloc(cells, sizeof(*reach));
    uint64_t *red   = mu_calloc(cells, sizeof(*red));
    uint64_t *covered = mu_calloc(words, sizeof(*covered));
    MatrixStatus status = GRAPH_OK;

    if (topo == NULL || pos == NULL || succ == NULL || npred == NULL || cur == NULL ||
        best == NULL || pick == NULL || reach == NULL || red == NULL || covered == NULL) {
        status = GRAPH_ERR_OUT_OF_MEMORY;
        goto done;
    }
    if (!matrix_topological_sort(g, topo)) {
        status = matrix_has_cycle(g) ? GRAPH_ERR_CYCLE : GRAPH_ERR_OUT_OF_MEMORY;
        goto done;
    }
    for (size_t i = 0; i < live; ++i) {
        pos[topo[i]] = (int)i;
    }

    /* (1a) reach[v] = strict descendants of v, filled in reverse topological order. */
    for (size_t i = live; i-- > 0;) {
        const int v = topo[i];
        uint64_t *rv = reach + (size_t)v * words;
        size_t c = 0;
        int    s;
        while (matrix_graph_next_neighbor(g, v, &c, &s, NULL)) {
            mu_bit_set(rv, (size_t)s);
            const uint64_t *rs = reach + (size_t)s * words;
            for (size_t w = 0; w < words; ++w) {
                rv[w] |= rs[w];
            }
        }
    }
    /* (1b) keep arc u -> s only if no earlier (in topological order) successor
     *      of u already reaches s. */
    for (size_t ui = 0; ui < n; ++ui) {
        const int u = (int)ui;
        if (!matrix_graph_is_alive(g, u)) {
            continue;
        }
        size_t k = 0, c = 0;
        int    s;
        while (matrix_graph_next_neighbor(g, u, &c, &s, NULL)) {
            succ[k++] = s;
        }
        if (!ma_sort_idx(succ, k, ma_cmp_pos, pos)) {
            status = GRAPH_ERR_OUT_OF_MEMORY;
            goto done;
        }
        memset(covered, 0, words * sizeof(*covered));
        for (size_t i = 0; i < k; ++i) {
            if (mu_bit_test(covered, (size_t)succ[i])) {
                continue;
            }
            mu_bit_set(red + ui * words, (size_t)succ[i]);
            ++npred[succ[i]];
            const uint64_t *rs = reach + (size_t)succ[i] * words;
            for (size_t w = 0; w < words; ++w) {
                covered[w] |= rs[w];
            }
        }
    }

    /* (2) labelling. */
    for (int next_label = 1; (size_t)next_label <= live; ++next_label) {
        int    best_v = -1;
        size_t best_len = 0;
        for (size_t v = 0; v < n; ++v) {
            if (!matrix_graph_is_alive(g, (int)v) || label[v] != -1) {
                continue;
            }
            size_t len = 0;
            bool   ready = true;
            const uint64_t *rv = red + v * words;
            for (size_t w = 0; w < words && ready; ++w) {
                uint64_t x = rv[w];
                while (x != 0) {
                    const size_t s = (w << 6) + (size_t)ma_ctz64(x);
                    x &= x - 1;
                    if (label[s] == -1) {
                        ready = false;
                        break;
                    }
                    cur[len++] = label[s];
                }
            }
            if (!ready) {
                continue;
            }
            qsort(cur, len, sizeof(*cur), ma_cmp_int_desc);
            if (best_v < 0 || ma_lex_less(cur, len, best, best_len)) {
                memcpy(best, cur, len * sizeof(*cur));
                best_len = len;
                best_v   = (int)v;
            }
        }
        if (best_v < 0) {
            status = GRAPH_ERR_CYCLE;       /* cannot happen for a DAG */
            goto done;
        }
        label[best_v] = next_label;
    }

    /* (3) schedule by decreasing label. */
    size_t scheduled = 0, step = 0;
    while (scheduled < live) {
        size_t picked = 0;
        while (picked < processors) {
            int best_v = -1;
            for (size_t v = 0; v < n; ++v) {
                if (matrix_graph_is_alive(g, (int)v) && slot[v] == -1 && npred[v] == 0 &&
                    (best_v < 0 || label[v] > label[best_v])) {
                    best_v = (int)v;
                }
            }
            if (best_v < 0) {
                break;
            }
            slot[best_v]  = (int)step;
            pick[picked++] = best_v;
        }
        if (picked == 0) {
            status = GRAPH_ERR_CYCLE;
            goto done;
        }
        for (size_t i = 0; i < picked; ++i) {      /* successors become ready next step */
            const uint64_t *rv = red + (size_t)pick[i] * words;
            for (size_t w = 0; w < words; ++w) {
                uint64_t x = rv[w];
                while (x != 0) {
                    --npred[(w << 6) + (size_t)ma_ctz64(x)];
                    x &= x - 1;
                }
            }
        }
        scheduled += picked;
        ++step;
    }
    *makespan = step;

done:
    free(topo);
    free(pos);
    free(succ);
    free(npred);
    free(cur);
    free(best);
    free(pick);
    free(reach);
    free(red);
    free(covered);
    return status;
}

/* ================================================================== */
/*  9. Crossing number (heuristic)                                      */
/* ================================================================== */

/** Do chords (a,b) and (c,d) of a circle cross?  Arguments are positions. */
static bool ma_chords_cross(size_t a, size_t b, size_t c, size_t d)
{
    if (a == c || a == d || b == c || b == d) {
        return false;
    }
    const size_t lo1 = a < b ? a : b, hi1 = a < b ? b : a;
    const size_t lo2 = c < d ? c : d, hi2 = c < d ? d : c;
    return (lo1 < lo2 && lo2 < hi1 && hi1 < hi2) || (lo2 < lo1 && lo1 < hi2 && hi2 < hi1);
}

/**
 * matrix_adv_crossing_number_heuristic - an UPPER BOUND on the crossing number
 * of the underlying simple undirected graph (directions, self-loops and
 * parallel arcs ignored).  Computing the exact crossing number is NP-hard.
 *
 * Heuristic: draw the graph with all vertices on a circle and edges as
 * straight chords, so the crossings of a drawing are the interleaving chord
 * pairs.  Start from a BFS order (adjacent vertices end up close), then run
 * passes of adjacent-position swaps, accepting a swap only if it strictly
 * lowers the crossing count (the change is computed locally from the
 * edges at the two swapped vertices).  Each accepted swap lowers an integer,
 * so the search terminates; max_passes bounds it (0 = until no swap helps).
 * Because this is an actual drawing, the result is a valid upper bound, and 0
 * proves the graph is outerplanar-drawable but a positive value does NOT prove
 * non-planarity.
 *
 * *crossings receives the count; order (optional, vertex_count entries) the
 * circular order of the alive vertices.
 * Time O(E^2) for the initial count plus O(passes * V * D^2); memory O(V + E).
 */
MatrixStatus matrix_adv_crossing_number_heuristic(const MatrixGraph *g, size_t max_passes,
                                                  long long *crossings, int *order)
{
    if (g == NULL || crossings == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    *crossings = 0;
    const size_t n  = matrix_graph_id_bound(g);
    const size_t vc = matrix_graph_vertex_count(g);
    if (vc == 0) {
        return GRAPH_OK;
    }
    MuSym sym;
    if (!mu_sym_build_ex(g, &sym, false)) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    int    *ord   = ma_malloc(vc, sizeof(*ord));
    size_t *pos   = ma_malloc(n, sizeof(*pos));
    bool   *seen  = mu_calloc(n, sizeof(*seen));
    int    *eu    = ma_malloc(sym.off[n] / 2, sizeof(*eu));
    int    *ev    = ma_malloc(sym.off[n] / 2, sizeof(*ev));
    MatrixStatus status = GRAPH_OK;
    if (ord == NULL || pos == NULL || seen == NULL || eu == NULL || ev == NULL) {
        status = GRAPH_ERR_OUT_OF_MEMORY;
        goto done;
    }

    /* BFS order over every component, using the symmetrised adjacency. */
    size_t placed = 0;
    for (size_t s = 0; s < n; ++s) {
        if (!matrix_graph_is_alive(g, (int)s) || seen[s]) {
            continue;
        }
        size_t head = placed;
        seen[s] = true;
        ord[placed++] = (int)s;
        while (head < placed) {
            const int u = ord[head++];
            for (size_t k = sym.off[u]; k < sym.off[u + 1]; ++k) {
                const int v = sym.to[k];
                if (!seen[v]) {
                    seen[v] = true;
                    ord[placed++] = v;
                }
            }
        }
    }
    for (size_t i = 0; i < vc; ++i) {
        pos[ord[i]] = i;
    }
    size_t m = 0;
    for (size_t u = 0; u < n; ++u) {
        for (size_t k = sym.off[u]; k < sym.off[u + 1]; ++k) {
            if (sym.to[k] > (int)u) {
                eu[m] = (int)u;
                ev[m] = sym.to[k];
                ++m;
            }
        }
    }
    long long total = 0;
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = i + 1; j < m; ++j) {
            if (ma_chords_cross(pos[eu[i]], pos[ev[i]], pos[eu[j]], pos[ev[j]])) {
                ++total;
            }
        }
    }

    bool improved = vc >= 4;
    for (size_t pass = 0; improved && (max_passes == 0 || pass < max_passes); ++pass) {
        improved = false;
        for (size_t i = 0; i < vc; ++i) {
            const int u = ord[i];
            const int w = ord[(i + 1) % vc];
            long long before = 0, after = 0;
            for (size_t a = sym.off[u]; a < sym.off[u + 1]; ++a) {
                const int x = sym.to[a];
                if (x == w) {
                    continue;
                }
                for (size_t b = sym.off[w]; b < sym.off[w + 1]; ++b) {
                    const int y = sym.to[b];
                    if (y == u || y == x) {
                        continue;
                    }
                    before += ma_chords_cross(pos[u], pos[x], pos[w], pos[y]);
                    after  += ma_chords_cross(pos[w], pos[x], pos[u], pos[y]);
                }
            }
            if (after < before) {
                const size_t t = pos[u];
                pos[u] = pos[w];
                pos[w] = t;
                ord[pos[u]] = u;
                ord[pos[w]] = w;
                total -= before - after;
                improved = true;
            }
        }
    }
    *crossings = total;
    if (order != NULL) {
        memcpy(order, ord, vc * sizeof(*ord));
    }
done:
    free(ord);
    free(pos);
    free(seen);
    free(eu);
    free(ev);
    mu_sym_free(&sym);
    return status;
}

/* ================================================================== */
/*  10. Colour refinement (1-dimensional Weisfeiler–Leman)              */
/* ================================================================== */

typedef struct {
    const int    *sig;
    const size_t *off;
} MaSigCtx;

static int ma_cmp_sig(const void *ctx, int a, int b)
{
    const MaSigCtx *c = ctx;
    const size_t la = c->off[a + 1] - c->off[a];
    const size_t lb = c->off[b + 1] - c->off[b];
    const int *pa = c->sig + c->off[a];
    const int *pb = c->sig + c->off[b];
    const size_t m = la < lb ? la : lb;
    for (size_t i = 0; i < m; ++i) {
        if (pa[i] != pb[i]) {
            return pa[i] < pb[i] ? -1 : 1;
        }
    }
    return (la > lb) - (la < lb);
}

/**
 * ma_wl - colour refinement of graph a, or of the disjoint union of a and b
 * (so the colour names are comparable between them).  A vertex's new colour
 * is the rank of its signature (old colour, number of out/in neighbours,
 * sorted multiset of out-neighbour colours, sorted multiset of in-neighbour
 * colours) among all signatures, so the naming depends only on the isomorphism
 * type.  Edge weights are ignored.  O(rounds * (N + E) log N).
 */
static MatrixStatus ma_wl(const MatrixGraph *a, const MatrixGraph *b, int *ca, int *cb,
                          size_t *ncolors, size_t *rounds)
{
    const size_t na = matrix_graph_id_bound(a);
    const size_t nb = b != NULL ? matrix_graph_id_bound(b) : 0;
    const size_t N  = na + nb;
    *ncolors = 0;
    if (rounds != NULL) {
        *rounds = 0;
    }
    if (N == 0) {
        return GRAPH_OK;
    }
    size_t arcs = 0, alive_n = 0;
    for (size_t u = 0; u < na; ++u) {
        if (matrix_graph_is_alive(a, (int)u)) {
            arcs += matrix_graph_degree(a, (int)u);
            ++alive_n;
        }
    }
    for (size_t u = 0; u < nb; ++u) {
        if (matrix_graph_is_alive(b, (int)u)) {
            arcs += matrix_graph_degree(b, (int)u);
            ++alive_n;
        }
    }
    if (alive_n == 0) {
        for (size_t i = 0; i < na; ++i) {
            ca[i] = -1;
        }
        for (size_t i = 0; i < nb && cb != NULL; ++i) {
            cb[i] = -1;
        }
        return GRAPH_OK;
    }

    int    *src   = ma_malloc(arcs, sizeof(*src));
    int    *dst   = ma_malloc(arcs, sizeof(*dst));
    int    *col   = ma_malloc(N, sizeof(*col));
    int    *ncol  = ma_malloc(N, sizeof(*ncol));
    int    *idx   = ma_malloc(N, sizeof(*idx));
    size_t *ooff  = mu_calloc(N + 2, sizeof(*ooff));
    size_t *ioff  = mu_calloc(N + 2, sizeof(*ioff));
    size_t *sigoff = mu_calloc(N + 2, sizeof(*sigoff));
    int    *oto   = ma_malloc(arcs, sizeof(*oto));
    int    *ito   = ma_malloc(arcs, sizeof(*ito));
    size_t *ofill = ma_malloc(N, sizeof(*ofill));
    size_t *ifill = ma_malloc(N, sizeof(*ifill));
    int    *sig   = NULL;
    MatrixStatus status = GRAPH_OK;

    if (src == NULL || dst == NULL || col == NULL || ncol == NULL || idx == NULL ||
        ooff == NULL || ioff == NULL || sigoff == NULL || oto == NULL || ito == NULL ||
        ofill == NULL || ifill == NULL) {
        status = GRAPH_ERR_OUT_OF_MEMORY;
        goto done;
    }
    size_t k = 0;
    for (int part = 0; part < 2; ++part) {
        const MatrixGraph *gp = part == 0 ? a : b;
        const size_t base = part == 0 ? 0 : na;
        const size_t np   = part == 0 ? na : nb;
        for (size_t u = 0; gp != NULL && u < np; ++u) {
            col[base + u] = -1;
            if (!matrix_graph_is_alive(gp, (int)u)) {
                continue;
            }
            col[base + u] = 0;
            size_t cur = 0;
            int    v;
            while (matrix_graph_next_neighbor(gp, (int)u, &cur, &v, NULL)) {
                src[k] = (int)(base + u);
                dst[k] = (int)(base + (size_t)v);
                ++ooff[base + u + 1];
                ++ioff[base + (size_t)v + 1];
                ++k;
            }
        }
    }
    for (size_t i = 0; i < N; ++i) {
        ooff[i + 1] += ooff[i];
        ioff[i + 1] += ioff[i];
        ofill[i] = ooff[i];
        ifill[i] = ioff[i];
    }
    for (size_t e = 0; e < arcs; ++e) {
        oto[ofill[src[e]]++] = dst[e];
        ito[ifill[dst[e]]++] = src[e];
    }
    for (size_t v = 0; v < N; ++v) {
        sigoff[v + 1] = sigoff[v] +
            (col[v] < 0 ? 0 : 3 + (ooff[v + 1] - ooff[v]) + (ioff[v + 1] - ioff[v]));
    }
    sig = ma_malloc(sigoff[N], sizeof(*sig));
    if (sig == NULL) {
        status = GRAPH_ERR_OUT_OF_MEMORY;
        goto done;
    }

    size_t count_old = 1, round = 0;
    for (;;) {
        size_t alive_idx = 0;
        for (size_t v = 0; v < N; ++v) {
            if (col[v] < 0) {
                continue;
            }
            const size_t od = ooff[v + 1] - ooff[v];
            const size_t id = ioff[v + 1] - ioff[v];
            int *p = sig + sigoff[v];
            p[0] = col[v];
            p[1] = (int)od;
            p[2] = (int)id;
            for (size_t j = 0; j < od; ++j) {
                p[3 + j] = col[oto[ooff[v] + j]];
            }
            for (size_t j = 0; j < id; ++j) {
                p[3 + od + j] = col[ito[ioff[v] + j]];
            }
            qsort(p + 3, od, sizeof(*p), ma_cmp_int_asc);
            qsort(p + 3 + od, id, sizeof(*p), ma_cmp_int_asc);
            idx[alive_idx++] = (int)v;
        }
        MaSigCtx ctx = {sig, sigoff};
        if (!ma_sort_idx(idx, alive_idx, ma_cmp_sig, &ctx)) {
            status = GRAPH_ERR_OUT_OF_MEMORY;
            goto done;
        }
        size_t count = 0;
        for (size_t i = 0; i < alive_idx; ++i) {
            if (i > 0 && ma_cmp_sig(&ctx, idx[i - 1], idx[i]) != 0) {
                ++count;
            }
            ncol[idx[i]] = (int)count;
        }
        ++count;
        for (size_t v = 0; v < N; ++v) {
            if (col[v] >= 0) {
                col[v] = ncol[v];
            }
        }
        if (count == count_old) {
            count_old = count;
            break;
        }
        count_old = count;
        ++round;
    }
    *ncolors = count_old;
    if (rounds != NULL) {
        *rounds = round;
    }
    for (size_t i = 0; i < na; ++i) {
        ca[i] = col[i];
    }
    for (size_t i = 0; cb != NULL && i < nb; ++i) {
        cb[i] = col[na + i];
    }
done:
    free(src);
    free(dst);
    free(col);
    free(ncol);
    free(idx);
    free(ooff);
    free(ioff);
    free(sigoff);
    free(oto);
    free(ito);
    free(ofill);
    free(ifill);
    free(sig);
    return status;
}

/**
 * matrix_adv_color_refinement - 1-WL colour refinement to the stable colouring.
 * All vertices start with one colour; each round splits colours by the
 * multiset of neighbour colours (out- and in-neighbours for directed graphs).
 * color[v] is in 0..*num_colors-1 (-1 for dead IDs) and numbered canonically
 * (by sorted signature), so isomorphic graphs get identical colour histograms.
 * *rounds (optional) is the number of rounds that split a class.
 * Time O(rounds * (N + E) log N); rounds <= V.
 */
MatrixStatus matrix_adv_color_refinement(const MatrixGraph *g, int *color, size_t *num_colors,
                                         size_t *rounds)
{
    if (g == NULL || color == NULL || num_colors == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    return ma_wl(g, NULL, color, NULL, num_colors, rounds);
}

/* ================================================================== */
/*  11. Color-coding (k-vertex simple path)                             */
/* ================================================================== */

static uint64_t ma_rng_next(uint64_t *s)
{
    uint64_t x = *s;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    *s = x;
    return x * UINT64_C(2685821657736338717);
}

/**
 * matrix_adv_color_coding_path - Alon–Yuster–Zwick colour coding for a simple
 * path on exactly k vertices (k-1 edges).  Each trial colours the vertices
 * uniformly with k colours; a "colourful" path (all colours distinct) is then
 * found exactly by DP over colour subsets: reach[S] = vertices at which a
 * colourful path using exactly the colours S can end.  A fixed k-path is
 * colourful with probability k!/k^k >= e^-k, so the failure probability after
 * T trials is at most (1 - e^-k)^T.
 *
 * trials == 0 picks ceil(e^k * ln 100) (99 % success for an existing path),
 * capped at 100000.  seed makes the run reproducible.  One-sided error: a
 * returned path is always genuine; *found == false only means none was seen.
 * path (optional, k entries) receives the vertices in order.  Works for
 * directed (arcs followed forward) and undirected graphs.  k must be in
 * 1..20; the DP table needs 2^k * N / 8 bytes (TOO_LARGE above 1 GiB).
 * Time O(trials * 2^k * E); memory O(2^k * N / 64) words.
 */
MatrixStatus matrix_adv_color_coding_path(const MatrixGraph *g, size_t k, size_t trials,
                                          uint64_t seed, bool *found, int *path)
{
    if (g == NULL || found == NULL || k == 0 || k > MA_CC_MAX_K) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    *found = false;
    const size_t n  = matrix_graph_id_bound(g);
    const size_t vc = matrix_graph_vertex_count(g);
    if (vc < k) {
        return GRAPH_OK;
    }
    if (trials == 0) {
        const double t = ceil(exp((double)k) * log(100.0));
        trials = t > 100000.0 ? 100000 : (size_t)t;
    }
    const size_t words = mu_words(n);
    const size_t masks = (size_t)1 << k;
    const size_t table = mu_mul(masks, words);
    if (table == 0 || table > ((size_t)1 << 27)) {
        return GRAPH_ERR_TOO_LARGE;
    }
    uint64_t *dp    = ma_malloc(table, sizeof(*dp));
    int      *color = ma_malloc(n, sizeof(*color));
    if (dp == NULL || color == NULL) {
        free(dp);
        free(color);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    uint64_t rng = seed != 0 ? seed : UINT64_C(0x9E3779B97F4A7C15);
    const size_t full = masks - 1;

    for (size_t trial = 0; trial < trials && !*found; ++trial) {
        memset(dp, 0, table * sizeof(*dp));
        for (size_t v = 0; v < n; ++v) {
            color[v] = -1;
            if (matrix_graph_is_alive(g, (int)v)) {
                color[v] = (int)(ma_rng_next(&rng) % k);
                mu_bit_set(dp + ((size_t)1 << color[v]) * words, v);
            }
        }
        for (size_t mask = 1; mask < full; ++mask) {
            const uint64_t *row = dp + mask * words;
            for (size_t w = 0; w < words; ++w) {
                uint64_t x = row[w];
                while (x != 0) {
                    const int u = (int)((w << 6) + (size_t)ma_ctz64(x));
                    x &= x - 1;
                    size_t cur = 0;
                    int    v;
                    while (matrix_graph_next_neighbor(g, u, &cur, &v, NULL)) {
                        const size_t bit = (size_t)1 << color[v];
                        if ((mask & bit) == 0) {
                            mu_bit_set(dp + (mask | bit) * words, (size_t)v);
                        }
                    }
                }
            }
        }
        const uint64_t *last = dp + full * words;
        int end = -1;
        for (size_t w = 0; w < words && end < 0; ++w) {
            if (last[w] != 0) {
                end = (int)((w << 6) + (size_t)ma_ctz64(last[w]));
            }
        }
        if (k == 1) {                        /* full == 1: the base case is the answer */
            end = -1;
            for (size_t v = 0; v < n && end < 0; ++v) {
                if (color[v] >= 0) {
                    end = (int)v;
                }
            }
        }
        if (end < 0) {
            continue;
        }
        *found = true;
        if (path != NULL) {
            size_t mask = full;
            int    cur  = end;
            path[k - 1] = cur;
            for (size_t i = k - 1; i > 0; --i) {
                const size_t prev_mask = mask ^ ((size_t)1 << color[cur]);
                const uint64_t *prow = dp + prev_mask * words;
                int pred = -1;
                for (size_t u = 0; u < n && pred < 0; ++u) {
                    if (mu_bit_test(prow, u) && matrix_graph_has_edge(g, (int)u, cur, NULL)) {
                        pred = (int)u;
                    }
                }
                path[i - 1] = pred;
                cur  = pred;
                mask = prev_mask;
            }
        }
    }
    free(dp);
    free(color);
    return GRAPH_OK;
}

/* ================================================================== */
/*  12. Centrality measures                                             */
/* ================================================================== */

/**
 * matrix_adv_harmonic_centrality - H(v) = (1/(V-1)) * sum over u != v of
 * 1/d(v,u), with d the hop distance along out-arcs and 1/inf = 0.  Unlike
 * closeness it is well defined on disconnected graphs.  Result in [0,1];
 * dead IDs get 0.  Time O(V * (V + E)) via BFS.
 */
MatrixStatus matrix_adv_harmonic_centrality(const MatrixGraph *g, double *out)
{
    if (g == NULL || out == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n  = matrix_graph_id_bound(g);
    const size_t vc = matrix_graph_vertex_count(g);
    int *dist  = ma_malloc(n, sizeof(*dist));
    int *queue = ma_malloc(n, sizeof(*queue));
    if (dist == NULL || queue == NULL) {
        free(dist);
        free(queue);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    for (size_t s = 0; s < n; ++s) {
        out[s] = 0.0;
        if (!matrix_graph_is_alive(g, (int)s) || vc < 2) {
            continue;
        }
        const size_t reached = mu_bfs_hops(g, (int)s, dist, queue);
        double sum = 0.0;
        for (size_t i = 1; i < reached; ++i) {
            sum += 1.0 / (double)dist[queue[i]];
        }
        out[s] = sum / (double)(vc - 1);
    }
    free(dist);
    free(queue);
    return GRAPH_OK;
}

/**
 * matrix_adv_eigenvector_centrality - dominant eigenvector of the adjacency
 * matrix by power iteration on (A + I) (the shift keeps bipartite graphs from
 * oscillating), scores flowing along arcs so a vertex is central when central
 * vertices point to it.  Edge weights are ignored.  The result has unit
 * Euclidean length; dead IDs get 0.  On a disconnected graph the score
 * concentrates on the component with the largest eigenvalue.
 * max_iter == 0 uses 1000; tol == 0 uses 1e-9.  *iterations (optional) is
 * the number of rounds performed.  Time O(iterations * (V + E)).
 */
MatrixStatus matrix_adv_eigenvector_centrality(const MatrixGraph *g, size_t max_iter, double tol,
                                               double *out, size_t *iterations)
{
    if (g == NULL || out == NULL || !(tol >= 0.0)) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    if (max_iter == 0) {
        max_iter = 1000;
    }
    if (tol == 0.0) {
        tol = 1e-9;
    }
    const size_t n  = matrix_graph_id_bound(g);
    const size_t vc = matrix_graph_vertex_count(g);
    if (iterations != NULL) {
        *iterations = 0;
    }
    mu_fill_double(out, n, 0.0);
    if (vc == 0) {
        return GRAPH_OK;
    }
    double *next = ma_malloc(n, sizeof(*next));
    if (next == NULL) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    for (size_t i = 0; i < n; ++i) {
        out[i] = matrix_graph_is_alive(g, (int)i) ? 1.0 / sqrt((double)vc) : 0.0;
    }
    size_t it = 0;
    while (it < max_iter) {
        for (size_t i = 0; i < n; ++i) {
            next[i] = out[i];
        }
        for (size_t ui = 0; ui < n; ++ui) {
            if (!matrix_graph_is_alive(g, (int)ui)) {
                continue;
            }
            size_t cur = 0;
            int    v;
            while (matrix_graph_next_neighbor(g, (int)ui, &cur, &v, NULL)) {
                next[v] += out[ui];
            }
        }
        double norm = 0.0;
        for (size_t i = 0; i < n; ++i) {
            norm += next[i] * next[i];
        }
        norm = sqrt(norm);
        double diff = 0.0;
        for (size_t i = 0; i < n; ++i) {
            next[i] /= norm;
            diff = fmax(diff, fabs(next[i] - out[i]));
            out[i] = next[i];
        }
        ++it;
        if (diff < tol) {
            break;
        }
    }
    if (iterations != NULL) {
        *iterations = it;
    }
    free(next);
    return GRAPH_OK;
}

/**
 * matrix_adv_centrality - one entry point for the measures.  DEGREE, CLOSENESS
 * and BETWEENNESS delegate to matrix_degree_centrality,
 * matrix_closeness_centrality and matrix_betweenness_centrality (already in
 * matrix_analysis.c); HARMONIC and EIGENVECTOR are the new ones above
 * (EIGENVECTOR with default iteration limits).  out needs id_bound entries.
 */
MatrixStatus matrix_adv_centrality(const MatrixGraph *g, MatrixAdvCentralityKind kind, double *out)
{
    switch (kind) {
    case MATRIX_ADV_CENTRALITY_DEGREE:
        return matrix_degree_centrality(g, out);
    case MATRIX_ADV_CENTRALITY_CLOSENESS:
        return matrix_closeness_centrality(g, out);
    case MATRIX_ADV_CENTRALITY_BETWEENNESS:
        return matrix_betweenness_centrality(g, out);
    case MATRIX_ADV_CENTRALITY_HARMONIC:
        return matrix_adv_harmonic_centrality(g, out);
    case MATRIX_ADV_CENTRALITY_EIGENVECTOR:
        return matrix_adv_eigenvector_centrality(g, 0, 0.0, out, NULL);
    }
    return GRAPH_ERR_INVALID_ARGUMENT;
}

/* ================================================================== */
/*  13. A* search                                                       */
/* ================================================================== */

/** Heuristic value for v; false if the callback returns a negative or NaN value. */
static bool ma_heuristic(MatrixAdvHeuristic h, void *ctx, int v, int goal, double *out)
{
    const double x = h != NULL ? h(v, goal, ctx) : 0.0;
    if (isnan(x) || x < 0.0) {
        return false;
    }
    *out = x;
    return true;
}

/**
 * matrix_adv_astar - A* shortest path from source to target.  Entries are
 * ordered by f = g + h in a lazy-deletion binary heap, and a vertex is
 * re-expanded whenever a cheaper route to it is found, so the result is
 * optimal for any admissible heuristic (consistency is not required).
 * h == NULL degrades to Dijkstra.  The callback must return values >= 0
 * (otherwise INVALID_ARGUMENT) and is evaluated at most once per vertex.
 *
 * path receives the route (free with matrix_adv_path_free); *expanded
 * (optional) the number of expansions.  Weights must be >= 0
 * (NEGATIVE_WEIGHT otherwise); NOT_FOUND if target is unreachable.
 * Time O(E log V) with a consistent heuristic; may re-expand vertices (up to
 * exponentially often in contrived cases) with an inconsistent one.
 */
MatrixStatus matrix_adv_astar(const MatrixGraph *g, int source, int target,
                              MatrixAdvHeuristic h, void *ctx, MatrixAdvPath *path,
                              size_t *expanded)
{
    if (g == NULL || path == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    memset(path, 0, sizeof(*path));
    if (expanded != NULL) {
        *expanded = 0;
    }
    if (!matrix_graph_is_alive(g, source) || !matrix_graph_is_alive(g, target)) {
        return GRAPH_ERR_INVALID_VERTEX;
    }
    const size_t n = matrix_graph_id_bound(g);
    double *gs = ma_malloc(n, sizeof(*gs));
    double *hv = ma_malloc(n, sizeof(*hv));
    int    *par = ma_malloc(n, sizeof(*par));
    unsigned char *known = mu_calloc(n, sizeof(*known));
    MuHeap heap;
    const bool heap_ok = mu_heap_init(&heap, 64);
    MatrixStatus status = GRAPH_OK;
    size_t count = 0;

    if (gs == NULL || hv == NULL || par == NULL || known == NULL || !heap_ok) {
        status = GRAPH_ERR_OUT_OF_MEMORY;
        goto done;
    }
    mu_fill_double(gs, n, INFINITY);
    mu_fill_int(par, n, -1);
    gs[source] = 0.0;
    if (!ma_heuristic(h, ctx, source, target, &hv[source])) {
        status = GRAPH_ERR_INVALID_ARGUMENT;
        goto done;
    }
    known[source] = 1;
    if (!mu_heap_push(&heap, hv[source], source)) {
        status = GRAPH_ERR_OUT_OF_MEMORY;
        goto done;
    }
    status = GRAPH_ERR_NOT_FOUND;
    while (heap.size > 0) {
        double f;
        int    u;
        mu_heap_pop(&heap, &f, &u);
        if (f > gs[u] + hv[u]) {
            continue;                    /* stale entry: a cheaper route was found later */
        }
        if (u == target) {
            status = ma_path_from_parents(par, n, source, target, gs[target], path);
            break;
        }
        ++count;
        size_t cur = 0;
        int    v;
        double w;
        while (matrix_graph_next_neighbor(g, u, &cur, &v, &w)) {
            if (w < 0.0) {
                status = GRAPH_ERR_NEGATIVE_WEIGHT;
                goto done;
            }
            const double cand = gs[u] + w;
            if (cand < gs[v]) {
                if (!known[v]) {
                    if (!ma_heuristic(h, ctx, v, target, &hv[v])) {
                        status = GRAPH_ERR_INVALID_ARGUMENT;
                        goto done;
                    }
                    known[v] = 1;
                }
                gs[v]  = cand;
                par[v] = u;
                if (!mu_heap_push(&heap, cand + hv[v], v)) {
                    status = GRAPH_ERR_OUT_OF_MEMORY;
                    goto done;
                }
            }
        }
    }
done:
    if (expanded != NULL) {
        *expanded = count;
    }
    if (heap_ok) {
        mu_heap_free(&heap);
    }
    free(gs);
    free(hv);
    free(par);
    free(known);
    return status;
}

/* ================================================================== */
/*  14. Flooding                                                        */
/* ================================================================== */

/**
 * matrix_adv_flood - synchronous network flooding from `source`.  In round r
 * every vertex that first received the message in round r-1 transmits it over
 * every out-arc except the one back to the vertex it came from; a vertex that
 * has already seen the message ignores later copies (duplicate suppression),
 * but those copies still count as transmissions.
 *
 * round[v] is the round in which v was reached (0 for the source, -1 if never
 * reached or dead); parent[v] the first sender (-1 for the source / unreached;
 * parent may be NULL).  ttl is the hop limit (SIZE_MAX = unlimited; ttl == 0
 * reaches only the source).  *messages (optional) counts all transmissions,
 * *reached (optional) the informed vertices including the source.  For ttl
 * unlimited the reach set equals BFS reachability and round[] equals the
 * BFS level; the extra output is the message count.  Self-loops are not
 * transmitted.  Time O(V + E).
 */
MatrixStatus matrix_adv_flood(const MatrixGraph *g, int source, size_t ttl, int *round,
                              int *parent, size_t *messages, size_t *reached)
{
    if (g == NULL || round == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    if (!matrix_graph_is_alive(g, source)) {
        return GRAPH_ERR_INVALID_VERTEX;
    }
    const size_t n = matrix_graph_id_bound(g);
    int *frontier = ma_malloc(n, sizeof(*frontier));
    int *next     = ma_malloc(n, sizeof(*next));
    int *par      = parent != NULL ? parent : ma_malloc(n, sizeof(*par));
    if (frontier == NULL || next == NULL || par == NULL) {
        free(frontier);
        free(next);
        if (parent == NULL) {
            free(par);
        }
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    mu_fill_int(round, n, -1);
    mu_fill_int(par, n, -1);
    size_t sent = 0, informed = 1, fn = 1;
    frontier[0]  = source;
    round[source] = 0;
    for (size_t r = 1; fn > 0 && r <= ttl; ++r) {
        size_t nn = 0;
        for (size_t i = 0; i < fn; ++i) {
            const int u = frontier[i];
            size_t cur = 0;
            int    v;
            while (matrix_graph_next_neighbor(g, u, &cur, &v, NULL)) {
                if (v == u || v == par[u]) {
                    continue;
                }
                ++sent;
                if (round[v] < 0) {
                    round[v] = (int)r;
                    par[v]   = u;
                    next[nn++] = v;
                    ++informed;
                }
            }
        }
        int *t = frontier;
        frontier = next;
        next = t;
        fn = nn;
    }
    if (messages != NULL) {
        *messages = sent;
    }
    if (reached != NULL) {
        *reached = informed;
    }
    free(frontier);
    free(next);
    if (parent == NULL) {
        free(par);
    }
    return GRAPH_OK;
}

/* ================================================================== */
/*  15. Graph drawing (force-directed)                                  */
/* ================================================================== */

/**
 * matrix_adv_force_layout - Fruchterman–Reingold force-directed layout in the
 * unit square.  Every pair of vertices repels with force k^2/d, every edge
 * attracts with d^2/k (k = sqrt(1/V)); a temperature that cools linearly from
 * 0.1 to 0 caps the step length, and positions are clamped to [0,1]^2.  The
 * start is a deterministic circle, so the result is reproducible.  Direction,
 * weights and self-loops are ignored (the symmetrised graph is drawn);
 * components repel each other and settle in separate regions.
 *
 * x[v], y[v] receive the coordinates; dead IDs get (0,0); a single vertex is
 * placed at the centre.  iterations == 0 returns the initial circle.
 * Time O(iterations * (V^2 + E)); memory O(N + E).
 */
MatrixStatus matrix_adv_force_layout(const MatrixGraph *g, size_t iterations, double *x, double *y)
{
    if (g == NULL || x == NULL || y == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n  = matrix_graph_id_bound(g);
    const size_t vc = matrix_graph_vertex_count(g);
    mu_fill_double(x, n, 0.0);
    mu_fill_double(y, n, 0.0);
    if (vc == 0) {
        return GRAPH_OK;
    }
    if (vc == 1) {
        for (size_t v = 0; v < n; ++v) {
            if (matrix_graph_is_alive(g, (int)v)) {
                x[v] = y[v] = 0.5;
            }
        }
        return GRAPH_OK;
    }
    MuSym sym;
    if (!mu_sym_build_ex(g, &sym, false)) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    double *dx = ma_malloc(n, sizeof(*dx));
    double *dy = ma_malloc(n, sizeof(*dy));
    int    *ids = ma_malloc(vc, sizeof(*ids));
    if (dx == NULL || dy == NULL || ids == NULL) {
        free(dx);
        free(dy);
        free(ids);
        mu_sym_free(&sym);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    size_t m = 0;
    for (size_t v = 0; v < n; ++v) {
        if (matrix_graph_is_alive(g, (int)v)) {
            const double a = 2.0 * MA_PI * (double)m / (double)vc;
            x[v] = 0.5 + 0.4 * cos(a);
            y[v] = 0.5 + 0.4 * sin(a);
            ids[m++] = (int)v;
        }
    }
    const double k = sqrt(1.0 / (double)vc);
    for (size_t it = 0; it < iterations; ++it) {
        const double temp = 0.1 * (1.0 - (double)it / (double)iterations);
        for (size_t i = 0; i < vc; ++i) {
            dx[ids[i]] = dy[ids[i]] = 0.0;
        }
        for (size_t i = 0; i < vc; ++i) {
            for (size_t j = i + 1; j < vc; ++j) {
                const int a = ids[i], b = ids[j];
                double ddx = x[a] - x[b];
                double ddy = y[a] - y[b];
                double d   = sqrt(ddx * ddx + ddy * ddy);
                if (d < 1e-9) {              /* coincident: deterministic nudge */
                    ddx = ((i + j) & 1) ? 1e-6 : -1e-6;
                    ddy = 1e-6;
                    d   = sqrt(ddx * ddx + ddy * ddy);
                }
                const double f = k * k / d / d;
                dx[a] += ddx * f;
                dy[a] += ddy * f;
                dx[b] -= ddx * f;
                dy[b] -= ddy * f;
            }
        }
        for (size_t i = 0; i < vc; ++i) {
            const int a = ids[i];
            for (size_t e = sym.off[a]; e < sym.off[a + 1]; ++e) {
                const int b = sym.to[e];
                if (b < a) {
                    continue;
                }
                const double ddx = x[a] - x[b];
                const double ddy = y[a] - y[b];
                const double d   = sqrt(ddx * ddx + ddy * ddy);
                if (d < 1e-12) {
                    continue;
                }
                const double f = d / k;          /* (d^2/k) / d, applied to the vector */
                dx[a] -= ddx * f;
                dy[a] -= ddy * f;
                dx[b] += ddx * f;
                dy[b] += ddy * f;
            }
        }
        for (size_t i = 0; i < vc; ++i) {
            const int v = ids[i];
            const double len = sqrt(dx[v] * dx[v] + dy[v] * dy[v]);
            if (len > 0.0) {
                const double step = len < temp ? len : temp;
                x[v] += dx[v] / len * step;
                y[v] += dy[v] / len * step;
            }
            x[v] = fmin(1.0, fmax(0.0, x[v]));
            y[v] = fmin(1.0, fmax(0.0, y[v]));
        }
    }
    free(dx);
    free(dy);
    free(ids);
    mu_sym_free(&sym);
    return GRAPH_OK;
}

/* ================================================================== */
/*  16. Floyd–Warshall                                                  */
/* ================================================================== */

/** n * n, or 0 on overflow (n > 0). */
static size_t ma_square(size_t n)
{
    return mu_mul(n, n);
}

/**
 * matrix_adv_floyd_warshall - all-pairs shortest paths.  dist is N x N,
 * row-major, dist[u*N+v] = shortest distance (INFINITY if unreachable, and for
 * every row/column of a dead ID; 0 on the diagonal of alive vertices).  pred
 * (optional, N x N) holds the vertex before v on the path from u (-1 if none);
 * use matrix_adv_apsp_path to rebuild paths.
 *
 * Negative weights are allowed.  If a negative cycle exists the distances are
 * meaningless and GRAPH_ERR_NEGATIVE_CYCLE is returned.  An
 * undirected graph with a negative edge always has one (u-v-u).  Self-loops
 * with negative weight are negative cycles; positive ones are ignored.
 * Errors: TOO_LARGE if N*N overflows size_t.  Time O(N^3), memory O(N^2).
 */
MatrixStatus matrix_adv_floyd_warshall(const MatrixGraph *g, double *dist, int *pred)
{
    if (g == NULL || dist == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n = matrix_graph_id_bound(g);
    if (n == 0) {
        return GRAPH_OK;
    }
    if (ma_square(n) == 0) {
        return GRAPH_ERR_TOO_LARGE;
    }
    mu_fill_double(dist, n * n, INFINITY);
    if (pred != NULL) {
        mu_fill_int(pred, n * n, -1);
    }
    for (size_t ui = 0; ui < n; ++ui) {
        const int u = (int)ui;
        if (!matrix_graph_is_alive(g, u)) {
            continue;
        }
        dist[ui * n + ui] = 0.0;
        size_t cur = 0;
        int    v;
        double w;
        while (matrix_graph_next_neighbor(g, u, &cur, &v, &w)) {
            if (v == u) {
                if (w < 0.0) {
                    dist[ui * n + ui] = w;
                }
                continue;
            }
            dist[ui * n + (size_t)v] = w;
            if (pred != NULL) {
                pred[ui * n + (size_t)v] = u;
            }
        }
    }
    for (size_t k = 0; k < n; ++k) {
        for (size_t i = 0; i < n; ++i) {
            const double dik = dist[i * n + k];
            if (!isfinite(dik)) {
                continue;
            }
            for (size_t j = 0; j < n; ++j) {
                const double dkj = dist[k * n + j];
                if (isfinite(dkj) && dik + dkj < dist[i * n + j]) {
                    dist[i * n + j] = dik + dkj;
                    if (pred != NULL) {
                        pred[i * n + j] = pred[k * n + j];
                    }
                }
            }
        }
    }
    for (size_t i = 0; i < n; ++i) {
        if (dist[i * n + i] < 0.0) {
            return MA_ERR_NEG_CYCLE;
        }
    }
    return GRAPH_OK;
}

/* ================================================================== */
/*  17. Graph neural network (GCN forward pass)                         */
/* ================================================================== */

static double ma_activate(MatrixAdvActivation a, double x)
{
    switch (a) {
    case MATRIX_ADV_ACT_RELU:
        return x > 0.0 ? x : 0.0;
    case MATRIX_ADV_ACT_TANH:
        return tanh(x);
    case MATRIX_ADV_ACT_NONE:
        break;
    }
    return x;
}

/**
 * matrix_adv_gnn_forward - inference through a stack of graph-convolution
 * layers (Kipf & Welling GCN), i.e. basic message passing:
 *
 *     H' = act( A_hat * H * W + b ),   A_hat = D^-1/2 (A + I) D^-1/2
 *
 * A is the symmetrised, unweighted adjacency (directions and weights are
 * ignored, self-loops dropped, then I added), D its degree matrix.  Each
 * vertex therefore averages its own and its neighbours' feature rows before
 * the dense transform.  There is no training: the caller supplies weights.
 *
 * features is N x in_dim (row-major; rows of dead IDs are ignored); output is
 * N x layers[num_layers-1].out_dim, with zero rows for dead IDs.  Layer i's
 * in_dim must equal the previous out_dim (or the features' in_dim).
 * Errors: INVALID_ARGUMENT for NULL pointers, zero sizes or mismatched
 * dimensions; TOO_LARGE on size overflow.
 * Time O(L * (E * d + V * d_in * d_out)); memory O(N * d_max).
 */
MatrixStatus matrix_adv_gnn_forward(const MatrixGraph *g, const double *features, size_t in_dim,
                                    const MatrixAdvGnnLayer *layers, size_t num_layers,
                                    double *output)
{
    if (g == NULL || features == NULL || layers == NULL || output == NULL || in_dim == 0 ||
        num_layers == 0) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    size_t dim = in_dim, maxdim = in_dim;
    for (size_t l = 0; l < num_layers; ++l) {
        if (layers[l].in_dim != dim || layers[l].out_dim == 0 || layers[l].weights == NULL) {
            return GRAPH_ERR_INVALID_ARGUMENT;
        }
        dim = layers[l].out_dim;
        if (dim > maxdim) {
            maxdim = dim;
        }
    }
    const size_t n = matrix_graph_id_bound(g);
    if (n == 0) {
        return GRAPH_OK;
    }
    const size_t cells = mu_mul(n, maxdim);
    if (cells == 0) {
        return GRAPH_ERR_TOO_LARGE;
    }
    MuSym sym;
    if (!mu_sym_build_ex(g, &sym, false)) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    double *cur = ma_malloc(cells, sizeof(*cur));
    double *agg = ma_malloc(cells, sizeof(*agg));
    double *nxt = ma_malloc(cells, sizeof(*nxt));
    double *inv = ma_malloc(n, sizeof(*inv));
    MatrixStatus status = GRAPH_OK;
    if (cur == NULL || agg == NULL || nxt == NULL || inv == NULL) {
        status = GRAPH_ERR_OUT_OF_MEMORY;
        goto done;
    }
    for (size_t v = 0; v < n; ++v) {
        inv[v] = 1.0 / sqrt(1.0 + (double)mu_sym_degree(&sym, (int)v));
        for (size_t c = 0; c < in_dim; ++c) {
            cur[v * in_dim + c] = matrix_graph_is_alive(g, (int)v) ? features[v * in_dim + c] : 0.0;
        }
    }
    dim = in_dim;
    for (size_t l = 0; l < num_layers; ++l) {
        const MatrixAdvGnnLayer *L = &layers[l];
        for (size_t v = 0; v < n; ++v) {
            double *row = agg + v * dim;
            if (!matrix_graph_is_alive(g, (int)v)) {
                for (size_t c = 0; c < dim; ++c) {
                    row[c] = 0.0;
                }
                continue;
            }
            const double self = inv[v] * inv[v];
            for (size_t c = 0; c < dim; ++c) {
                row[c] = self * cur[v * dim + c];
            }
            for (size_t e = sym.off[v]; e < sym.off[v + 1]; ++e) {
                const size_t u = (size_t)sym.to[e];
                const double wgt = inv[v] * inv[u];
                for (size_t c = 0; c < dim; ++c) {
                    row[c] += wgt * cur[u * dim + c];
                }
            }
        }
        for (size_t v = 0; v < n; ++v) {
            for (size_t o = 0; o < L->out_dim; ++o) {
                double acc = 0.0;
                if (matrix_graph_is_alive(g, (int)v)) {
                    acc = L->bias != NULL ? L->bias[o] : 0.0;
                    for (size_t c = 0; c < dim; ++c) {
                        acc += agg[v * dim + c] * L->weights[c * L->out_dim + o];
                    }
                    acc = ma_activate(L->activation, acc);
                }
                nxt[v * L->out_dim + o] = acc;
            }
        }
        double *t = cur;
        cur = nxt;
        nxt = t;
        dim = L->out_dim;
    }
    memcpy(output, cur, n * dim * sizeof(*output));
done:
    free(cur);
    free(agg);
    free(nxt);
    free(inv);
    mu_sym_free(&sym);
    return status;
}

/* ================================================================== */
/*  18. Graph isomorphism                                               */
/* ================================================================== */

/**
 * matrix_adv_graph_isomorphism - are a and b isomorphic?  (Same directedness,
 * same structure of arcs and self-loops, and — if match_weights — equal
 * weights.)  Cheap necessary conditions first (vertex and edge counts), then
 * colour refinement on the disjoint union (identical colour histograms are
 * required), then exact backtracking restricted to vertices of equal colour
 * (ma_match with induced = true).  Refinement alone decides most non-
 * isomorphic pairs and prunes the search drastically on the rest.
 *
 * *isomorphic is the answer.  mapping (optional, id_bound(a) entries) receives
 * an isomorphism a -> b when one exists (-1 for dead IDs).  Graphs of
 * different directedness are simply not isomorphic.
 * Time O((N + E) log N * rounds) for refinement; exponential worst case for
 * the search (e.g. strongly regular graphs); memory O(N + E).
 */
MatrixStatus matrix_adv_graph_isomorphism(const MatrixGraph *a, const MatrixGraph *b,
                                          bool match_weights, bool *isomorphic, int *mapping)
{
    if (a == NULL || b == NULL || isomorphic == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    *isomorphic = false;
    if (mapping != NULL) {
        mu_fill_int(mapping, matrix_graph_id_bound(a), -1);
    }
    if (matrix_graph_is_directed(a) != matrix_graph_is_directed(b) ||
        matrix_graph_vertex_count(a) != matrix_graph_vertex_count(b) ||
        matrix_graph_edge_count(a) != matrix_graph_edge_count(b)) {
        return GRAPH_OK;
    }
    const size_t na = matrix_graph_id_bound(a);
    const size_t nb = matrix_graph_id_bound(b);
    int *ca = ma_malloc(na, sizeof(*ca));
    int *cb = ma_malloc(nb, sizeof(*cb));
    size_t *ha = NULL, *hb = NULL;
    size_t ncolors = 0;
    MatrixStatus status = GRAPH_OK;
    if (ca == NULL || cb == NULL) {
        status = GRAPH_ERR_OUT_OF_MEMORY;
        goto done;
    }
    status = ma_wl(a, b, ca, cb, &ncolors, NULL);
    if (status != GRAPH_OK) {
        goto done;
    }
    ha = mu_calloc(ncolors, sizeof(*ha));
    hb = mu_calloc(ncolors, sizeof(*hb));
    if (ha == NULL || hb == NULL) {
        status = GRAPH_ERR_OUT_OF_MEMORY;
        goto done;
    }
    for (size_t i = 0; i < na; ++i) {
        if (ca[i] >= 0) {
            ++ha[ca[i]];
        }
    }
    for (size_t i = 0; i < nb; ++i) {
        if (cb[i] >= 0) {
            ++hb[cb[i]];
        }
    }
    if (ncolors > 0 && memcmp(ha, hb, ncolors * sizeof(*ha)) != 0) {
        goto done;                         /* refinement already separates them */
    }
    status = ma_match(a, b, true, match_weights, ca, cb, mapping, isomorphic);
done:
    free(ca);
    free(cb);
    free(ha);
    free(hb);
    return status;
}

/* ================================================================== */
/*  19. Iterative-deepening A*                                          */
/* ================================================================== */

/**
 * matrix_adv_ida_star - IDA*: depth-first searches bounded by f = g + h,
 * the bound rising each round to the smallest f that exceeded it.  Uses
 * O(depth) memory instead of A*'s open list.  Vertices already on the current
 * path are not revisited, so zero-weight cycles cannot loop; there is no
 * transposition table, so a vertex may be re-expanded many times.
 *
 * The result is optimal when h is admissible (never overestimates);
 * h == NULL gives iterative deepening on path cost.  The callback must return
 * values >= 0 (else INVALID_ARGUMENT).  Weights must be >= 0.
 * node_limit (0 = none) caps the number of generated nodes; exceeding it
 * returns TOO_LARGE, since the running time is exponential in the worst case.
 * NOT_FOUND if target is unreachable.
 * Time O(b^d) worst case (b branching factor, d solution depth); memory
 * O(V).  path is freed with matrix_adv_path_free.
 */
MatrixStatus matrix_adv_ida_star(const MatrixGraph *g, int source, int target,
                                 MatrixAdvHeuristic h, void *ctx, size_t node_limit,
                                 MatrixAdvPath *path)
{
    if (g == NULL || path == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    memset(path, 0, sizeof(*path));
    if (!matrix_graph_is_alive(g, source) || !matrix_graph_is_alive(g, target)) {
        return GRAPH_ERR_INVALID_VERTEX;
    }
    const size_t n = matrix_graph_id_bound(g);
    int      *stack  = ma_malloc(n + 1, sizeof(*stack));
    size_t   *cursor = ma_malloc(n + 1, sizeof(*cursor));
    double   *gs     = ma_malloc(n + 1, sizeof(*gs));
    uint64_t *onpath = mu_calloc(mu_words(n), sizeof(*onpath));
    MatrixStatus status = GRAPH_OK;
    size_t nodes = 0;
    if (stack == NULL || cursor == NULL || gs == NULL || onpath == NULL) {
        status = GRAPH_ERR_OUT_OF_MEMORY;
        goto done;
    }
    if (source == target) {
        int *vs = malloc(sizeof(*vs));
        if (vs == NULL) {
            status = GRAPH_ERR_OUT_OF_MEMORY;
        } else {
            vs[0] = source;
            path->vertices = vs;
            path->length   = 1;
            path->cost     = 0.0;
        }
        goto done;
    }
    double bound;
    if (!ma_heuristic(h, ctx, source, target, &bound)) {
        status = GRAPH_ERR_INVALID_ARGUMENT;
        goto done;
    }
    for (;;) {
        double next_bound = INFINITY;
        size_t depth = 0;
        stack[0]  = source;
        cursor[0] = 0;
        gs[0]     = 0.0;
        mu_bit_set(onpath, (size_t)source);
        for (;;) {
            const int u = stack[depth];
            bool descended = false;
            size_t *cur = &cursor[depth];
            int    v;
            double w;
            while (matrix_graph_next_neighbor(g, u, cur, &v, &w)) {
                if (w < 0.0) {
                    status = GRAPH_ERR_NEGATIVE_WEIGHT;
                    goto done;
                }
                if (mu_bit_test(onpath, (size_t)v)) {
                    continue;
                }
                double hval;
                if (!ma_heuristic(h, ctx, v, target, &hval)) {
                    status = GRAPH_ERR_INVALID_ARGUMENT;
                    goto done;
                }
                const double gv = gs[depth] + w;
                const double f  = gv + hval;
                if (node_limit != 0 && ++nodes > node_limit) {
                    status = GRAPH_ERR_TOO_LARGE;
                    goto done;
                }
                if (f > bound) {
                    if (f < next_bound) {
                        next_bound = f;
                    }
                    continue;
                }
                if (v == target) {
                    int *vs = malloc((depth + 2) * sizeof(*vs));
                    if (vs == NULL) {
                        status = GRAPH_ERR_OUT_OF_MEMORY;
                        goto done;
                    }
                    memcpy(vs, stack, (depth + 1) * sizeof(*vs));
                    vs[depth + 1]  = v;
                    path->vertices = vs;
                    path->length   = depth + 2;
                    path->cost     = gv;
                    goto done;
                }
                ++depth;
                stack[depth]  = v;
                cursor[depth] = 0;
                gs[depth]     = gv;
                mu_bit_set(onpath, (size_t)v);
                descended = true;
                break;
            }
            if (descended) {
                continue;
            }
            mu_bit_clear(onpath, (size_t)u);
            if (depth == 0) {
                break;
            }
            --depth;
        }
        if (next_bound == INFINITY) {
            status = GRAPH_ERR_NOT_FOUND;
            goto done;
        }
        bound = next_bound;
    }
done:
    free(stack);
    free(cursor);
    free(gs);
    free(onpath);
    return status;
}

/* ================================================================== */
/*  20. Iterative-deepening DFS                                         */
/* ================================================================== */

/**
 * matrix_adv_iddfs - iterative-deepening depth-first search for a path with
 * the fewest edges from source to target.  Depth-limited DFS is run with
 * limits 0, 1, 2, ...; within one round a vertex is expanded again only if it
 * is reached at a smaller depth than before (so each round is O(V + E) rather
 * than exponential).  The search ends when the target is found, when a round
 * is never cut off by the limit (target unreachable), or at max_depth
 * (SIZE_MAX = unbounded; the limit never needs to exceed V - 1).
 *
 * path->cost is the sum of the edge weights along the path found; weights do
 * not influence which path is chosen.  NOT_FOUND if no path within max_depth.
 * Time O(min(max_depth, V) * (V + E)) — the usual IDDFS re-visit overhead;
 * memory O(V).
 */
MatrixStatus matrix_adv_iddfs(const MatrixGraph *g, int source, int target, size_t max_depth,
                              MatrixAdvPath *path)
{
    if (g == NULL || path == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    memset(path, 0, sizeof(*path));
    if (!matrix_graph_is_alive(g, source) || !matrix_graph_is_alive(g, target)) {
        return GRAPH_ERR_INVALID_VERTEX;
    }
    const size_t n    = matrix_graph_id_bound(g);
    const size_t live = matrix_graph_vertex_count(g);
    const size_t top_limit = (max_depth < live - 1) ? max_depth : live - 1;
    int    *stack  = ma_malloc(n, sizeof(*stack));
    size_t *cursor = ma_malloc(n, sizeof(*cursor));
    size_t *bestd  = ma_malloc(n, sizeof(*bestd));
    double *gs     = ma_malloc(n, sizeof(*gs));
    MatrixStatus status = GRAPH_ERR_NOT_FOUND;
    if (stack == NULL || cursor == NULL || bestd == NULL || gs == NULL) {
        status = GRAPH_ERR_OUT_OF_MEMORY;
        goto done;
    }
    for (size_t limit = 0; limit <= top_limit && status == GRAPH_ERR_NOT_FOUND; ++limit) {
        for (size_t i = 0; i < n; ++i) {
            bestd[i] = SIZE_MAX;
        }
        bool   cutoff = false;
        size_t depth  = 0;
        stack[0]  = source;
        cursor[0] = 0;
        gs[0]     = 0.0;
        bestd[source] = 0;
        if (source == target) {
            int *vs = malloc(sizeof(*vs));
            if (vs == NULL) {
                status = GRAPH_ERR_OUT_OF_MEMORY;
            } else {
                vs[0] = source;
                path->vertices = vs;
                path->length   = 1;
                path->cost     = 0.0;
                status = GRAPH_OK;
            }
            break;
        }
        for (;;) {
            const int u = stack[depth];
            bool descended = false;
            if (depth == limit) {
                if (matrix_graph_degree(g, u) > 0) {
                    cutoff = true;
                }
            } else {
                int    v;
                double w;
                while (matrix_graph_next_neighbor(g, u, &cursor[depth], &v, &w)) {
                    if (bestd[v] <= depth + 1) {
                        continue;
                    }
                    bestd[v] = depth + 1;
                    if (v == target) {
                        int *vs = malloc((depth + 2) * sizeof(*vs));
                        if (vs == NULL) {
                            status = GRAPH_ERR_OUT_OF_MEMORY;
                        } else {
                            memcpy(vs, stack, (depth + 1) * sizeof(*vs));
                            vs[depth + 1]  = v;
                            path->vertices = vs;
                            path->length   = depth + 2;
                            path->cost     = gs[depth] + w;
                            status = GRAPH_OK;
                        }
                        goto round_done;
                    }
                    ++depth;
                    stack[depth]  = v;
                    cursor[depth] = 0;
                    gs[depth]     = gs[depth - 1] + w;
                    descended = true;
                    break;
                }
            }
            if (descended) {
                continue;
            }
            if (depth == 0) {
                break;
            }
            --depth;
        }
round_done:
        if (status == GRAPH_ERR_NOT_FOUND && !cutoff) {
            break;                      /* nothing was cut off: the target is unreachable */
        }
    }
done:
    free(stack);
    free(cursor);
    free(bestd);
    free(gs);
    return status;
}

/* ================================================================== */
/*  21. Jump point search                                               */
/* ================================================================== */

typedef struct {
    const MatrixGraph *g;
    size_t             w;       /* grid width; cell (x, y) is vertex y*w + x */
    size_t             n;
} MaGrid;

static bool ma_walk(const MaGrid *gr, long x, long y)
{
    if (x < 0 || y < 0 || (size_t)x >= gr->w || (size_t)y > gr->n / gr->w) {
        return false;
    }
    const size_t idx = (size_t)y * gr->w + (size_t)x;
    return idx < gr->n && matrix_graph_is_alive(gr->g, (int)idx);
}

/** Jump from (x,y) in direction (dx,dy); true with *jx,*jy if a jump point is found. */
static bool ma_jump(const MaGrid *gr, long x, long y, int dx, int dy, long gx, long gy,
                    long *jx, long *jy)
{
    for (;;) {
        x += dx;
        y += dy;
        if (!ma_walk(gr, x, y)) {
            return false;
        }
        if (x == gx && y == gy) {
            *jx = x;
            *jy = y;
            return true;
        }
        bool jump_point = false;
        if (dx != 0 && dy != 0) {
            jump_point = (!ma_walk(gr, x - dx, y) && ma_walk(gr, x - dx, y + dy)) ||
                         (!ma_walk(gr, x, y - dy) && ma_walk(gr, x + dx, y - dy));
            long tx, ty;
            if (!jump_point && (ma_jump(gr, x, y, dx, 0, gx, gy, &tx, &ty) ||
                                ma_jump(gr, x, y, 0, dy, gx, gy, &tx, &ty))) {
                jump_point = true;
            }
        } else if (dx != 0) {
            jump_point = (!ma_walk(gr, x, y + 1) && ma_walk(gr, x + dx, y + 1)) ||
                         (!ma_walk(gr, x, y - 1) && ma_walk(gr, x + dx, y - 1));
        } else {
            jump_point = (!ma_walk(gr, x + 1, y) && ma_walk(gr, x + 1, y + dy)) ||
                         (!ma_walk(gr, x - 1, y) && ma_walk(gr, x - 1, y + dy));
        }
        if (jump_point) {
            *jx = x;
            *jy = y;
            return true;
        }
    }
}

static double ma_octile(long dx, long dy)
{
    const double ax = (double)(dx < 0 ? -dx : dx);
    const double ay = (double)(dy < 0 ? -dy : dy);
    return ax > ay ? ax + (MA_SQRT2 - 1.0) * ay : ay + (MA_SQRT2 - 1.0) * ax;
}

static int ma_sign(long v)
{
    return (v > 0) - (v < 0);
}

/**
 * matrix_adv_jump_point_search - Jump Point Search (Harabor & Grastien) for
 * shortest paths on a uniform-cost 8-connected grid.  JPS is a grid-specific
 * algorithm, so the matrix graph is *interpreted* as a grid instead of adding
 * a new representation: vertex id = y * width + x, and a cell is walkable iff
 * its vertex is alive (remove a vertex to place an obstacle).  Edges and edge
 * weights are NOT used; moves cost 1 (orthogonal) or sqrt(2) (diagonal), and
 * a diagonal step is allowed whenever its target cell is walkable.  IDs at or
 * beyond id_bound are treated as blocked, so id_bound need not be a multiple
 * of width.
 *
 * A* runs over jump points only, with the octile-distance heuristic; the
 * returned path lists every cell, not only the jump points; path->cost is the
 * octile length.  Errors: INVALID_ARGUMENT for width == 0 or > id_bound,
 * INVALID_VERTEX for a blocked source/target, NOT_FOUND if unreachable.
 * Time O(N log N) worst case, usually far fewer expansions than plain A*;
 * memory O(N).
 */
MatrixStatus matrix_adv_jump_point_search(const MatrixGraph *g, size_t width, int source,
                                          int target, MatrixAdvPath *path)
{
    if (g == NULL || path == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    memset(path, 0, sizeof(*path));
    const size_t n = matrix_graph_id_bound(g);
    if (width == 0 || width > n) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    if (!matrix_graph_is_alive(g, source) || !matrix_graph_is_alive(g, target)) {
        return GRAPH_ERR_INVALID_VERTEX;
    }
    const MaGrid gr = {g, width, n};
    const long w  = (long)width;
    const long gx = target % w, gy = target / w;

    double        *gs = ma_malloc(n, sizeof(*gs));
    int           *par = ma_malloc(n, sizeof(*par));
    unsigned char *closed = mu_calloc(n, sizeof(*closed));
    MuHeap heap;
    const bool heap_ok = mu_heap_init(&heap, 64);
    MatrixStatus status = GRAPH_ERR_NOT_FOUND;
    if (gs == NULL || par == NULL || closed == NULL || !heap_ok) {
        status = GRAPH_ERR_OUT_OF_MEMORY;
        goto done;
    }
    mu_fill_double(gs, n, INFINITY);
    mu_fill_int(par, n, -1);
    gs[source] = 0.0;
    if (!mu_heap_push(&heap, ma_octile(gx - source % w, gy - source / w), source)) {
        status = GRAPH_ERR_OUT_OF_MEMORY;
        goto done;
    }
    while (heap.size > 0) {
        double f;
        int    u;
        mu_heap_pop(&heap, &f, &u);
        if (closed[u]) {
            continue;
        }
        closed[u] = 1;
        if (u == target) {
            size_t count = 1;
            for (int v = target; par[v] >= 0; v = par[v]) {
                const long ddx = v % w - par[v] % w, ddy = v / w - par[v] / w;
                count += (size_t)((ddx < 0 ? -ddx : ddx) > (ddy < 0 ? -ddy : ddy)
                                      ? (ddx < 0 ? -ddx : ddx) : (ddy < 0 ? -ddy : ddy));
            }
            int *vs = malloc(count * sizeof(*vs));
            if (vs == NULL) {
                status = GRAPH_ERR_OUT_OF_MEMORY;
                break;
            }
            size_t at = count;
            long cx = gx, cy = gy;
            vs[--at] = target;
            for (int v = target; par[v] >= 0; v = par[v]) {
                const long px = par[v] % w, py = par[v] / w;
                const int sx = ma_sign(px - cx), sy = ma_sign(py - cy);
                while (cx != px || cy != py) {
                    cx += sx;
                    cy += sy;
                    vs[--at] = (int)(cy * w + cx);
                }
            }
            path->vertices = vs;
            path->length   = count;
            path->cost     = gs[target];
            status = GRAPH_OK;
            break;
        }
        const long x = u % w, y = u / w;
        int dirs[8][2];
        size_t nd = 0;
        if (par[u] < 0) {
            for (int ddy = -1; ddy <= 1; ++ddy) {
                for (int ddx = -1; ddx <= 1; ++ddx) {
                    if (ddx != 0 || ddy != 0) {
                        dirs[nd][0] = ddx;
                        dirs[nd][1] = ddy;
                        ++nd;
                    }
                }
            }
        } else {
            const int dx = ma_sign(x - par[u] % w);
            const int dy = ma_sign(y - par[u] / w);
            if (dx != 0 && dy != 0) {
                const int c[3][2] = {{dx, dy}, {dx, 0}, {0, dy}};
                for (int i = 0; i < 3; ++i) {
                    dirs[nd][0] = c[i][0];
                    dirs[nd][1] = c[i][1];
                    ++nd;
                }
                if (!ma_walk(&gr, x - dx, y) && ma_walk(&gr, x - dx, y + dy)) {
                    dirs[nd][0] = -dx;
                    dirs[nd][1] = dy;
                    ++nd;
                }
                if (!ma_walk(&gr, x, y - dy) && ma_walk(&gr, x + dx, y - dy)) {
                    dirs[nd][0] = dx;
                    dirs[nd][1] = -dy;
                    ++nd;
                }
            } else if (dx != 0) {
                dirs[nd][0] = dx;
                dirs[nd][1] = 0;
                ++nd;
                if (!ma_walk(&gr, x, y + 1) && ma_walk(&gr, x + dx, y + 1)) {
                    dirs[nd][0] = dx;
                    dirs[nd][1] = 1;
                    ++nd;
                }
                if (!ma_walk(&gr, x, y - 1) && ma_walk(&gr, x + dx, y - 1)) {
                    dirs[nd][0] = dx;
                    dirs[nd][1] = -1;
                    ++nd;
                }
            } else {
                dirs[nd][0] = 0;
                dirs[nd][1] = dy;
                ++nd;
                if (!ma_walk(&gr, x + 1, y) && ma_walk(&gr, x + 1, y + dy)) {
                    dirs[nd][0] = 1;
                    dirs[nd][1] = dy;
                    ++nd;
                }
                if (!ma_walk(&gr, x - 1, y) && ma_walk(&gr, x - 1, y + dy)) {
                    dirs[nd][0] = -1;
                    dirs[nd][1] = dy;
                    ++nd;
                }
            }
        }
        for (size_t i = 0; i < nd; ++i) {
            long jx, jy;
            if (!ma_jump(&gr, x, y, dirs[i][0], dirs[i][1], gx, gy, &jx, &jy)) {
                continue;
            }
            const int j = (int)(jy * w + jx);
            const double cost = gs[u] + ma_octile(jx - x, jy - y);
            if (!closed[j] && cost < gs[j]) {
                gs[j]  = cost;
                par[j] = u;
                if (!mu_heap_push(&heap, cost + ma_octile(gx - jx, gy - jy), j)) {
                    status = GRAPH_ERR_OUT_OF_MEMORY;
                    goto done;
                }
            }
        }
    }
done:
    if (heap_ok) {
        mu_heap_free(&heap);
    }
    free(gs);
    free(par);
    free(closed);
    return status;
}

/* ================================================================== */
/*  22. Johnson's algorithm                                             */
/* ================================================================== */

/**
 * matrix_adv_johnson - all-pairs shortest paths for sparse graphs with
 * negative weights.  Bellman-Ford from a virtual source yields potentials h;
 * every weight is replaced by w + h(u) - h(v) >= 0 (clamped at 0 against
 * rounding), which lets the existing matrix_dijkstra run from every source;
 * distances are then un-shifted: d(u,v) = d'(u,v) - h(u) + h(v).  If no weight
 * is negative the original graph is used directly.
 *
 * Output layout is the same as matrix_adv_floyd_warshall: dist is N x N
 * (INFINITY for unreachable and dead IDs, 0 on the diagonal of alive
 * vertices), pred (optional) N x N with pred[s*N+v] the vertex before v on
 * the shortest path from s; use matrix_adv_apsp_path.  Negative cycles return
 * GRAPH_ERR_NEGATIVE_CYCLE; an undirected graph with a negative edge
 * has one.  Distances through reweighted edges can differ from exact
 * arithmetic by rounding error.
 * Errors: TOO_LARGE if N*N overflows.
 * Time O(V*E + V*(N*N/64 + E) log V) with matrix_dijkstra's row scans;
 * memory O(N^2) for the output.
 */
MatrixStatus matrix_adv_johnson(const MatrixGraph *g, double *dist, int *pred)
{
    if (g == NULL || dist == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n = matrix_graph_id_bound(g);
    if (n == 0) {
        return GRAPH_OK;
    }
    if (ma_square(n) == 0) {
        return GRAPH_ERR_TOO_LARGE;
    }
    mu_fill_double(dist, n * n, INFINITY);
    if (pred != NULL) {
        mu_fill_int(pred, n * n, -1);
    }
    double *pot = ma_malloc(n, sizeof(*pot));
    int    *par = ma_malloc(n, sizeof(*par));
    int    *scratch = pred == NULL ? ma_malloc(n, sizeof(*scratch)) : NULL;
    MatrixGraph *rew = NULL;
    MatrixStatus status = GRAPH_OK;
    if (pot == NULL || par == NULL || (pred == NULL && scratch == NULL)) {
        status = GRAPH_ERR_OUT_OF_MEMORY;
        goto done;
    }
    int witness;
    status = ma_potentials(g, 0.0, pot, par, &witness);
    if (status != GRAPH_OK) {
        goto done;
    }
    const MatrixGraph *work = g;
    if (ma_has_negative_weight(g)) {
        rew = ma_clone(g, pot);
        if (rew == NULL) {
            status = GRAPH_ERR_OUT_OF_MEMORY;
            goto done;
        }
        work = rew;
    } else {
        mu_fill_double(pot, n, 0.0);
    }
    for (size_t s = 0; s < n; ++s) {
        if (!matrix_graph_is_alive(g, (int)s)) {
            continue;
        }
        double *row = dist + s * n;
        if (!matrix_dijkstra(work, (int)s, row, pred != NULL ? pred + s * n : scratch)) {
            status = GRAPH_ERR_OUT_OF_MEMORY;     /* inputs were validated above */
            goto done;
        }
        for (size_t v = 0; v < n; ++v) {
            if (isfinite(row[v])) {
                row[v] += pot[v] - pot[s];
            }
        }
        row[s] = 0.0;
    }
done:
    matrix_graph_destroy(rew);
    free(pot);
    free(par);
    free(scratch);
    return status;
}

/* ================================================================== */
/*  23. K shortest paths (Yen)                                          */
/* ================================================================== */

/**
 * ma_restricted_dijkstra - shortest src -> dst path that avoids the vertices
 * in blocked_v and the arcs (be_u[i], be_v[i]) (both directions if the graph
 * is undirected).  Uses the caller's scratch arrays and heap.
 */
static MatrixStatus ma_restricted_dijkstra(const MatrixGraph *g, int src, int dst,
                                           const uint64_t *blocked_v, const int *be_u,
                                           const int *be_v, size_t n_be, double *dist, int *par,
                                           uint64_t *settled, MuHeap *heap, MatrixAdvPath *out)
{
    const size_t n = matrix_graph_id_bound(g);
    const bool   directed = matrix_graph_is_directed(g);
    mu_fill_double(dist, n, INFINITY);
    mu_fill_int(par, n, -1);
    memset(settled, 0, mu_words(n) * sizeof(*settled));
    heap->size = 0;
    dist[src]  = 0.0;
    if (!mu_heap_push(heap, 0.0, src)) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    while (heap->size > 0) {
        double d;
        int    u;
        mu_heap_pop(heap, &d, &u);
        if (mu_bit_test(settled, (size_t)u)) {
            continue;
        }
        mu_bit_set(settled, (size_t)u);
        if (u == dst) {
            break;
        }
        size_t cur = 0;
        int    v;
        double w;
        while (matrix_graph_next_neighbor(g, u, &cur, &v, &w)) {
            if (w < 0.0) {
                return GRAPH_ERR_NEGATIVE_WEIGHT;
            }
            if (mu_bit_test(blocked_v, (size_t)v)) {
                continue;
            }
            bool blocked = false;
            for (size_t i = 0; i < n_be && !blocked; ++i) {
                blocked = (be_u[i] == u && be_v[i] == v) ||
                          (!directed && be_u[i] == v && be_v[i] == u);
            }
            if (blocked) {
                continue;
            }
            if (d + w < dist[v]) {
                dist[v] = d + w;
                par[v]  = u;
                if (!mu_heap_push(heap, dist[v], v)) {
                    return GRAPH_ERR_OUT_OF_MEMORY;
                }
            }
        }
    }
    if (!isfinite(dist[dst])) {
        return GRAPH_ERR_NOT_FOUND;
    }
    return ma_path_from_parents(par, n, src, dst, dist[dst], out);
}

static bool ma_prefix_eq(const int *a, const int *b, size_t len)
{
    return memcmp(a, b, len * sizeof(*a)) == 0;
}

static bool ma_path_eq(const MatrixAdvPath *a, const MatrixAdvPath *b)
{
    return a->length == b->length && ma_prefix_eq(a->vertices, b->vertices, a->length);
}

/**
 * matrix_adv_k_shortest_paths - Yen's algorithm for the k shortest LOOPLESS
 * paths from source to target.  The first path is a Dijkstra path; each next
 * path is the cheapest "deviation" of an accepted one: for every vertex of the
 * last accepted path as spur node, the root prefix and the arcs already used
 * by accepted paths with the same prefix are blocked and Dijkstra completes
 * the route from the spur node.  Candidates are kept in a pool until chosen.
 *
 * paths must have room for k entries; on success *found (<= k) entries are
 * filled in non-decreasing cost order and each must be released with
 * matrix_adv_path_free.  Fewer than k are returned when fewer exist.
 * source == target yields the single path [source].  Weights must be >= 0
 * (NEGATIVE_WEIGHT).  NOT_FOUND (with *found == 0) if target is unreachable.
 * On any error nothing is left allocated and *found is 0.  Undirected edges
 * may be used in either direction.
 * Time O(k * V * (E log V)) (V spur computations per path); memory O(k * V).
 */
MatrixStatus matrix_adv_k_shortest_paths(const MatrixGraph *g, int source, int target, size_t k,
                                         MatrixAdvPath *paths, size_t *found)
{
    if (g == NULL || paths == NULL || found == NULL || k == 0) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    *found = 0;
    if (!matrix_graph_is_alive(g, source) || !matrix_graph_is_alive(g, target)) {
        return GRAPH_ERR_INVALID_VERTEX;
    }
    const size_t n     = matrix_graph_id_bound(g);
    const size_t words = mu_words(n);
    double   *dist    = ma_malloc(n, sizeof(*dist));
    int      *par     = ma_malloc(n, sizeof(*par));
    uint64_t *blocked = mu_calloc(words, sizeof(*blocked));
    uint64_t *settled = mu_calloc(words, sizeof(*settled));
    int      *be_u    = ma_malloc(k, sizeof(*be_u));
    int      *be_v    = ma_malloc(k, sizeof(*be_v));
    MuHeap heap;
    const bool heap_ok = mu_heap_init(&heap, 64);
    MatrixAdvPath *cand = NULL;
    size_t ncand = 0, capcand = 0;
    MatrixStatus status = GRAPH_OK;

    if (dist == NULL || par == NULL || blocked == NULL || settled == NULL || be_u == NULL ||
        be_v == NULL || !heap_ok) {
        status = GRAPH_ERR_OUT_OF_MEMORY;
        goto done;
    }
    memset(&paths[0], 0, sizeof(paths[0]));
    status = ma_restricted_dijkstra(g, source, target, blocked, be_u, be_v, 0, dist, par,
                                    settled, &heap, &paths[0]);
    if (status != GRAPH_OK) {
        goto done;
    }
    *found = 1;

    while (*found < k) {
        const MatrixAdvPath *prev = &paths[*found - 1];
        double rootcost = 0.0;
        for (size_t i = 0; i + 1 < prev->length; ++i) {
            size_t n_be = 0;
            for (size_t j = 0; j < *found; ++j) {
                const MatrixAdvPath *p = &paths[j];
                if (p->length > i + 1 && ma_prefix_eq(p->vertices, prev->vertices, i + 1)) {
                    be_u[n_be] = p->vertices[i];
                    be_v[n_be] = p->vertices[i + 1];
                    ++n_be;
                }
            }
            memset(blocked, 0, words * sizeof(*blocked));
            for (size_t j = 0; j < i; ++j) {
                mu_bit_set(blocked, (size_t)prev->vertices[j]);
            }
            MatrixAdvPath spur;
            memset(&spur, 0, sizeof(spur));
            const MatrixStatus st = ma_restricted_dijkstra(g, prev->vertices[i], target, blocked,
                                                           be_u, be_v, n_be, dist, par, settled,
                                                           &heap, &spur);
            if (st == GRAPH_OK) {
                MatrixAdvPath total;
                total.length   = i + spur.length;
                total.cost     = rootcost + spur.cost;
                total.vertices = malloc(total.length * sizeof(*total.vertices));
                if (total.vertices == NULL) {
                    matrix_adv_path_free(&spur);
                    status = GRAPH_ERR_OUT_OF_MEMORY;
                    goto done;
                }
                memcpy(total.vertices, prev->vertices, i * sizeof(int));
                memcpy(total.vertices + i, spur.vertices, spur.length * sizeof(int));
                matrix_adv_path_free(&spur);
                bool dup = false;
                for (size_t j = 0; j < ncand && !dup; ++j) {
                    dup = ma_path_eq(&cand[j], &total);
                }
                for (size_t j = 0; j < *found && !dup; ++j) {
                    dup = ma_path_eq(&paths[j], &total);
                }
                if (dup) {
                    matrix_adv_path_free(&total);
                } else {
                    MatrixAdvPath *c = ma_grow(cand, &capcand, ncand + 1, sizeof(*cand));
                    if (c == NULL) {
                        matrix_adv_path_free(&total);
                        status = GRAPH_ERR_OUT_OF_MEMORY;
                        goto done;
                    }
                    cand = c;
                    cand[ncand++] = total;
                }
            } else if (st != GRAPH_ERR_NOT_FOUND) {
                status = st;
                goto done;
            }
            double w;
            if (matrix_graph_has_edge(g, prev->vertices[i], prev->vertices[i + 1], &w)) {
                rootcost += w;
            }
        }
        if (ncand == 0) {
            break;
        }
        size_t bi = 0;
        for (size_t j = 1; j < ncand; ++j) {
            if (cand[j].cost < cand[bi].cost) {
                bi = j;
            }
        }
        paths[(*found)++] = cand[bi];
        cand[bi] = cand[--ncand];
    }
done:
    for (size_t j = 0; j < ncand; ++j) {
        matrix_adv_path_free(&cand[j]);
    }
    free(cand);
    if (status != GRAPH_OK) {
        for (size_t j = 0; j < *found; ++j) {
            matrix_adv_path_free(&paths[j]);
        }
        *found = 0;
    }
    if (heap_ok) {
        mu_heap_free(&heap);
    }
    free(dist);
    free(par);
    free(blocked);
    free(settled);
    free(be_u);
    free(be_v);
    return status;
}

/* ================================================================== */
/*  24. Widest path                                                     */
/* ================================================================== */

/**
 * matrix_adv_widest_path - maximum-bottleneck ("widest") paths from source:
 * the path whose smallest edge weight is as large as possible.  Dijkstra with
 * the roles swapped: the vertex with the largest bottleneck so far is settled
 * first, and an arc (u,v,w) offers min(width[u], w) to v.  Weights are read as
 * capacities and may be any finite value (including negative or 0).  Directed
 * and undirected graphs are handled; self-loops never help.
 *
 * width[v] is the best bottleneck (+INFINITY for the source, -INFINITY if v is
 * unreachable or dead); parent[v] the tree edge (-1 where undefined; may be
 * NULL).  Rebuild a path with matrix_adv_path_from_parents, passing width[v]
 * as the cost.  Time O(E log V).
 */
MatrixStatus matrix_adv_widest_path(const MatrixGraph *g, int source, double *width, int *parent)
{
    if (g == NULL || width == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    if (!matrix_graph_is_alive(g, source)) {
        return GRAPH_ERR_INVALID_VERTEX;
    }
    const size_t n = matrix_graph_id_bound(g);
    uint64_t *settled = mu_calloc(mu_words(n), sizeof(*settled));
    MuHeap heap;
    const bool heap_ok = mu_heap_init(&heap, 64);
    if (settled == NULL || !heap_ok) {
        free(settled);
        if (heap_ok) {
            mu_heap_free(&heap);
        }
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    mu_fill_double(width, n, -INFINITY);
    if (parent != NULL) {
        mu_fill_int(parent, n, -1);
    }
    MatrixStatus status = GRAPH_OK;
    width[source] = INFINITY;
    if (!mu_heap_push(&heap, -INFINITY, source)) {
        status = GRAPH_ERR_OUT_OF_MEMORY;
    }
    while (status == GRAPH_OK && heap.size > 0) {
        double key;
        int    u;
        mu_heap_pop(&heap, &key, &u);        /* key = -width: smallest key = widest */
        if (mu_bit_test(settled, (size_t)u)) {
            continue;
        }
        mu_bit_set(settled, (size_t)u);
        size_t cur = 0;
        int    v;
        double w;
        while (matrix_graph_next_neighbor(g, u, &cur, &v, &w)) {
            const double cand = width[u] < w ? width[u] : w;
            if (cand > width[v]) {
                width[v] = cand;
                if (parent != NULL) {
                    parent[v] = u;
                }
                if (!mu_heap_push(&heap, -cand, v)) {
                    status = GRAPH_ERR_OUT_OF_MEMORY;
                    break;
                }
            }
        }
    }
    mu_heap_free(&heap);
    free(settled);
    return status;
}

/* ================================================================== */
/*  25. Lexicographic BFS                                               */
/* ================================================================== */

/**
 * matrix_adv_lex_bfs - lexicographic breadth-first search (Rose, Tarjan,
 * Lueker) by partition refinement.  Vertices carry a class number; each step
 * visits an unvisited vertex of the highest class (smallest ID on ties) and
 * splits every class into "neighbours of the visited vertex" (placed ahead)
 * and the rest, which is equivalent to comparing the lexicographic labels of
 * the original formulation.  The reverse of the visit order is a perfect
 * elimination ordering exactly when the graph is chordal.
 *
 * Directions and self-loops are ignored (the underlying simple undirected
 * graph is searched).  start >= 0 forces that vertex to be visited first
 * (INVALID_VERTEX if dead); start < 0 lets the algorithm choose.  order
 * receives the vertex_count vertices in visit order; disconnected graphs are
 * covered component by component.
 * Time O(V^2 + E) (counting-sort renumbering in each of the V steps);
 * memory O(N + E).
 */
MatrixStatus matrix_adv_lex_bfs(const MatrixGraph *g, int start, int *order)
{
    if (g == NULL || order == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    if (start >= 0 && !matrix_graph_is_alive(g, start)) {
        return GRAPH_ERR_INVALID_VERTEX;
    }
    const size_t n  = matrix_graph_id_bound(g);
    const size_t vc = matrix_graph_vertex_count(g);
    if (vc == 0) {
        return GRAPH_OK;
    }
    MuSym sym;
    if (!mu_sym_build_ex(g, &sym, false)) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    int *cls   = ma_malloc(n, sizeof(*cls));
    int *rank  = mu_calloc(2 * vc + 2, sizeof(*rank));
    unsigned char *done = mu_calloc(n, sizeof(*done));
    unsigned char *mark = mu_calloc(n, sizeof(*mark));
    if (cls == NULL || rank == NULL || done == NULL || mark == NULL) {
        free(cls);
        free(rank);
        free(done);
        free(mark);
        mu_sym_free(&sym);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    for (size_t v = 0; v < n; ++v) {
        cls[v] = 0;
    }
    for (size_t step = 0; step < vc; ++step) {
        int pick = -1;
        if (step == 0 && start >= 0) {
            pick = start;
        } else {
            for (size_t v = 0; v < n; ++v) {
                if (!done[v] && matrix_graph_is_alive(g, (int)v) &&
                    (pick < 0 || cls[v] > cls[pick])) {
                    pick = (int)v;
                }
            }
        }
        done[pick]  = 1;
        order[step] = pick;
        for (size_t e = sym.off[pick]; e < sym.off[pick + 1]; ++e) {
            if (!done[sym.to[e]]) {
                mark[sym.to[e]] = 1;
            }
        }
        /* New key 2*class + (is neighbour); compress the keys in order. */
        memset(rank, 0, (2 * vc + 2) * sizeof(*rank));
        for (size_t v = 0; v < n; ++v) {
            if (!done[v] && matrix_graph_is_alive(g, (int)v)) {
                rank[2 * (size_t)cls[v] + mark[v]] = 1;
            }
        }
        int next = 0;
        for (size_t key = 0; key < 2 * vc + 2; ++key) {
            if (rank[key] != 0) {
                rank[key] = next++;
            }
        }
        for (size_t v = 0; v < n; ++v) {
            if (!done[v] && matrix_graph_is_alive(g, (int)v)) {
                cls[v] = rank[2 * (size_t)cls[v] + mark[v]];
            }
        }
        for (size_t e = sym.off[pick]; e < sym.off[pick + 1]; ++e) {
            mark[sym.to[e]] = 0;
        }
    }
    free(cls);
    free(rank);
    free(done);
    free(mark);
    mu_sym_free(&sym);
    return GRAPH_OK;
}

/* ================================================================== */
/*  26. Crossing count for a given vertex order                         */
/* ================================================================== */

/**
 * matrix_adv_crossing_count - number of edge crossings when the alive
 * vertices are placed on a circle in the given order and every edge is drawn
 * as a straight chord (the exact value, not a bound, for THIS drawing; the
 * crossing number itself is the minimum over all drawings and is NP-hard).
 * Directions, self-loops and parallel arcs are ignored.
 *
 * order must list every alive vertex exactly once (order_len == vertex_count).
 * Errors: INVALID_ARGUMENT for NULL or a wrong length, INVALID_VERTEX for a
 * dead/out-of-range/duplicate entry.  Time O(E^2); memory O(V + E).
 */
MatrixStatus matrix_adv_crossing_count(const MatrixGraph *g, const int *order,
                                       size_t order_len, long long *crossings)
{
    if (g == NULL || order == NULL || crossings == NULL ||
        order_len != matrix_graph_vertex_count(g)) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    *crossings = 0;
    const size_t n = matrix_graph_id_bound(g);
    size_t *pos = ma_malloc(n, sizeof(*pos));
    bool *seen = mu_calloc(n, sizeof(*seen));
    if (pos == NULL || seen == NULL) {
        free(pos);
        free(seen);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    MatrixStatus status = GRAPH_OK;
    for (size_t i = 0; i < order_len; ++i) {
        const int v = order[i];
        if (v < 0 || (size_t)v >= n || !matrix_graph_is_alive(g, v) || seen[v]) {
            status = GRAPH_ERR_INVALID_VERTEX;
            goto done;
        }
        seen[v] = true;
        pos[v] = i;
    }
    MuSym sym;
    if (!mu_sym_build_ex(g, &sym, false)) {
        status = GRAPH_ERR_OUT_OF_MEMORY;
        goto done;
    }
    const size_t m = sym.off[n] / 2;
    int *eu = ma_malloc(m, sizeof(*eu));
    int *ev = ma_malloc(m, sizeof(*ev));
    if (eu == NULL || ev == NULL) {
        status = GRAPH_ERR_OUT_OF_MEMORY;
    } else {
        size_t k = 0;
        for (size_t u = 0; u < n; ++u) {
            for (size_t e = sym.off[u]; e < sym.off[u + 1]; ++e) {
                if (sym.to[e] > (int)u) {
                    eu[k] = (int)u;
                    ev[k] = sym.to[e];
                    ++k;
                }
            }
        }
        long long total = 0;
        for (size_t i = 0; i < k; ++i) {
            for (size_t j = i + 1; j < k; ++j) {
                total += ma_chords_cross(pos[eu[i]], pos[ev[i]], pos[eu[j]], pos[ev[j]]);
            }
        }
        *crossings = total;
    }
    free(eu);
    free(ev);
    mu_sym_free(&sym);
done:
    free(pos);
    free(seen);
    return status;
}

/* ================================================================== */
/*  27. Circular layout                                                 */
/* ================================================================== */

/**
 * matrix_adv_circle_layout - place the alive vertices evenly on a circle of
 * the given radius centred at the origin, in increasing ID order (vertex i of
 * V at angle 2*pi*i/V).  Deterministic and independent of the edges.
 *
 * x[v], y[v] receive the coordinates; dead IDs get (0,0).  radius must be
 * finite and >= 0.  Time O(N); memory O(1).
 */
MatrixStatus matrix_adv_circle_layout(const MatrixGraph *g, double radius, double *x, double *y)
{
    if (g == NULL || x == NULL || y == NULL || !(radius >= 0.0) || isinf(radius)) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n  = matrix_graph_id_bound(g);
    const size_t vc = matrix_graph_vertex_count(g);
    mu_fill_double(x, n, 0.0);
    mu_fill_double(y, n, 0.0);
    size_t rank = 0;
    for (size_t v = 0; v < n; ++v) {
        if (matrix_graph_is_alive(g, (int)v)) {
            const double a = 2.0 * MA_PI * (double)rank / (double)vc;
            x[v] = radius * cos(a);
            y[v] = radius * sin(a);
            ++rank;
        }
    }
    return GRAPH_OK;
}

/* ================================================================== */
/*  28. Dynamic connectivity                                            */
/* ================================================================== */

struct MatrixAdvDynConn {
    const MatrixGraph *g;       /* borrowed; must outlive the tracker */
    MuUF               uf;
    size_t             n;       /* id_bound the union-find was built for */
    size_t             components;
    uint64_t           version; /* matrix_graph_version at the last rebuild */
};

/** Rebuilds the union-find from the graph's current edges.  False on OOM. */
static bool ma_dc_rebuild(MatrixAdvDynConn *dc)
{
    const size_t n = matrix_graph_id_bound(dc->g);
    mu_uf_free(&dc->uf);
    dc->uf.parent = dc->uf.rank = NULL;
    dc->n = 0;
    if (!mu_uf_init(&dc->uf, n)) {
        dc->uf.parent = dc->uf.rank = NULL;
        return false;
    }
    dc->n = n;
    size_t comps = matrix_graph_vertex_count(dc->g);
    for (size_t u = 0; u < n; ++u) {
        size_t cur = 0;
        int to;
        while (matrix_graph_next_neighbor(dc->g, (int)u, &cur, &to, NULL)) {
            if (mu_uf_union(&dc->uf, (int)u, to)) {
                --comps;
            }
        }
    }
    dc->components = comps;
    dc->version = matrix_graph_version(dc->g);
    return true;
}

/** Brings the structure up to date with the graph (lazy; O(V + E) after a change). */
static MatrixStatus ma_dc_sync(MatrixAdvDynConn *dc)
{
    if (dc->n != 0 || dc->uf.parent != NULL) {
        if (dc->version == matrix_graph_version(dc->g)) {
            return GRAPH_OK;
        }
    }
    return ma_dc_rebuild(dc) ? GRAPH_OK : GRAPH_ERR_OUT_OF_MEMORY;
}

/**
 * matrix_adv_dynconn_create - connectivity tracker for a graph that changes.
 * Connectivity is WEAK: edge directions are ignored.  The tracker borrows g
 * (it must outlive the tracker) and never modifies it.  After the graph
 * changes (edges or vertices added or removed), the next query notices the new
 * matrix_graph_version and rebuilds in O(V + E) before answering; queries on an
 * unchanged graph cost O(alpha(V)) each.  Returns NULL on NULL or allocation
 * failure.  Free with matrix_adv_dynconn_destroy.
 */
MatrixAdvDynConn *matrix_adv_dynconn_create(const MatrixGraph *g)
{
    if (g == NULL) {
        return NULL;
    }
    MatrixAdvDynConn *dc = calloc(1, sizeof(*dc));
    if (dc == NULL) {
        return NULL;
    }
    dc->g = g;
    if (!ma_dc_rebuild(dc)) {
        free(dc);
        return NULL;
    }
    return dc;
}

/** Frees the tracker (not the graph).  Safe with NULL. */
void matrix_adv_dynconn_destroy(MatrixAdvDynConn *dc)
{
    if (dc != NULL) {
        mu_uf_free(&dc->uf);
        free(dc);
    }
}

/**
 * matrix_adv_dynconn_connected - are u and v in the same weakly connected
 * component of the graph as it is NOW?  Errors: INVALID_ARGUMENT for NULL,
 * INVALID_VERTEX if u or v is dead, OUT_OF_MEMORY if a rebuild fails.
 */
MatrixStatus matrix_adv_dynconn_connected(MatrixAdvDynConn *dc, int u, int v, bool *connected)
{
    if (dc == NULL || connected == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const MatrixStatus st = ma_dc_sync(dc);
    if (st != GRAPH_OK) {
        return st;
    }
    if (!mu_ok(dc->g, u) || !mu_ok(dc->g, v)) {
        return GRAPH_ERR_INVALID_VERTEX;
    }
    *connected = mu_uf_find(&dc->uf, u) == mu_uf_find(&dc->uf, v);
    return GRAPH_OK;
}

/** Number of weakly connected components of the graph as it is now. */
MatrixStatus matrix_adv_dynconn_components(MatrixAdvDynConn *dc, size_t *count)
{
    if (dc == NULL || count == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const MatrixStatus st = ma_dc_sync(dc);
    if (st != GRAPH_OK) {
        return st;
    }
    *count = dc->components;
    return GRAPH_OK;
}
