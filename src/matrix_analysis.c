/*
 * matrix_analysis.c — shared result helpers, traversal callbacks, centrality,
 * clustering and other whole-graph measures.  See matrix_algorithms.h.
 */

#include "matrix_util.h"

/* ------------------------------------------------------------------ */
/*  Shared result helpers                                               */
/* ------------------------------------------------------------------ */

void matrix_groups_free(MatrixGroups *groups)
{
    if (groups == NULL) {
        return;
    }
    free(groups->offsets);
    free(groups->members);
    memset(groups, 0, sizeof(*groups));
}

void matrix_path_free(MatrixPath *path)
{
    if (path == NULL) {
        return;
    }
    free(path->vertices);
    memset(path, 0, sizeof(*path));
}

/* ------------------------------------------------------------------ */
/*  Traversal with callbacks                                            */
/* ------------------------------------------------------------------ */

size_t matrix_bfs_visit(const MatrixGraph *g, int start, MatrixVisitFn fn, void *ctx)
{
    if (!mu_ok(g, start) || fn == NULL) {
        return 0;
    }
    const size_t n = matrix_graph_id_bound(g);
    int *queue  = malloc(n * sizeof(int));
    int *parent = malloc(n * sizeof(int));
    uint64_t *seen = mu_calloc(mu_words(n), sizeof(uint64_t));
    if (queue == NULL || parent == NULL || seen == NULL) {
        free(queue);
        free(parent);
        free(seen);
        return 0;
    }
    size_t head = 0, tail = 0, visited = 0;
    queue[tail++]  = start;
    parent[start]  = -1;
    mu_bit_set(seen, (size_t)start);
    while (head < tail) {
        const int u = queue[head++];
        ++visited;
        if (!fn(u, parent[u], ctx)) {
            break;
        }
        size_t cur = 0;
        int    v;
        while (matrix_graph_next_neighbor(g, u, &cur, &v, NULL)) {
            if (!mu_bit_test(seen, (size_t)v)) {
                mu_bit_set(seen, (size_t)v);
                parent[v]     = u;
                queue[tail++] = v;
            }
        }
    }
    free(queue);
    free(parent);
    free(seen);
    return visited;
}

size_t matrix_dfs_visit(const MatrixGraph *g, int start, MatrixVisitFn fn, void *ctx)
{
    if (!mu_ok(g, start) || fn == NULL) {
        return 0;
    }
    const size_t n = matrix_graph_id_bound(g);
    int    *stack  = malloc(n * sizeof(int));
    size_t *cursor = malloc(n * sizeof(size_t));
    int    *parent = malloc(n * sizeof(int));
    uint64_t *seen = mu_calloc(mu_words(n), sizeof(uint64_t));
    if (stack == NULL || cursor == NULL || parent == NULL || seen == NULL) {
        free(stack);
        free(cursor);
        free(parent);
        free(seen);
        return 0;
    }
    size_t depth = 0, visited = 1;
    stack[0]  = start;
    cursor[0] = 0;
    parent[start] = -1;
    mu_bit_set(seen, (size_t)start);
    bool go = fn(start, -1, ctx);
    while (go) {
        const int u = stack[depth];
        int    v;
        bool   descended = false;
        while (matrix_graph_next_neighbor(g, u, &cursor[depth], &v, NULL)) {
            if (!mu_bit_test(seen, (size_t)v)) {
                mu_bit_set(seen, (size_t)v);
                parent[v] = u;
                ++visited;
                if (!fn(v, u, ctx)) {
                    go = false;
                    break;
                }
                ++depth;
                stack[depth]  = v;
                cursor[depth] = 0;
                descended = true;
                break;
            }
        }
        if (!go) {
            break;
        }
        if (!descended) {
            if (depth == 0) {
                break;
            }
            --depth;
        }
    }
    free(stack);
    free(cursor);
    free(parent);
    free(seen);
    return visited;
}

/* ------------------------------------------------------------------ */
/*  Centrality                                                          */
/* ------------------------------------------------------------------ */

MatrixStatus matrix_degree_centrality(const MatrixGraph *g, double *out)
{
    if (g == NULL || out == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n = matrix_graph_id_bound(g);
    const size_t v = matrix_graph_vertex_count(g);
    for (size_t i = 0; i < n; ++i) {
        out[i] = (v > 1 && matrix_graph_is_alive(g, (int)i))
                     ? (double)matrix_graph_degree(g, (int)i) / (double)(v - 1)
                     : 0.0;
    }
    return GRAPH_OK;
}

MatrixStatus matrix_closeness_centrality(const MatrixGraph *g, double *out)
{
    if (g == NULL || out == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n = matrix_graph_id_bound(g);
    const size_t vc = matrix_graph_vertex_count(g);
    int *dist  = malloc((n > 0 ? n : 1) * sizeof(int));
    int *queue = malloc((n > 0 ? n : 1) * sizeof(int));
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
        const size_t reached = mu_bfs_hops(g, (int)s, dist, queue) - 1;
        double sum = 0.0;
        for (size_t i = 1; i <= reached; ++i) {
            sum += dist[queue[i]];
        }
        if (reached > 0 && sum > 0.0) {
            out[s] = ((double)reached / (double)(vc - 1)) * ((double)reached / sum);
        }
    }
    free(dist);
    free(queue);
    return GRAPH_OK;
}

MatrixStatus matrix_betweenness_centrality(const MatrixGraph *g, double *out)
{
    if (g == NULL || out == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n = matrix_graph_id_bound(g);
    int    *dist  = malloc((n > 0 ? n : 1) * sizeof(int));
    int    *order = malloc((n > 0 ? n : 1) * sizeof(int));
    double *sigma = malloc((n > 0 ? n : 1) * sizeof(double));
    double *delta = malloc((n > 0 ? n : 1) * sizeof(double));
    if (dist == NULL || order == NULL || sigma == NULL || delta == NULL) {
        free(dist);
        free(order);
        free(sigma);
        free(delta);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    mu_fill_double(out, n, 0.0);

    for (size_t si = 0; si < n; ++si) {
        const int s = (int)si;
        if (!matrix_graph_is_alive(g, s)) {
            continue;
        }
        /* Forward BFS counting shortest paths. */
        mu_fill_int(dist, n, -1);
        mu_fill_double(sigma, n, 0.0);
        size_t head = 0, tail = 0;
        dist[s]  = 0;
        sigma[s] = 1.0;
        order[tail++] = s;
        while (head < tail) {
            const int u = order[head++];
            size_t cur = 0;
            int    v;
            while (matrix_graph_next_neighbor(g, u, &cur, &v, NULL)) {
                if (dist[v] < 0) {
                    dist[v]       = dist[u] + 1;
                    order[tail++] = v;
                }
                if (dist[v] == dist[u] + 1) {
                    sigma[v] += sigma[u];
                }
            }
        }
        /* Back-propagate dependencies in reverse BFS order, using successors. */
        for (size_t i = tail; i-- > 0;) {
            const int u = order[i];
            double d = 0.0;
            size_t cur = 0;
            int    v;
            while (matrix_graph_next_neighbor(g, u, &cur, &v, NULL)) {
                if (dist[v] == dist[u] + 1) {
                    d += sigma[u] / sigma[v] * (1.0 + delta[v]);
                }
            }
            delta[u] = d;
            if (u != s) {
                out[u] += d;
            }
        }
    }
    if (!matrix_graph_is_directed(g)) {
        for (size_t i = 0; i < n; ++i) {
            out[i] /= 2.0;
        }
    }
    free(dist);
    free(order);
    free(sigma);
    free(delta);
    return GRAPH_OK;
}

/* ------------------------------------------------------------------ */
/*  Triangles and clustering                                            */
/* ------------------------------------------------------------------ */

/** Neighbours of v other than v itself, ascending, into nb; returns the count. */
static size_t simple_neighbors(const MatrixGraph *g, int v, int *nb)
{
    size_t k = 0, cur = 0;
    int    w;
    while (matrix_graph_next_neighbor(g, v, &cur, &w, NULL)) {
        if (w != v) {
            nb[k++] = w;
        }
    }
    return k;
}

MatrixStatus matrix_triangle_count(const MatrixGraph *g, long long *total)
{
    if (g == NULL || total == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    if (matrix_graph_is_directed(g)) {
        return GRAPH_ERR_UNSUPPORTED;
    }
    const size_t n = matrix_graph_id_bound(g);
    int *nb = malloc((n > 0 ? n : 1) * sizeof(int));
    if (nb == NULL) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    long long count = 0;
    for (size_t ui = 0; ui < n; ++ui) {
        const int u = (int)ui;
        if (!matrix_graph_is_alive(g, u)) {
            continue;
        }
        const size_t k = simple_neighbors(g, u, nb);
        /* Only pairs of neighbours both above u, so each triangle counts once. */
        size_t first = 0;
        while (first < k && nb[first] < u) {
            ++first;
        }
        for (size_t i = first; i < k; ++i) {
            for (size_t j = i + 1; j < k; ++j) {
                if (matrix_graph_has_edge(g, nb[i], nb[j], NULL)) {
                    ++count;
                }
            }
        }
    }
    free(nb);
    *total = count;
    return GRAPH_OK;
}

MatrixStatus matrix_clustering_coefficient(const MatrixGraph *g, double *local, double *average)
{
    if (g == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    if (matrix_graph_is_directed(g)) {
        return GRAPH_ERR_UNSUPPORTED;
    }
    const size_t n = matrix_graph_id_bound(g);
    int *nb = malloc((n > 0 ? n : 1) * sizeof(int));
    if (nb == NULL) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    double sum = 0.0;
    size_t alive = 0;
    for (size_t ui = 0; ui < n; ++ui) {
        const int u = (int)ui;
        double c = 0.0;
        if (matrix_graph_is_alive(g, u)) {
            ++alive;
            const size_t k = simple_neighbors(g, u, nb);
            if (k >= 2) {
                size_t links = 0;
                for (size_t i = 0; i < k; ++i) {
                    for (size_t j = i + 1; j < k; ++j) {
                        if (matrix_graph_has_edge(g, nb[i], nb[j], NULL)) {
                            ++links;
                        }
                    }
                }
                c = 2.0 * (double)links / ((double)k * (double)(k - 1));
            }
            sum += c;
        }
        if (local != NULL) {
            local[ui] = c;
        }
    }
    if (average != NULL) {
        *average = alive > 0 ? sum / (double)alive : 0.0;
    }
    free(nb);
    return GRAPH_OK;
}

/* ------------------------------------------------------------------ */
/*  Bipartite                                                           */
/* ------------------------------------------------------------------ */

MatrixStatus matrix_is_bipartite(const MatrixGraph *g, bool *result, int *color)
{
    if (g == NULL || result == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n = matrix_graph_id_bound(g);
    MuSym sym;
    if (!mu_sym_build_ex(g, &sym, false)) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    int *col   = malloc((n > 0 ? n : 1) * sizeof(int));
    int *queue = malloc((n > 0 ? n : 1) * sizeof(int));
    if (col == NULL || queue == NULL) {
        free(col);
        free(queue);
        mu_sym_free(&sym);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    mu_fill_int(col, n, -1);
    bool ok = true;
    for (size_t s = 0; s < n && ok; ++s) {
        if (!matrix_graph_is_alive(g, (int)s)) {
            continue;
        }
        if (matrix_graph_has_edge(g, (int)s, (int)s, NULL)) {
            ok = false;     /* a self-loop is an odd cycle */
            break;
        }
        if (col[s] >= 0) {
            continue;
        }
        size_t head = 0, tail = 0;
        col[s] = 0;
        queue[tail++] = (int)s;
        while (head < tail && ok) {
            const int u = queue[head++];
            for (size_t k = sym.off[u]; k < sym.off[u + 1]; ++k) {
                const int v = sym.to[k];
                if (col[v] < 0) {
                    col[v] = 1 - col[u];
                    queue[tail++] = v;
                } else if (col[v] == col[u]) {
                    ok = false;
                    break;
                }
            }
        }
    }
    *result = ok;
    if (color != NULL) {
        memcpy(color, col, n * sizeof(int));
    }
    free(col);
    free(queue);
    mu_sym_free(&sym);
    return GRAPH_OK;
}

/* ------------------------------------------------------------------ */
/*  Eccentricity, diameter, radius                                      */
/* ------------------------------------------------------------------ */

MatrixStatus matrix_eccentricity(const MatrixGraph *g, int *ecc)
{
    if (g == NULL || ecc == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n = matrix_graph_id_bound(g);
    int *dist  = malloc((n > 0 ? n : 1) * sizeof(int));
    int *queue = malloc((n > 0 ? n : 1) * sizeof(int));
    if (dist == NULL || queue == NULL) {
        free(dist);
        free(queue);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    for (size_t s = 0; s < n; ++s) {
        if (!matrix_graph_is_alive(g, (int)s)) {
            ecc[s] = -1;
            continue;
        }
        const size_t reached = mu_bfs_hops(g, (int)s, dist, queue);
        ecc[s] = dist[queue[reached - 1]];     /* BFS order: last is farthest */
    }
    free(dist);
    free(queue);
    return GRAPH_OK;
}

MatrixStatus matrix_diameter_radius(const MatrixGraph *g, int *diameter, int *radius)
{
    if (g == NULL || diameter == NULL || radius == NULL || matrix_graph_vertex_count(g) == 0) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n  = matrix_graph_id_bound(g);
    const size_t vc = matrix_graph_vertex_count(g);
    int *dist  = malloc((n > 0 ? n : 1) * sizeof(int));
    int *queue = malloc((n > 0 ? n : 1) * sizeof(int));
    if (dist == NULL || queue == NULL) {
        free(dist);
        free(queue);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    int diam = 0, rad = -1;
    MatrixStatus st = GRAPH_OK;
    for (size_t s = 0; s < n; ++s) {
        if (!matrix_graph_is_alive(g, (int)s)) {
            continue;
        }
        const size_t reached = mu_bfs_hops(g, (int)s, dist, queue);
        if (reached != vc) {
            st = GRAPH_ERR_UNSUPPORTED;
            break;
        }
        const int e = dist[queue[reached - 1]];
        if (e > diam) {
            diam = e;
        }
        if (rad < 0 || e < rad) {
            rad = e;
        }
    }
    free(dist);
    free(queue);
    if (st == GRAPH_OK) {
        *diameter = diam;
        *radius   = rad;
    }
    return st;
}

/* ------------------------------------------------------------------ */
/*  Transpose, neighbourhood similarity, summary statistics             */
/* ------------------------------------------------------------------ */

MatrixStatus matrix_transpose(const MatrixGraph *g, MatrixGraph **out)
{
    if (g == NULL || out == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    if (!matrix_graph_is_directed(g)) {
        *out = matrix_graph_copy(g);
        return *out != NULL ? GRAPH_OK : GRAPH_ERR_OUT_OF_MEMORY;
    }
    const size_t n = matrix_graph_id_bound(g);
    MatrixGraph *t = matrix_graph_create(n, true);
    if (t == NULL) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    for (size_t u = 0; u < n; ++u) {
        if (!matrix_graph_is_alive(g, (int)u) && !matrix_graph_remove_vertex(t, (int)u)) {
            matrix_graph_destroy(t);
            return GRAPH_ERR_OUT_OF_MEMORY;
        }
    }
    for (size_t u = 0; u < n; ++u) {
        size_t cur = 0;
        int    v;
        double w;
        while (matrix_graph_is_alive(g, (int)u) &&
               matrix_graph_next_neighbor(g, (int)u, &cur, &v, &w)) {
            if (!matrix_graph_add_edge(t, v, (int)u, w)) {
                matrix_graph_destroy(t);
                return GRAPH_ERR_OUT_OF_MEMORY;
            }
        }
    }
    *out = t;
    return GRAPH_OK;
}

MatrixStatus matrix_common_neighbors(const MatrixGraph *g, int u, int v, int *out, size_t *count)
{
    if (!mu_ok(g, u) || !mu_ok(g, v) || out == NULL || count == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    size_t k = 0, cur = 0;
    int    w;
    while (matrix_graph_next_neighbor(g, u, &cur, &w, NULL)) {
        if (matrix_graph_has_edge(g, v, w, NULL)) {
            out[k++] = w;
        }
    }
    *count = k;
    return GRAPH_OK;
}

MatrixStatus matrix_jaccard(const MatrixGraph *g, int u, int v, double *score)
{
    if (!mu_ok(g, u) || !mu_ok(g, v) || score == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    size_t inter = 0, cur = 0;
    int    w;
    while (matrix_graph_next_neighbor(g, u, &cur, &w, NULL)) {
        if (matrix_graph_has_edge(g, v, w, NULL)) {
            ++inter;
        }
    }
    const size_t uni = matrix_graph_degree(g, u) + matrix_graph_degree(g, v) - inter;
    *score = uni == 0 ? 0.0 : (double)inter / (double)uni;
    return GRAPH_OK;
}

MatrixStatus matrix_network_stats(const MatrixGraph *g, MatrixNetworkStats *out)
{
    if (g == NULL || out == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n  = matrix_graph_id_bound(g);
    const size_t vc = matrix_graph_vertex_count(g);
    const bool directed = matrix_graph_is_directed(g);
    size_t arcs = 0, reciprocal = 0, maxdeg = 0;
    for (size_t ui = 0; ui < n; ++ui) {
        const int u = (int)ui;
        if (!matrix_graph_is_alive(g, u)) {
            continue;
        }
        const size_t d = matrix_graph_degree(g, u);
        if (d > maxdeg) {
            maxdeg = d;
        }
        size_t cur = 0;
        int    v;
        while (matrix_graph_next_neighbor(g, u, &cur, &v, NULL)) {
            if (v != u) {              /* self-loops count towards degree only */
                ++arcs;
                if (directed && matrix_graph_has_edge(g, v, u, NULL)) {
                    ++reciprocal;
                }
            }
        }
    }
    const double edges = directed ? (double)arcs : (double)arcs / 2.0;
    const double possible = directed ? (double)vc * (double)(vc > 0 ? vc - 1 : 0)
                                     : (double)vc * (double)(vc > 0 ? vc - 1 : 0) / 2.0;
    out->density     = possible > 0.0 ? edges / possible : 0.0;
    out->reciprocity = !directed ? 1.0 : (arcs > 0 ? (double)reciprocal / (double)arcs : 0.0);
    out->max_degree  = maxdeg;
    return GRAPH_OK;
}

MatrixStatus matrix_degree_histogram(const MatrixGraph *g, size_t *hist, size_t len)
{
    if (g == NULL || hist == NULL || len == 0) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    memset(hist, 0, len * sizeof(size_t));
    const size_t n = matrix_graph_id_bound(g);
    for (size_t i = 0; i < n; ++i) {
        if (!matrix_graph_is_alive(g, (int)i)) {
            continue;
        }
        const size_t d = matrix_graph_degree(g, (int)i);
        ++hist[d >= len ? len - 1 : d];
    }
    return GRAPH_OK;
}
