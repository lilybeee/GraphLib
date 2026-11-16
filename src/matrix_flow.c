/*
 * matrix_flow.c — Prim's MST, maximum flow / minimum cut, bipartite matching
 * and the Hungarian assignment.  See matrix_algorithms.h.
 */

#include "matrix_util.h"

/* ------------------------------------------------------------------ */
/*  Prim                                                                */
/* ------------------------------------------------------------------ */

double matrix_mst_prim(const MatrixGraph *g, MatrixEdge *edges, size_t *edge_count)
{
    if (g == NULL || edge_count == NULL) {
        return INFINITY;
    }
    *edge_count = 0;
    if (matrix_graph_is_directed(g)) {
        return INFINITY;
    }
    const size_t n  = matrix_graph_id_bound(g);
    const size_t vc = matrix_graph_vertex_count(g);
    if (vc <= 1) {
        return 0.0;
    }
    if (edges == NULL) {
        return INFINITY;
    }
    double *key    = malloc(n * sizeof(double));
    int    *from   = malloc(n * sizeof(int));
    bool   *intree = calloc(n, sizeof(bool));
    if (key == NULL || from == NULL || intree == NULL) {
        free(key);
        free(from);
        free(intree);
        return INFINITY;
    }
    mu_fill_double(key, n, INFINITY);
    mu_fill_int(from, n, -1);

    int start = -1;
    for (size_t i = 0; i < n && start < 0; ++i) {
        if (matrix_graph_is_alive(g, (int)i)) {
            start = (int)i;
        }
    }
    key[start] = 0.0;
    double total = 0.0;
    size_t added = 0;
    bool   ok    = true;
    for (size_t step = 0; step < vc && ok; ++step) {
        int best = -1;
        for (size_t i = 0; i < n; ++i) {
            if (!intree[i] && key[i] < INFINITY && (best < 0 || key[i] < key[best])) {
                best = (int)i;
            }
        }
        if (best < 0) {
            ok = false;     /* the rest is unreachable: disconnected */
            break;
        }
        intree[best] = true;
        if (from[best] >= 0) {
            edges[added].u      = from[best] < best ? from[best] : best;
            edges[added].v      = from[best] < best ? best : from[best];
            edges[added].weight = key[best];
            total += key[best];
            ++added;
        }
        size_t cur = 0;
        int    v;
        double w;
        while (matrix_graph_next_neighbor(g, best, &cur, &v, &w)) {
            if (!intree[v] && w < key[v]) {
                key[v]  = w;
                from[v] = best;
            }
        }
    }
    free(key);
    free(from);
    free(intree);
    if (!ok) {
        return INFINITY;
    }
    *edge_count = added;
    return total;
}

/* ------------------------------------------------------------------ */
/*  Residual network                                                    */
/* ------------------------------------------------------------------ */

typedef struct {
    size_t  n;
    int    *head;     /* first arc out of each vertex, -1 if none */
    int    *nxt;
    int    *to;
    double *cap;      /* residual capacity; arc a's reverse is a ^ 1 */
    size_t  arcs;
    double  eps;
} Residual;

static void res_free(Residual *r)
{
    free(r->head);
    free(r->nxt);
    free(r->to);
    free(r->cap);
    memset(r, 0, sizeof(*r));
}

static void res_add(Residual *r, int u, int v, double fwd, double back)
{
    const int a = (int)r->arcs;
    r->to[a] = v;   r->cap[a] = fwd;   r->nxt[a] = r->head[u];   r->head[u] = a;
    r->to[a + 1] = u;   r->cap[a + 1] = back;   r->nxt[a + 1] = r->head[v];   r->head[v] = a + 1;
    r->arcs += 2;
}

static MatrixStatus res_build(const MatrixGraph *g, Residual *r)
{
    memset(r, 0, sizeof(*r));
    const size_t n = matrix_graph_id_bound(g);
    const size_t e = matrix_graph_edge_count(g);
    const bool directed = matrix_graph_is_directed(g);
    r->n    = n;
    r->head = malloc((n > 0 ? n : 1) * sizeof(int));
    r->nxt  = malloc((2 * e + 2) * sizeof(int));
    r->to   = malloc((2 * e + 2) * sizeof(int));
    r->cap  = malloc((2 * e + 2) * sizeof(double));
    if (r->head == NULL || r->nxt == NULL || r->to == NULL || r->cap == NULL) {
        res_free(r);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    mu_fill_int(r->head, n, -1);
    double maxcap = 1.0;
    for (size_t ui = 0; ui < n; ++ui) {
        if (!matrix_graph_is_alive(g, (int)ui)) {
            continue;
        }
        size_t cur = 0;
        int    v;
        double w;
        while (matrix_graph_next_neighbor(g, (int)ui, &cur, &v, &w)) {
            if (w < 0.0) {
                res_free(r);
                return GRAPH_ERR_INVALID_ARGUMENT;
            }
            if (v == (int)ui || (!directed && v < (int)ui)) {
                continue;       /* loops carry no flow; undirected edges once */
            }
            if (w > maxcap) {
                maxcap = w;
            }
            res_add(r, (int)ui, v, w, directed ? 0.0 : w);
        }
    }
    r->eps = 1e-12 * maxcap;
    return GRAPH_OK;
}

/** BFS levels over arcs with residual capacity; returns true if t is reachable. */
static bool res_levels(const Residual *r, int s, int t, int *level, int *queue)
{
    mu_fill_int(level, r->n, -1);
    size_t head = 0, tail = 0;
    level[s] = 0;
    queue[tail++] = s;
    while (head < tail) {
        const int u = queue[head++];
        for (int a = r->head[u]; a >= 0; a = r->nxt[a]) {
            if (r->cap[a] > r->eps && level[r->to[a]] < 0) {
                level[r->to[a]] = level[u] + 1;
                queue[tail++] = r->to[a];
            }
        }
    }
    return level[t] >= 0;
}

static double dinic(Residual *r, int s, int t, int *level, int *queue)
{
    int *it    = malloc(r->n * sizeof(int));
    int *stack = malloc((r->n + 1) * sizeof(int));     /* arcs on the current path */
    if (it == NULL || stack == NULL) {
        free(it);
        free(stack);
        return -1.0;
    }
    double flow = 0.0;
    while (res_levels(r, s, t, level, queue)) {
        memcpy(it, r->head, r->n * sizeof(int));
        size_t sp = 0;
        int    u  = s;
        for (;;) {
            if (u == t) {
                double push = INFINITY;
                for (size_t i = 0; i < sp; ++i) {
                    if (r->cap[stack[i]] < push) {
                        push = r->cap[stack[i]];
                    }
                }
                for (size_t i = 0; i < sp; ++i) {
                    r->cap[stack[i]]     -= push;
                    r->cap[stack[i] ^ 1] += push;
                }
                flow += push;
                sp = 0;
                u  = s;
                continue;
            }
            bool advanced = false;
            while (it[u] >= 0) {
                const int a = it[u];
                if (r->cap[a] > r->eps && level[r->to[a]] == level[u] + 1) {
                    stack[sp++] = a;
                    u = r->to[a];
                    advanced = true;
                    break;
                }
                it[u] = r->nxt[a];
            }
            if (advanced) {
                continue;
            }
            if (u == s) {
                break;              /* blocking flow reached */
            }
            const int a = stack[--sp];
            u = r->to[a ^ 1];
            it[u] = r->nxt[it[u]];  /* that arc leads nowhere useful */
        }
    }
    free(it);
    free(stack);
    return flow;
}

static double edmonds_karp(Residual *r, int s, int t, int *via, int *queue)
{
    double flow = 0.0;
    for (;;) {
        mu_fill_int(via, r->n, -2);
        size_t head = 0, tail = 0;
        via[s] = -1;
        queue[tail++] = s;
        while (head < tail && via[t] == -2) {
            const int u = queue[head++];
            for (int a = r->head[u]; a >= 0; a = r->nxt[a]) {
                if (r->cap[a] > r->eps && via[r->to[a]] == -2) {
                    via[r->to[a]] = a;
                    queue[tail++] = r->to[a];
                }
            }
        }
        if (via[t] == -2) {
            return flow;
        }
        double push = INFINITY;
        for (int v = t; v != s; v = r->to[via[v] ^ 1]) {
            if (r->cap[via[v]] < push) {
                push = r->cap[via[v]];
            }
        }
        for (int v = t; v != s; v = r->to[via[v] ^ 1]) {
            r->cap[via[v]]     -= push;
            r->cap[via[v] ^ 1] += push;
        }
        flow += push;
    }
}

MatrixStatus matrix_max_flow(const MatrixGraph *g, int source, int sink, MatrixFlowAlgo algo,
                             double *value, bool *source_side)
{
    if (!mu_ok(g, source) || !mu_ok(g, sink) || source == sink || value == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    Residual r;
    MatrixStatus st = res_build(g, &r);
    if (st != GRAPH_OK) {
        return st;
    }
    int *level = malloc((r.n > 0 ? r.n : 1) * sizeof(int));
    int *queue = malloc((r.n > 0 ? r.n : 1) * sizeof(int));
    if (level == NULL || queue == NULL) {
        free(level);
        free(queue);
        res_free(&r);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    const double f = algo == MATRIX_FLOW_EDMONDS_KARP ? edmonds_karp(&r, source, sink, level, queue)
                                                      : dinic(&r, source, sink, level, queue);
    if (f < 0.0) {
        st = GRAPH_ERR_OUT_OF_MEMORY;
    } else {
        *value = f;
        if (source_side != NULL) {
            (void)res_levels(&r, source, sink, level, queue);
            for (size_t i = 0; i < r.n; ++i) {
                source_side[i] = level[i] >= 0;
            }
        }
    }
    free(level);
    free(queue);
    res_free(&r);
    return st;
}

MatrixStatus matrix_min_cut(const MatrixGraph *g, int source, int sink, double *value,
                            MatrixEdge *cut, size_t *count)
{
    if (!mu_ok(g, source) || !mu_ok(g, sink) || value == NULL || cut == NULL || count == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n = matrix_graph_id_bound(g);
    bool *side = malloc((n > 0 ? n : 1) * sizeof(bool));
    if (side == NULL) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    MatrixStatus st = matrix_max_flow(g, source, sink, MATRIX_FLOW_DINIC, value, side);
    if (st == GRAPH_OK) {
        size_t k = 0;
        for (size_t ui = 0; ui < n; ++ui) {
            if (!side[ui] || !matrix_graph_is_alive(g, (int)ui)) {
                continue;
            }
            size_t cur = 0;
            int    v;
            double w;
            while (matrix_graph_next_neighbor(g, (int)ui, &cur, &v, &w)) {
                if (!side[v]) {
                    cut[k].u      = (int)ui;
                    cut[k].v      = v;
                    cut[k].weight = w;
                    ++k;
                }
            }
        }
        *count = k;
    }
    free(side);
    return st;
}

/* ------------------------------------------------------------------ */
/*  Hopcroft-Karp                                                       */
/* ------------------------------------------------------------------ */

MatrixStatus matrix_bipartite_matching(const MatrixGraph *g, int *match, size_t *size)
{
    if (g == NULL || match == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n = matrix_graph_id_bound(g);
    int *color = malloc((n > 0 ? n : 1) * sizeof(int));
    if (color == NULL) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    bool bip = false;
    MatrixStatus st = matrix_is_bipartite(g, &bip, color);
    if (st != GRAPH_OK || !bip) {
        free(color);
        return st != GRAPH_OK ? st : GRAPH_ERR_UNSUPPORTED;
    }
    MuSym sym;
    if (!mu_sym_build_ex(g, &sym, false)) {
        free(color);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    const int inf = (int)(n + 1);
    int    *dist  = malloc((n > 0 ? n : 1) * sizeof(int));
    int    *queue = malloc((n > 0 ? n : 1) * sizeof(int));
    int    *stack = malloc((n + 1) * sizeof(int));
    int    *via   = malloc((n + 1) * sizeof(int));
    size_t *cur   = malloc((n > 0 ? n : 1) * sizeof(size_t));
    if (dist == NULL || queue == NULL || stack == NULL || via == NULL || cur == NULL) {
        free(color); free(dist); free(queue); free(stack); free(via); free(cur);
        mu_sym_free(&sym);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    mu_fill_int(match, n, -1);
    size_t matched = 0;

    for (;;) {
        /* BFS from every free left vertex, layering by alternating paths. */
        size_t head = 0, tail = 0;
        for (size_t u = 0; u < n; ++u) {
            dist[u] = inf;
            if (color[u] == 0 && match[u] < 0) {
                dist[u] = 0;
                queue[tail++] = (int)u;
            }
        }
        bool reachable_free = false;
        while (head < tail) {
            const int u = queue[head++];
            for (size_t k = sym.off[u]; k < sym.off[u + 1]; ++k) {
                const int m = match[sym.to[k]];
                if (m < 0) {
                    reachable_free = true;
                } else if (dist[m] == inf) {
                    dist[m] = dist[u] + 1;
                    queue[tail++] = m;
                }
            }
        }
        if (!reachable_free) {
            break;
        }
        /* DFS: vertex-disjoint shortest augmenting paths. */
        for (size_t u = 0; u < n; ++u) {
            cur[u] = sym.off[u];
        }
        size_t grew = 0;
        for (size_t s = 0; s < n; ++s) {
            if (color[s] != 0 || match[s] >= 0) {
                continue;
            }
            size_t sp = 0;
            stack[sp++] = (int)s;
            while (sp > 0) {
                const int x = stack[sp - 1];
                if (cur[x] < sym.off[x + 1]) {
                    const int v = sym.to[cur[x]++];
                    const int m = match[v];
                    via[sp - 1] = v;
                    if (m < 0) {
                        for (size_t i = sp; i-- > 0;) {
                            match[stack[i]] = via[i];
                            match[via[i]]   = stack[i];
                        }
                        ++grew;
                        break;
                    }
                    if (dist[m] == dist[x] + 1) {
                        stack[sp++] = m;
                    }
                } else {
                    dist[x] = inf;      /* dead end for this phase */
                    --sp;
                }
            }
        }
        if (grew == 0) {
            break;
        }
        matched += grew;
    }
    if (size != NULL) {
        *size = matched;
    }
    free(color); free(dist); free(queue); free(stack); free(via); free(cur);
    mu_sym_free(&sym);
    return GRAPH_OK;
}

/* ------------------------------------------------------------------ */
/*  Hungarian                                                           */
/* ------------------------------------------------------------------ */

MatrixStatus matrix_hungarian(const double *cost, size_t rows, size_t cols, int *assignment,
                              double *total)
{
    if (cost == NULL || assignment == NULL || rows == 0 || cols < rows) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    for (size_t i = 0; i < rows * cols; ++i) {
        if (!isfinite(cost[i])) {
            return GRAPH_ERR_INVALID_ARGUMENT;
        }
    }
    const size_t n = rows, m = cols;
    double *u = calloc(n + 1, sizeof(double));
    double *v = calloc(m + 1, sizeof(double));
    int    *p = calloc(m + 1, sizeof(int));
    int    *way = calloc(m + 1, sizeof(int));
    double *minv = malloc((m + 1) * sizeof(double));
    bool   *used = malloc((m + 1) * sizeof(bool));
    if (u == NULL || v == NULL || p == NULL || way == NULL || minv == NULL || used == NULL) {
        free(u); free(v); free(p); free(way); free(minv); free(used);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    for (size_t i = 1; i <= n; ++i) {
        p[0] = (int)i;
        size_t j0 = 0;
        for (size_t j = 0; j <= m; ++j) {
            minv[j] = INFINITY;
            used[j] = false;
        }
        do {
            used[j0] = true;
            const size_t i0 = (size_t)p[j0];
            double delta = INFINITY;
            size_t j1 = 0;
            for (size_t j = 1; j <= m; ++j) {
                if (used[j]) {
                    continue;
                }
                const double cur = cost[(i0 - 1) * m + (j - 1)] - u[i0] - v[j];
                if (cur < minv[j]) {
                    minv[j] = cur;
                    way[j]  = (int)j0;
                }
                if (minv[j] < delta) {
                    delta = minv[j];
                    j1    = j;
                }
            }
            for (size_t j = 0; j <= m; ++j) {
                if (used[j]) {
                    u[p[j]] += delta;
                    v[j]    -= delta;
                } else {
                    minv[j] -= delta;
                }
            }
            j0 = j1;
        } while (p[j0] != 0);
        do {
            const size_t j1 = (size_t)way[j0];
            p[j0] = p[j1];
            j0 = j1;
        } while (j0 != 0);
    }
    double sum = 0.0;
    for (size_t j = 1; j <= m; ++j) {
        if (p[j] != 0) {
            assignment[p[j] - 1] = (int)(j - 1);
            sum += cost[(size_t)(p[j] - 1) * m + (j - 1)];
        }
    }
    if (total != NULL) {
        *total = sum;
    }
    free(u); free(v); free(p); free(way); free(minv); free(used);
    return GRAPH_OK;
}
