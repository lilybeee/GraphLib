/*
 * list_paths.c — shortest-path algorithms.  See list_algorithms.h.
 */

#include "list_util.h"

#define FW_MAX_VERTICES 8192u

/* ------------------------------------------------------------------ */
/*  Path reconstruction                                                 */
/* ------------------------------------------------------------------ */

size_t list_path_from_parent(const int *parent, size_t n, int src, int dst, int *path, size_t cap)
{
    if (parent == NULL || path == NULL || src < 0 || dst < 0 || (size_t)src >= n ||
        (size_t)dst >= n) {
        return 0;
    }
    size_t len = 1;
    int    v   = dst;
    while (v != src) {
        v = parent[v];
        if (v < 0 || len > n) {
            return 0;
        }
        ++len;
    }
    if (len > cap) {
        return 0;
    }
    v = dst;
    for (size_t i = len; i-- > 0;) {
        path[i] = v;
        v = parent[v];
    }
    return len;
}

size_t list_path_from_next(const int *next, size_t n, int u, int v, int *path, size_t cap)
{
    if (next == NULL || path == NULL || u < 0 || v < 0 || (size_t)u >= n || (size_t)v >= n) {
        return 0;
    }
    size_t len = 0;
    int    x   = u;
    for (;;) {
        if (len >= cap || len > n) {
            return 0;
        }
        path[len++] = x;
        if (x == v) {
            return len;
        }
        x = next[(size_t)x * n + (size_t)v];
        if (x < 0) {
            return 0;
        }
    }
}

/* ------------------------------------------------------------------ */
/*  Bellman-Ford                                                        */
/* ------------------------------------------------------------------ */

ListStatus list_bellman_ford(const ListGraph *g, int src, double *distance, int *parent)
{
    if (!mu_ok(g, src) || distance == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n  = list_graph_id_bound(g);
    const size_t vc = list_graph_vertex_count(g);
    mu_fill_double(distance, n, INFINITY);
    if (parent != NULL) {
        mu_fill_int(parent, n, -1);
    }
    distance[src] = 0.0;

    bool   changed = true;
    size_t rounds  = 0;
    while (changed && rounds < vc) {
        changed = false;
        for (size_t ui = 0; ui < n; ++ui) {
            if (isinf(distance[ui]) || !list_graph_is_alive(g, (int)ui)) {
                continue;
            }
            size_t cur = 0;
            int    v;
            double w;
            while (list_graph_next_neighbor(g, (int)ui, &cur, &v, &w)) {
                if (distance[ui] + w < distance[v]) {
                    distance[v] = distance[ui] + w;
                    if (parent != NULL) {
                        parent[v] = (int)ui;
                    }
                    changed = true;
                }
            }
        }
        ++rounds;
    }
    return changed ? GRAPH_ERR_NEGATIVE_CYCLE : GRAPH_OK;
}

/* ------------------------------------------------------------------ */
/*  All pairs                                                           */
/* ------------------------------------------------------------------ */

ListStatus list_floyd_warshall(const ListGraph *g, double *distance, int *next)
{
    if (g == NULL || distance == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n = list_graph_id_bound(g);
    if (n > FW_MAX_VERTICES) {
        return GRAPH_ERR_TOO_LARGE;
    }
    mu_fill_double(distance, n * n, INFINITY);
    if (next != NULL) {
        mu_fill_int(next, n * n, -1);
    }
    for (size_t u = 0; u < n; ++u) {
        if (!list_graph_is_alive(g, (int)u)) {
            continue;
        }
        distance[u * n + u] = 0.0;
        size_t cur = 0;
        int    v;
        double w;
        while (list_graph_next_neighbor(g, (int)u, &cur, &v, &w)) {
            if (w < distance[u * n + (size_t)v]) {
                distance[u * n + (size_t)v] = w;
                if (next != NULL) {
                    next[u * n + (size_t)v] = v;
                }
            }
        }
    }
    for (size_t k = 0; k < n; ++k) {
        if (!list_graph_is_alive(g, (int)k)) {
            continue;
        }
        const double *rk = distance + k * n;
        for (size_t i = 0; i < n; ++i) {
            const double dik = distance[i * n + k];
            if (isinf(dik)) {
                continue;
            }
            double *ri = distance + i * n;
            for (size_t j = 0; j < n; ++j) {
                const double nd = dik + rk[j];
                if (nd < ri[j]) {
                    ri[j] = nd;
                    if (next != NULL) {
                        next[i * n + j] = next[i * n + k];
                    }
                }
            }
        }
    }
    for (size_t i = 0; i < n; ++i) {
        if (list_graph_is_alive(g, (int)i) && distance[i * n + i] < 0.0) {
            return GRAPH_ERR_NEGATIVE_CYCLE;
        }
    }
    return GRAPH_OK;
}

ListStatus list_johnson(const ListGraph *g, double *distance, int *next)
{
    if (g == NULL || distance == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n  = list_graph_id_bound(g);
    const size_t vc = list_graph_vertex_count(g);
    if (n > FW_MAX_VERTICES) {
        return GRAPH_ERR_TOO_LARGE;
    }
    double *h     = calloc(n > 0 ? n : 1, sizeof(double));
    double *d     = malloc((n > 0 ? n : 1) * sizeof(double));
    int    *first = malloc((n > 0 ? n : 1) * sizeof(int));
    MuHeap heap;
    if (h == NULL || d == NULL || first == NULL || !mu_heap_init(&heap, 64)) {
        free(h);
        free(d);
        free(first);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }

    /* Potentials: Bellman-Ford from a virtual source joined to every vertex by a 0 edge. */
    bool   changed = true;
    size_t rounds  = 0;
    while (changed && rounds < vc) {
        changed = false;
        for (size_t u = 0; u < n; ++u) {
            if (!list_graph_is_alive(g, (int)u)) {
                continue;
            }
            size_t cur = 0;
            int    v;
            double w;
            while (list_graph_next_neighbor(g, (int)u, &cur, &v, &w)) {
                if (h[u] + w < h[v]) {
                    h[v]    = h[u] + w;
                    changed = true;
                }
            }
        }
        ++rounds;
    }
    ListStatus st = changed ? GRAPH_ERR_NEGATIVE_CYCLE : GRAPH_OK;

    mu_fill_double(distance, n * n, INFINITY);
    if (next != NULL) {
        mu_fill_int(next, n * n, -1);
    }
    for (size_t s = 0; s < n && st == GRAPH_OK; ++s) {
        if (!list_graph_is_alive(g, (int)s)) {
            continue;
        }
        mu_fill_double(d, n, INFINITY);
        mu_fill_int(first, n, -1);
        d[s] = 0.0;
        heap.size = 0;
        if (!mu_heap_push(&heap, 0.0, (int)s)) {
            st = GRAPH_ERR_OUT_OF_MEMORY;
            break;
        }
        while (heap.size > 0) {
            double key;
            int    u;
            mu_heap_pop(&heap, &key, &u);
            if (key > d[u]) {
                continue;
            }
            size_t cur = 0;
            int    v;
            double w;
            while (list_graph_next_neighbor(g, u, &cur, &v, &w)) {
                double rw = w + h[u] - h[v];
                if (rw < 0.0) {
                    rw = 0.0;     /* rounding only; true reduced weights are >= 0 */
                }
                if (d[u] + rw < d[v]) {
                    d[v]     = d[u] + rw;
                    first[v] = (u == (int)s) ? v : first[u];
                    if (!mu_heap_push(&heap, d[v], v)) {
                        st = GRAPH_ERR_OUT_OF_MEMORY;
                        break;
                    }
                }
            }
            if (st != GRAPH_OK) {
                break;
            }
        }
        for (size_t v = 0; v < n && st == GRAPH_OK; ++v) {
            if (!isinf(d[v])) {
                distance[s * n + v] = d[v] - h[s] + h[v];
                if (next != NULL) {
                    next[s * n + v] = first[v];
                }
            }
        }
        if (st == GRAPH_OK) {
            distance[s * n + s] = 0.0;
            if (next != NULL) {
                next[s * n + s] = -1;
            }
        }
    }
    mu_heap_free(&heap);
    free(h);
    free(d);
    free(first);
    return st;
}

/* ------------------------------------------------------------------ */
/*  Dijkstra family                                                     */
/* ------------------------------------------------------------------ */

ListStatus list_dijkstra_multi(const ListGraph *g, const int *sources, size_t count,
                                   double *distance, int *parent, int *origin)
{
    if (g == NULL || sources == NULL || count == 0 || distance == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    for (size_t i = 0; i < count; ++i) {
        if (!list_graph_is_alive(g, sources[i])) {
            return GRAPH_ERR_INVALID_ARGUMENT;
        }
    }
    const size_t n = list_graph_id_bound(g);
    mu_fill_double(distance, n, INFINITY);
    if (parent != NULL) {
        mu_fill_int(parent, n, -1);
    }
    if (origin != NULL) {
        mu_fill_int(origin, n, -1);
    }
    MuHeap heap;
    if (!mu_heap_init(&heap, 64)) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    ListStatus st = GRAPH_OK;
    for (size_t i = 0; i < count; ++i) {
        distance[sources[i]] = 0.0;
        if (origin != NULL) {
            origin[sources[i]] = sources[i];
        }
        if (!mu_heap_push(&heap, 0.0, sources[i])) {
            st = GRAPH_ERR_OUT_OF_MEMORY;
        }
    }
    while (st == GRAPH_OK && heap.size > 0) {
        double key;
        int    u;
        mu_heap_pop(&heap, &key, &u);
        if (key > distance[u]) {
            continue;
        }
        size_t cur = 0;
        int    v;
        double w;
        while (list_graph_next_neighbor(g, u, &cur, &v, &w)) {
            if (w < 0.0) {
                st = GRAPH_ERR_INVALID_ARGUMENT;
                break;
            }
            if (distance[u] + w < distance[v]) {
                distance[v] = distance[u] + w;
                if (parent != NULL) {
                    parent[v] = u;
                }
                if (origin != NULL) {
                    origin[v] = origin[u];
                }
                if (!mu_heap_push(&heap, distance[v], v)) {
                    st = GRAPH_ERR_OUT_OF_MEMORY;
                    break;
                }
            }
        }
    }
    mu_heap_free(&heap);
    return st;
}

ListStatus list_a_star(const ListGraph *g, int src, int dst, ListHeuristic h, void *ctx,
                           ListPath *out)
{
    if (!mu_ok(g, src) || !mu_ok(g, dst) || out == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n = list_graph_id_bound(g);
    double *gc     = malloc((n > 0 ? n : 1) * sizeof(double));
    int    *parent = malloc((n > 0 ? n : 1) * sizeof(int));
    MuHeap heap;
    if (gc == NULL || parent == NULL || !mu_heap_init(&heap, 64)) {
        free(gc);
        free(parent);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    mu_fill_double(gc, n, INFINITY);
    mu_fill_int(parent, n, -1);
    gc[src] = 0.0;
    ListStatus st = GRAPH_ERR_NOT_FOUND;
    if (!mu_heap_push(&heap, h != NULL ? h(src, dst, ctx) : 0.0, src)) {
        st = GRAPH_ERR_OUT_OF_MEMORY;
    }
    while (st == GRAPH_ERR_NOT_FOUND && heap.size > 0) {
        double key;
        int    u;
        mu_heap_pop(&heap, &key, &u);
        const double hu = h != NULL ? h(u, dst, ctx) : 0.0;
        if (key > gc[u] + hu + 1e-12 * fmax(1.0, fabs(key))) {
            continue;       /* a better route to u was found after this entry was pushed */
        }
        if (u == dst) {
            st = GRAPH_OK;
            break;
        }
        size_t cur = 0;
        int    v;
        double w;
        while (list_graph_next_neighbor(g, u, &cur, &v, &w)) {
            if (w < 0.0) {
                st = GRAPH_ERR_INVALID_ARGUMENT;
                break;
            }
            if (gc[u] + w < gc[v]) {
                gc[v]     = gc[u] + w;
                parent[v] = u;
                if (!mu_heap_push(&heap, gc[v] + (h != NULL ? h(v, dst, ctx) : 0.0), v)) {
                    st = GRAPH_ERR_OUT_OF_MEMORY;
                    break;
                }
            }
        }
    }
    mu_heap_free(&heap);
    if (st == GRAPH_OK) {
        int *buf = malloc(n * sizeof(int));
        if (buf == NULL) {
            st = GRAPH_ERR_OUT_OF_MEMORY;
        } else {
            out->length   = list_path_from_parent(parent, n, src, dst, buf, n);
            out->vertices = buf;
            out->cost     = gc[dst];
        }
    }
    free(gc);
    free(parent);
    return st;
}

ListStatus list_zero_one_bfs(const ListGraph *g, int src, double *distance, int *parent)
{
    if (!mu_ok(g, src) || distance == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n = list_graph_id_bound(g);
    size_t arcs = 0;
    for (size_t u = 0; u < n; ++u) {
        if (!list_graph_is_alive(g, (int)u)) {
            continue;
        }
        arcs += list_graph_degree(g, (int)u);
        size_t cur = 0;
        int    v;
        double w;
        while (list_graph_next_neighbor(g, (int)u, &cur, &v, &w)) {
            if (w != 0.0 && w != 1.0) {
                return GRAPH_ERR_INVALID_ARGUMENT;
            }
        }
    }
    const size_t cap = 2 * (arcs + 2) + 2;
    int    *dv = malloc(cap * sizeof(int));
    double *dd = malloc(cap * sizeof(double));
    if (dv == NULL || dd == NULL) {
        free(dv);
        free(dd);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    mu_fill_double(distance, n, INFINITY);
    if (parent != NULL) {
        mu_fill_int(parent, n, -1);
    }
    size_t head = cap / 2, tail = cap / 2;
    distance[src] = 0.0;
    dv[tail] = src;
    dd[tail] = 0.0;
    ++tail;
    while (head < tail) {
        const int    u  = dv[head];
        const double du = dd[head];
        ++head;
        if (du > distance[u]) {
            continue;
        }
        size_t cur = 0;
        int    v;
        double w;
        while (list_graph_next_neighbor(g, u, &cur, &v, &w)) {
            if (du + w < distance[v]) {
                distance[v] = du + w;
                if (parent != NULL) {
                    parent[v] = u;
                }
                if (w == 0.0) {
                    --head;
                    dv[head] = v;
                    dd[head] = distance[v];
                } else {
                    dv[tail] = v;
                    dd[tail] = distance[v];
                    ++tail;
                }
            }
        }
    }
    free(dv);
    free(dd);
    return GRAPH_OK;
}

/* ------------------------------------------------------------------ */
/*  Yen's k shortest loopless paths                                     */
/* ------------------------------------------------------------------ */

/**
 * Dijkstra from src to dst that may not enter banned vertices and, when it
 * stands on `spur`, may not take an edge to any of banned_next[].  On success
 * fills parent[] and returns true.
 */
static bool restricted_dijkstra(const ListGraph *g, int src, int dst, const bool *banned,
                                int spur, const int *banned_next, size_t nb, double *dist,
                                int *parent, MuHeap *heap)
{
    const size_t n = list_graph_id_bound(g);
    mu_fill_double(dist, n, INFINITY);
    mu_fill_int(parent, n, -1);
    dist[src]  = 0.0;
    heap->size = 0;
    if (!mu_heap_push(heap, 0.0, src)) {
        return false;
    }
    while (heap->size > 0) {
        double key;
        int    u;
        mu_heap_pop(heap, &key, &u);
        if (key > dist[u]) {
            continue;
        }
        if (u == dst) {
            return true;
        }
        size_t cur = 0;
        int    v;
        double w;
        while (list_graph_next_neighbor(g, u, &cur, &v, &w)) {
            if (banned[v]) {
                continue;
            }
            if (u == spur) {
                bool skip = false;
                for (size_t i = 0; i < nb; ++i) {
                    if (banned_next[i] == v) {
                        skip = true;
                        break;
                    }
                }
                if (skip) {
                    continue;
                }
            }
            if (dist[u] + w < dist[v]) {
                dist[v]   = dist[u] + w;
                parent[v] = u;
                if (!mu_heap_push(heap, dist[v], v)) {
                    return false;
                }
            }
        }
    }
    return false;
}

static double path_cost(const ListGraph *g, const int *p, size_t len)
{
    double c = 0.0;
    for (size_t i = 0; i + 1 < len; ++i) {
        double w = 0.0;
        (void)list_graph_has_edge(g, p[i], p[i + 1], &w);
        c += w;
    }
    return c;
}

static bool same_path(const ListPath *a, const ListPath *b)
{
    return a->length == b->length &&
           memcmp(a->vertices, b->vertices, a->length * sizeof(int)) == 0;
}

typedef struct {
    ListPath *items;
    size_t      count;
    size_t      cap;
} PathList;

static bool plist_push(PathList *l, ListPath p)
{
    if (l->count == l->cap) {
        const size_t cap = l->cap > 0 ? l->cap * 2 : 8;
        ListPath *m = realloc(l->items, cap * sizeof(ListPath));
        if (m == NULL) {
            return false;
        }
        l->items = m;
        l->cap   = cap;
    }
    l->items[l->count++] = p;
    return true;
}

static bool plist_contains(const PathList *l, const ListPath *p)
{
    for (size_t i = 0; i < l->count; ++i) {
        if (same_path(&l->items[i], p)) {
            return true;
        }
    }
    return false;
}

ListStatus list_k_shortest_paths(const ListGraph *g, int src, int dst, size_t k,
                                     ListPath *paths, size_t *found)
{
    if (!mu_ok(g, src) || !mu_ok(g, dst) || k == 0 || paths == NULL || found == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    *found = 0;
    const size_t n = list_graph_id_bound(g);
    double *dist   = malloc((n > 0 ? n : 1) * sizeof(double));
    int    *parent = malloc((n > 0 ? n : 1) * sizeof(int));
    bool   *banned = calloc(n > 0 ? n : 1, sizeof(bool));
    int    *bnext  = malloc((k + 1) * sizeof(int));
    int    *buf    = malloc((n > 0 ? n : 1) * sizeof(int));
    MuHeap heap;
    PathList cand = { NULL, 0, 0 };
    ListStatus st = GRAPH_OK;
    if (dist == NULL || parent == NULL || banned == NULL || bnext == NULL || buf == NULL ||
        !mu_heap_init(&heap, 64)) {
        free(dist);
        free(parent);
        free(banned);
        free(bnext);
        free(buf);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }

    /* First path. */
    if (!restricted_dijkstra(g, src, dst, banned, -1, NULL, 0, dist, parent, &heap)) {
        st = heap.size == 0 ? GRAPH_ERR_NOT_FOUND : GRAPH_ERR_OUT_OF_MEMORY;
        goto done;
    }
    {
        const size_t len = list_path_from_parent(parent, n, src, dst, buf, n);
        ListPath p = { malloc(len * sizeof(int)), len, dist[dst] };
        if (p.vertices == NULL) {
            st = GRAPH_ERR_OUT_OF_MEMORY;
            goto done;
        }
        memcpy(p.vertices, buf, len * sizeof(int));
        paths[(*found)++] = p;
    }

    while (*found < k && st == GRAPH_OK) {
        const ListPath *prev = &paths[*found - 1];
        for (size_t i = 0; i + 1 < prev->length && st == GRAPH_OK; ++i) {
            const int spur = prev->vertices[i];
            /* Ban the root (except the spur) and every edge an accepted path takes after it. */
            for (size_t j = 0; j < i; ++j) {
                banned[prev->vertices[j]] = true;
            }
            size_t nb = 0;
            for (size_t a = 0; a < *found; ++a) {
                const ListPath *p = &paths[a];
                if (p->length > i + 1 && memcmp(p->vertices, prev->vertices, (i + 1) * sizeof(int)) == 0) {
                    bnext[nb++] = p->vertices[i + 1];
                }
            }
            const bool ok = restricted_dijkstra(g, spur, dst, banned, spur, bnext, nb, dist,
                                                parent, &heap);
            for (size_t j = 0; j < i; ++j) {
                banned[prev->vertices[j]] = false;
            }
            if (!ok) {
                if (heap.size != 0) {
                    st = GRAPH_ERR_OUT_OF_MEMORY;
                }
                continue;
            }
            const size_t slen = list_path_from_parent(parent, n, spur, dst, buf, n);
            ListPath c = { malloc((i + slen) * sizeof(int)), i + slen, 0.0 };
            if (c.vertices == NULL) {
                st = GRAPH_ERR_OUT_OF_MEMORY;
                break;
            }
            memcpy(c.vertices, prev->vertices, i * sizeof(int));
            memcpy(c.vertices + i, buf, slen * sizeof(int));
            c.cost = path_cost(g, c.vertices, c.length);
            bool dup = plist_contains(&cand, &c);
            for (size_t a = 0; a < *found && !dup; ++a) {
                dup = same_path(&paths[a], &c);
            }
            if (dup || !plist_push(&cand, c)) {
                free(c.vertices);
                if (!dup) {
                    st = GRAPH_ERR_OUT_OF_MEMORY;
                }
            }
        }
        if (st != GRAPH_OK || cand.count == 0) {
            break;
        }
        size_t best = 0;
        for (size_t a = 1; a < cand.count; ++a) {
            if (cand.items[a].cost < cand.items[best].cost) {
                best = a;
            }
        }
        paths[(*found)++]  = cand.items[best];
        cand.items[best]   = cand.items[--cand.count];
    }

done:
    for (size_t a = 0; a < cand.count; ++a) {
        free(cand.items[a].vertices);
    }
    free(cand.items);
    if (st != GRAPH_OK && st != GRAPH_ERR_NOT_FOUND) {
        for (size_t a = 0; a < *found; ++a) {
            list_path_free(&paths[a]);
        }
        *found = 0;
    }
    mu_heap_free(&heap);
    free(dist);
    free(parent);
    free(banned);
    free(bnext);
    free(buf);
    return st;
}
