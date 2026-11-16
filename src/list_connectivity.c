/*
 * list_connectivity.c — components, cut vertices, bridges, blocks, cycles,
 * transitive closure / reduction.  See list_algorithms.h.
 */

#include "list_util.h"

/* ------------------------------------------------------------------ */
/*  Simple predicates                                                   */
/* ------------------------------------------------------------------ */

bool list_is_connected(const ListGraph *g)
{
    if (g == NULL) {
        return false;
    }
    const size_t vc = list_graph_vertex_count(g);
    if (vc == 0) {
        return true;
    }
    const size_t n = list_graph_id_bound(g);
    int *comp = malloc(n * sizeof(int));
    if (comp == NULL) {
        return false;
    }
    const size_t count = list_connected_components(g, comp);
    free(comp);
    return count == 1;
}

bool list_is_strongly_connected(const ListGraph *g)
{
    if (g == NULL) {
        return false;
    }
    if (list_graph_vertex_count(g) == 0) {
        return true;
    }
    const size_t n = list_graph_id_bound(g);
    int *comp = malloc(n * sizeof(int));
    if (comp == NULL) {
        return false;
    }
    size_t count = 0;
    const ListStatus st = list_strongly_connected_components(g, comp, &count);
    free(comp);
    return st == GRAPH_OK && count == 1;
}

bool list_is_dag(const ListGraph *g)
{
    return g != NULL && list_graph_is_directed(g) && !list_has_cycle(g);
}

bool list_is_tree(const ListGraph *g)
{
    if (g == NULL || list_graph_is_directed(g)) {
        return false;
    }
    const size_t vc = list_graph_vertex_count(g);
    return vc > 0 && list_graph_edge_count(g) == vc - 1 && list_is_connected(g);
}

bool list_is_forest(const ListGraph *g)
{
    return g != NULL && !list_graph_is_directed(g) && !list_has_cycle(g);
}

/* ------------------------------------------------------------------ */
/*  Strongly connected components (Tarjan, iterative)                   */
/* ------------------------------------------------------------------ */

ListStatus list_strongly_connected_components(const ListGraph *g, int *component,
                                                  size_t *count)
{
    if (g == NULL || component == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n = list_graph_id_bound(g);
    int    *index  = malloc((n > 0 ? n : 1) * sizeof(int));
    int    *low    = malloc((n > 0 ? n : 1) * sizeof(int));
    bool   *onstk  = calloc(n > 0 ? n : 1, sizeof(bool));
    int    *stk    = malloc((n > 0 ? n : 1) * sizeof(int));
    int    *call   = malloc((n > 0 ? n : 1) * sizeof(int));
    size_t *cursor = malloc((n > 0 ? n : 1) * sizeof(size_t));
    if (index == NULL || low == NULL || onstk == NULL || stk == NULL || call == NULL ||
        cursor == NULL) {
        free(index);
        free(low);
        free(onstk);
        free(stk);
        free(call);
        free(cursor);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    mu_fill_int(index, n, -1);
    mu_fill_int(component, n, -1);

    int    counter = 0;
    size_t comps = 0, sp = 0, cp = 0;
    for (size_t si = 0; si < n; ++si) {
        if (!list_graph_is_alive(g, (int)si) || index[si] >= 0) {
            continue;
        }
        index[si] = low[si] = counter++;
        stk[sp++] = (int)si;
        onstk[si] = true;
        call[cp]   = (int)si;
        cursor[cp] = 0;
        ++cp;
        while (cp > 0) {
            const int v = call[cp - 1];
            int w;
            bool descended = false;
            while (list_graph_next_neighbor(g, v, &cursor[cp - 1], &w, NULL)) {
                if (index[w] < 0) {
                    index[w] = low[w] = counter++;
                    stk[sp++] = w;
                    onstk[w]  = true;
                    call[cp]   = w;
                    cursor[cp] = 0;
                    ++cp;
                    descended = true;
                    break;
                }
                if (onstk[w] && index[w] < low[v]) {
                    low[v] = index[w];
                }
            }
            if (descended) {
                continue;
            }
            if (low[v] == index[v]) {
                int w2;
                do {
                    w2 = stk[--sp];
                    onstk[w2]    = false;
                    component[w2] = (int)comps;
                } while (w2 != v);
                ++comps;
            }
            --cp;
            if (cp > 0) {
                const int p = call[cp - 1];
                if (low[v] < low[p]) {
                    low[p] = low[v];
                }
            }
        }
    }
    if (count != NULL) {
        *count = comps;
    }
    free(index);
    free(low);
    free(onstk);
    free(stk);
    free(call);
    free(cursor);
    return GRAPH_OK;
}

/* ------------------------------------------------------------------ */
/*  Lowlink DFS (undirected): cut vertices, bridges, blocks             */
/* ------------------------------------------------------------------ */

typedef struct {
    int    *disc;      /* discovery time, -1 if unseen */
    int    *low;
    int    *parent;
    int    *call;
    size_t *cursor;
    size_t  n;
} Lowlink;

static void lowlink_free(Lowlink *l)
{
    free(l->disc);
    free(l->low);
    free(l->parent);
    free(l->call);
    free(l->cursor);
    memset(l, 0, sizeof(*l));
}

static bool lowlink_init(Lowlink *l, size_t n)
{
    l->n      = n;
    l->disc   = malloc((n > 0 ? n : 1) * sizeof(int));
    l->low    = malloc((n > 0 ? n : 1) * sizeof(int));
    l->parent = malloc((n > 0 ? n : 1) * sizeof(int));
    l->call   = malloc((n > 0 ? n : 1) * sizeof(int));
    l->cursor = malloc((n > 0 ? n : 1) * sizeof(size_t));
    if (l->disc == NULL || l->low == NULL || l->parent == NULL || l->call == NULL ||
        l->cursor == NULL) {
        lowlink_free(l);
        return false;
    }
    mu_fill_int(l->disc, n, -1);
    mu_fill_int(l->parent, n, -1);
    return true;
}

/** Runs the lowlink DFS over every component, filling disc / low / parent. */
static void lowlink_run(const ListGraph *g, Lowlink *l)
{
    int timer = 0;
    for (size_t si = 0; si < l->n; ++si) {
        if (!list_graph_is_alive(g, (int)si) || l->disc[si] >= 0) {
            continue;
        }
        size_t cp = 0;
        l->disc[si] = l->low[si] = timer++;
        l->call[cp]   = (int)si;
        l->cursor[cp] = 0;
        ++cp;
        while (cp > 0) {
            const int v = l->call[cp - 1];
            int w;
            bool descended = false;
            while (list_graph_next_neighbor(g, v, &l->cursor[cp - 1], &w, NULL)) {
                if (w == v || w == l->parent[v]) {
                    continue;
                }
                if (l->disc[w] < 0) {
                    l->parent[w] = v;
                    l->disc[w] = l->low[w] = timer++;
                    l->call[cp]   = w;
                    l->cursor[cp] = 0;
                    ++cp;
                    descended = true;
                    break;
                }
                if (l->disc[w] < l->low[v]) {
                    l->low[v] = l->disc[w];
                }
            }
            if (descended) {
                continue;
            }
            --cp;
            if (cp > 0) {
                const int p = l->call[cp - 1];
                if (l->low[v] < l->low[p]) {
                    l->low[p] = l->low[v];
                }
            }
        }
    }
}

ListStatus list_articulation_points(const ListGraph *g, bool *is_cut, size_t *count)
{
    if (g == NULL || is_cut == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    if (list_graph_is_directed(g)) {
        return GRAPH_ERR_UNSUPPORTED;
    }
    const size_t n = list_graph_id_bound(g);
    Lowlink l;
    int *root_children = calloc(n > 0 ? n : 1, sizeof(int));
    if (root_children == NULL || !lowlink_init(&l, n)) {
        free(root_children);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    lowlink_run(g, &l);
    for (size_t i = 0; i < n; ++i) {
        is_cut[i] = false;
    }
    for (size_t c = 0; c < n; ++c) {
        const int p = l.parent[c];
        if (p < 0) {
            continue;
        }
        if (l.parent[p] < 0) {
            ++root_children[p];
        } else if (l.low[c] >= l.disc[p]) {
            is_cut[p] = true;
        }
    }
    size_t cuts = 0;
    for (size_t i = 0; i < n; ++i) {
        if (l.disc[i] >= 0 && l.parent[i] < 0 && root_children[i] >= 2) {
            is_cut[i] = true;
        }
        cuts += is_cut[i] ? 1u : 0u;
    }
    if (count != NULL) {
        *count = cuts;
    }
    lowlink_free(&l);
    free(root_children);
    return GRAPH_OK;
}

ListStatus list_bridges(const ListGraph *g, ListEdge *out, size_t *count)
{
    if (g == NULL || out == NULL || count == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    if (list_graph_is_directed(g)) {
        return GRAPH_ERR_UNSUPPORTED;
    }
    const size_t n = list_graph_id_bound(g);
    Lowlink l;
    if (!lowlink_init(&l, n)) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    lowlink_run(g, &l);
    size_t k = 0;
    for (size_t v = 0; v < n; ++v) {
        const int p = l.parent[v];
        if (p >= 0 && l.low[v] == l.disc[v]) {
            double w = 0.0;
            (void)list_graph_has_edge(g, p, (int)v, &w);
            out[k].u      = p < (int)v ? p : (int)v;
            out[k].v      = p < (int)v ? (int)v : p;
            out[k].weight = w;
            ++k;
        }
    }
    *count = k;
    lowlink_free(&l);
    return GRAPH_OK;
}

typedef struct {
    size_t *offsets;
    int    *members;
    size_t  groups;
    size_t  used;
    size_t  cap;
} GroupBuilder;

static bool gb_reserve(GroupBuilder *b, size_t more)
{
    if (b->used + more <= b->cap) {
        return true;
    }
    size_t cap = b->cap > 0 ? b->cap : 16;
    while (cap < b->used + more) {
        cap *= 2;
    }
    int *m = realloc(b->members, cap * sizeof(int));
    if (m == NULL) {
        return false;
    }
    b->members = m;
    b->cap     = cap;
    return true;
}

ListStatus list_biconnected_components(const ListGraph *g, ListGroups *out)
{
    if (g == NULL || out == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    if (list_graph_is_directed(g)) {
        return GRAPH_ERR_UNSUPPORTED;
    }
    const size_t n = list_graph_id_bound(g);
    const size_t e = list_graph_edge_count(g);
    Lowlink l;
    int *eu    = malloc((e + 1) * sizeof(int));
    int *ev    = malloc((e + 1) * sizeof(int));
    int *stamp = calloc(n > 0 ? n : 1, sizeof(int));
    size_t *offsets = malloc((e + 2) * sizeof(size_t));   /* at most e blocks */
    GroupBuilder b = { offsets, NULL, 0, 0, 0 };
    if (eu == NULL || ev == NULL || stamp == NULL || offsets == NULL || !lowlink_init(&l, n)) {
        free(eu);
        free(ev);
        free(stamp);
        free(offsets);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    offsets[0] = 0;
    size_t esp = 0;
    ListStatus st = GRAPH_OK;
    int timer = 0;

    for (size_t si = 0; si < n && st == GRAPH_OK; ++si) {
        if (!list_graph_is_alive(g, (int)si) || l.disc[si] >= 0) {
            continue;
        }
        size_t cp = 0;
        l.disc[si] = l.low[si] = timer++;
        l.call[cp]   = (int)si;
        l.cursor[cp] = 0;
        ++cp;
        while (cp > 0 && st == GRAPH_OK) {
            const int v = l.call[cp - 1];
            int w;
            bool descended = false;
            while (list_graph_next_neighbor(g, v, &l.cursor[cp - 1], &w, NULL)) {
                if (w == v || w == l.parent[v]) {
                    continue;
                }
                if (l.disc[w] < 0) {
                    eu[esp] = v;
                    ev[esp] = w;
                    ++esp;
                    l.parent[w] = v;
                    l.disc[w] = l.low[w] = timer++;
                    l.call[cp]   = w;
                    l.cursor[cp] = 0;
                    ++cp;
                    descended = true;
                    break;
                }
                if (l.disc[w] < l.disc[v]) {          /* back edge to an ancestor */
                    eu[esp] = v;
                    ev[esp] = w;
                    ++esp;
                    if (l.disc[w] < l.low[v]) {
                        l.low[v] = l.disc[w];
                    }
                }
            }
            if (descended) {
                continue;
            }
            --cp;
            if (cp == 0) {
                break;
            }
            const int p = l.call[cp - 1];
            if (l.low[v] < l.low[p]) {
                l.low[p] = l.low[v];
            }
            if (l.low[v] >= l.disc[p]) {
                /* p separates v's subtree: pop its edges into one block. */
                const int id = (int)b.groups + 1;
                const size_t start = b.used;
                for (;;) {
                    --esp;
                    const int a = eu[esp], c = ev[esp];
                    for (int x = 0; x < 2; ++x) {
                        const int vert = x == 0 ? a : c;
                        if (stamp[vert] != id) {
                            stamp[vert] = id;
                            if (!gb_reserve(&b, 1)) {
                                st = GRAPH_ERR_OUT_OF_MEMORY;
                                break;
                            }
                            b.members[b.used++] = vert;
                        }
                    }
                    if (st != GRAPH_OK || (a == p && c == v)) {
                        break;
                    }
                }
                (void)start;
                ++b.groups;
                offsets[b.groups] = b.used;
            }
        }
    }
    free(eu);
    free(ev);
    free(stamp);
    lowlink_free(&l);
    if (st != GRAPH_OK) {
        free(offsets);
        free(b.members);
        return st;
    }
    out->count   = b.groups;
    out->offsets = offsets;
    out->members = b.members != NULL ? b.members : malloc(sizeof(int));
    return GRAPH_OK;
}

/* ------------------------------------------------------------------ */
/*  Cycle extraction                                                    */
/* ------------------------------------------------------------------ */

size_t list_find_cycle(const ListGraph *g, int *cycle)
{
    if (g == NULL || cycle == NULL) {
        return 0;
    }
    const bool   directed = list_graph_is_directed(g);
    const size_t n = list_graph_id_bound(g);

    if (!directed) {
        for (size_t v = 0; v < n; ++v) {
            if (list_graph_is_alive(g, (int)v) && list_graph_has_edge(g, (int)v, (int)v, NULL)) {
                cycle[0] = (int)v;
                return 1;
            }
        }
    }

    unsigned char *state  = calloc(n > 0 ? n : 1, 1);   /* 0 new, 1 on stack, 2 done */
    int           *parent = malloc((n > 0 ? n : 1) * sizeof(int));
    size_t        *cursor = malloc((n > 0 ? n : 1) * sizeof(size_t));
    int           *call   = malloc((n > 0 ? n : 1) * sizeof(int));
    if (state == NULL || parent == NULL || cursor == NULL || call == NULL) {
        free(state);
        free(parent);
        free(cursor);
        free(call);
        return 0;
    }
    mu_fill_int(parent, n, -1);

    size_t found = 0;
    for (size_t si = 0; si < n && found == 0; ++si) {
        if (!list_graph_is_alive(g, (int)si) || state[si] != 0) {
            continue;
        }
        size_t cp = 0;
        state[si]  = 1;
        call[cp]   = (int)si;
        cursor[cp] = 0;
        ++cp;
        while (cp > 0 && found == 0) {
            const int v = call[cp - 1];
            int w;
            bool descended = false;
            while (list_graph_next_neighbor(g, v, &cursor[cp - 1], &w, NULL)) {
                if (!directed && (w == v || w == parent[v])) {
                    continue;
                }
                if (state[w] == 1) {
                    /* Back edge v -> w: the tree path w .. v plus this edge. */
                    size_t len = 0;
                    for (int x = v; x != w; x = parent[x]) {
                        cycle[len++] = x;
                    }
                    cycle[len++] = w;
                    for (size_t i = 0; i < len / 2; ++i) {
                        const int t = cycle[i];
                        cycle[i] = cycle[len - 1 - i];
                        cycle[len - 1 - i] = t;
                    }
                    found = len;
                    break;
                }
                if (state[w] == 0) {
                    state[w]  = 1;
                    parent[w] = v;
                    call[cp]   = w;
                    cursor[cp] = 0;
                    ++cp;
                    descended = true;
                    break;
                }
            }
            if (found != 0 || descended) {
                continue;
            }
            state[v] = 2;
            --cp;
        }
    }
    free(state);
    free(parent);
    free(cursor);
    free(call);
    return found;
}

/* ------------------------------------------------------------------ */
/*  Transitive closure and reduction                                    */
/* ------------------------------------------------------------------ */

/** A fresh graph over g's ID space (dead IDs removed), no edges. */
static ListGraph *empty_like(const ListGraph *g)
{
    const size_t n = list_graph_id_bound(g);
    ListGraph *r = list_graph_create(n, list_graph_is_directed(g));
    if (r == NULL) {
        return NULL;
    }
    for (size_t i = 0; i < n; ++i) {
        if (!list_graph_is_alive(g, (int)i) && !list_graph_remove_vertex(r, (int)i)) {
            list_graph_destroy(r);
            return NULL;
        }
    }
    return r;
}

ListStatus list_transitive_closure(const ListGraph *g, ListGraph **out)
{
    if (g == NULL || out == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n = list_graph_id_bound(g);
    const bool directed = list_graph_is_directed(g);
    ListGraph *r = empty_like(g);
    int *dist  = malloc((n > 0 ? n : 1) * sizeof(int));
    int *queue = malloc((n > 0 ? n : 1) * sizeof(int));
    if (r == NULL || dist == NULL || queue == NULL) {
        list_graph_destroy(r);
        free(dist);
        free(queue);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    for (size_t ui = 0; ui < n; ++ui) {
        if (!list_graph_is_alive(g, (int)ui)) {
            continue;
        }
        const size_t reached = mu_bfs_hops(g, (int)ui, dist, queue);
        for (size_t i = 1; i < reached; ++i) {
            const int v = queue[i];
            if (directed || (int)ui < v) {
                if (!list_graph_add_edge(r, (int)ui, v, 1.0)) {
                    list_graph_destroy(r);
                    free(dist);
                    free(queue);
                    return GRAPH_ERR_OUT_OF_MEMORY;
                }
            }
        }
    }
    free(dist);
    free(queue);
    *out = r;
    return GRAPH_OK;
}

ListStatus list_transitive_reduction(const ListGraph *g, ListGraph **out)
{
    if (g == NULL || out == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    if (!list_is_dag(g)) {
        return GRAPH_ERR_UNSUPPORTED;
    }
    const size_t n = list_graph_id_bound(g);
    const size_t vc = list_graph_vertex_count(g);
    const size_t words = mu_words(n);
    if (mu_mul(n, words) == 0 && n > 0) {
        return GRAPH_ERR_TOO_LARGE;
    }
    if (n > 65536) {
        return GRAPH_ERR_TOO_LARGE;
    }
    uint64_t *reach   = mu_calloc(n * words, sizeof(uint64_t));
    uint64_t *covered = mu_calloc(words, sizeof(uint64_t));
    int      *order   = malloc((vc > 0 ? vc : 1) * sizeof(int));
    ListGraph *r    = empty_like(g);
    if (reach == NULL || covered == NULL || order == NULL || r == NULL ||
        !list_topological_sort(g, order)) {
        free(reach);
        free(covered);
        free(order);
        list_graph_destroy(r);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    /* reach[u] = everything reachable from u by a path of >= 1 edge. */
    for (size_t i = vc; i-- > 0;) {
        const int u = order[i];
        uint64_t *ru = reach + (size_t)u * words;
        size_t cur = 0;
        int    v;
        while (list_graph_next_neighbor(g, u, &cur, &v, NULL)) {
            mu_bit_set(ru, (size_t)v);
            const uint64_t *rv = reach + (size_t)v * words;
            for (size_t k = 0; k < words; ++k) {
                ru[k] |= rv[k];
            }
        }
    }
    ListStatus st = GRAPH_OK;
    for (size_t ui = 0; ui < n && st == GRAPH_OK; ++ui) {
        const int u = (int)ui;
        if (!list_graph_is_alive(g, u)) {
            continue;
        }
        memset(covered, 0, words * sizeof(uint64_t));
        size_t cur = 0;
        int    v;
        while (list_graph_next_neighbor(g, u, &cur, &v, NULL)) {
            const uint64_t *rv = reach + (size_t)v * words;
            for (size_t k = 0; k < words; ++k) {
                covered[k] |= rv[k];
            }
        }
        cur = 0;
        double w;
        while (list_graph_next_neighbor(g, u, &cur, &v, &w)) {
            if (!mu_bit_test(covered, (size_t)v) && !list_graph_add_edge(r, u, v, w)) {
                st = GRAPH_ERR_OUT_OF_MEMORY;
                break;
            }
        }
    }
    free(reach);
    free(covered);
    free(order);
    if (st != GRAPH_OK) {
        list_graph_destroy(r);
        return st;
    }
    *out = r;
    return GRAPH_OK;
}
