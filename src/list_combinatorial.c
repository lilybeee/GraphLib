/*
 * list_combinatorial.c — Euler trails, Hamiltonian search, colouring,
 * cliques, covers, independent sets and DAG longest paths.
 * See list_algorithms.h.
 */

#include "list_util.h"

#define HAM_MAX_VERTICES 24u

/** Index of the lowest set bit; x must be non-zero. */
static unsigned lowest_bit(uint32_t x)
{
    unsigned i = 0;
    while ((x & 1u) == 0) {
        x >>= 1;
        ++i;
    }
    return i;
}

/* ------------------------------------------------------------------ */
/*  Euler trail (Hierholzer)                                            */
/* ------------------------------------------------------------------ */

ListStatus list_euler_trail(const ListGraph *g, bool closed, int *path, size_t *length)
{
    if (g == NULL || path == NULL || length == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n = list_graph_id_bound(g);
    const size_t e = list_graph_edge_count(g);
    const bool directed = list_graph_is_directed(g);
    *length = 0;
    if (e == 0) {
        return GRAPH_OK;
    }

    /* Adjacency as (neighbour, edge id) entries. */
    const size_t entries = directed ? e : 2 * e;
    size_t *off  = calloc(n + 2, sizeof(size_t));
    int    *eto  = malloc((entries + 1) * sizeof(int));
    int    *eid  = malloc((entries + 1) * sizeof(int));
    long   *bal  = calloc(n > 0 ? n : 1, sizeof(long));    /* out - in (directed) */
    bool   *used = calloc(e + 1, sizeof(bool));
    size_t *ptr  = malloc((n + 1) * sizeof(size_t));
    int    *stack = malloc((e + 2) * sizeof(int));
    int    *circ  = malloc((e + 2) * sizeof(int));
    if (off == NULL || eto == NULL || eid == NULL || bal == NULL || used == NULL ||
        ptr == NULL || stack == NULL || circ == NULL) {
        free(off); free(eto); free(eid); free(bal); free(used); free(ptr); free(stack); free(circ);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    for (int pass = 0; pass < 2; ++pass) {
        int next_id = 0;
        for (size_t ui = 0; ui < n; ++ui) {
            if (!list_graph_is_alive(g, (int)ui)) {
                continue;
            }
            size_t cur = 0;
            int    v;
            while (list_graph_next_neighbor(g, (int)ui, &cur, &v, NULL)) {
                if (!directed && v < (int)ui) {
                    continue;
                }
                const int id = next_id++;
                if (pass == 0) {
                    ++off[ui + 1];
                    if (directed) {
                        ++bal[ui];
                        --bal[v];
                    } else {
                        ++off[(size_t)v + 1];     /* loops count twice, as degree 2 */
                    }
                } else {
                    size_t k = ptr[ui]++;
                    eto[k] = v;
                    eid[k] = id;
                    if (!directed) {
                        k = ptr[v]++;
                        eto[k] = (int)ui;
                        eid[k] = id;
                    }
                }
            }
        }
        if (pass == 0) {
            for (size_t i = 0; i < n; ++i) {
                off[i + 1] += off[i];
            }
            for (size_t i = 0; i < n; ++i) {
                ptr[i] = off[i];
            }
        }
    }
    for (size_t i = 0; i < n; ++i) {
        ptr[i] = off[i];
    }

    /* Degree conditions decide where the trail must start. */
    int start = -1, odd = 0, plus = 0, minus = 0, plus_v = -1;
    bool feasible = true;
    for (size_t i = 0; i < n; ++i) {
        if (!list_graph_is_alive(g, (int)i)) {
            continue;
        }
        const size_t deg = off[i + 1] - off[i];
        if (start < 0 && deg > 0) {
            start = (int)i;
        }
        if (directed) {
            if (bal[i] == 1) {
                ++plus;
                plus_v = (int)i;
            } else if (bal[i] == -1) {
                ++minus;
            } else if (bal[i] != 0) {
                feasible = false;
            }
        } else if (deg % 2 == 1) {
            if (odd == 0) {
                plus_v = (int)i;
            }
            ++odd;
        }
    }
    if (directed) {
        if (!(plus == 0 && minus == 0) && !(plus == 1 && minus == 1 && !closed)) {
            feasible = false;
        }
        if (plus == 1) {
            start = plus_v;
        }
    } else {
        if (!(odd == 0 || (odd == 2 && !closed))) {
            feasible = false;
        }
        if (odd == 2) {
            start = plus_v;
        }
    }

    ListStatus st = GRAPH_ERR_NOT_FOUND;
    if (feasible && start >= 0) {
        size_t sp = 0, cn = 0;
        stack[sp++] = start;
        while (sp > 0) {
            const int v = stack[sp - 1];
            while (ptr[v] < off[v + 1] && used[eid[ptr[v]]]) {
                ++ptr[v];
            }
            if (ptr[v] < off[v + 1]) {
                const size_t k = ptr[v]++;
                used[eid[k]] = true;
                stack[sp++] = eto[k];
            } else {
                circ[cn++] = v;
                --sp;
            }
        }
        if (cn == e + 1) {
            for (size_t i = 0; i < cn; ++i) {
                path[i] = circ[cn - 1 - i];
            }
            *length = cn;
            st = GRAPH_OK;
        }
    }
    free(off); free(eto); free(eid); free(bal); free(used); free(ptr); free(stack); free(circ);
    return st;
}

/* ------------------------------------------------------------------ */
/*  Hamiltonian path / cycle                                            */
/* ------------------------------------------------------------------ */

ListStatus list_hamiltonian(const ListGraph *g, bool cycle, int *path)
{
    if (g == NULL || path == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n = list_graph_id_bound(g);
    const size_t m = list_graph_vertex_count(g);
    if (m == 0) {
        return GRAPH_ERR_NOT_FOUND;
    }
    if (m > HAM_MAX_VERTICES) {
        return GRAPH_ERR_TOO_LARGE;
    }
    const bool directed = list_graph_is_directed(g);
    int      ids[HAM_MAX_VERTICES];
    uint32_t adj[HAM_MAX_VERTICES];
    size_t k = 0;
    for (size_t i = 0; i < n; ++i) {
        if (list_graph_is_alive(g, (int)i)) {
            ids[k++] = (int)i;
        }
    }
    for (size_t i = 0; i < m; ++i) {
        adj[i] = 0;
        for (size_t j = 0; j < m; ++j) {
            if (i != j && list_graph_has_edge(g, ids[i], ids[j], NULL)) {
                adj[i] |= (uint32_t)1 << j;
            }
        }
    }
    if (m == 1) {
        if (cycle && !list_graph_has_edge(g, ids[0], ids[0], NULL)) {
            return GRAPH_ERR_NOT_FOUND;
        }
        path[0] = ids[0];
        return GRAPH_OK;
    }
    if (cycle && m < (directed ? 2u : 3u)) {
        return GRAPH_ERR_NOT_FOUND;
    }

    const size_t full = ((size_t)1 << m) - 1;
    uint32_t *dp = calloc(full + 1, sizeof(uint32_t));   /* dp[mask]: possible last vertices */
    if (dp == NULL) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    if (cycle) {
        dp[1] = 1;
    } else {
        for (size_t i = 0; i < m; ++i) {
            dp[(size_t)1 << i] = (uint32_t)1 << i;
        }
    }
    for (size_t mask = 1; mask <= full; ++mask) {
        uint32_t ends = dp[mask];
        while (ends != 0) {
            const unsigned i = lowest_bit(ends);
            ends &= ends - 1;
            uint32_t cand = adj[i] & ~(uint32_t)mask;
            while (cand != 0) {
                const unsigned j = lowest_bit(cand);
                cand &= cand - 1;
                dp[mask | ((size_t)1 << j)] |= (uint32_t)1 << j;
            }
        }
    }
    int end = -1;
    for (size_t i = 0; i < m && end < 0; ++i) {
        if ((dp[full] >> i) & 1u) {
            if (!cycle || (adj[i] & 1u)) {
                end = (int)i;
            }
        }
    }
    if (end < 0) {
        free(dp);
        return GRAPH_ERR_NOT_FOUND;
    }
    size_t mask = full;
    int    cur  = end;
    for (size_t pos = m; pos-- > 0;) {
        path[pos] = ids[cur];
        if (pos == 0) {
            break;
        }
        const size_t prev_mask = mask & ~((size_t)1 << cur);
        int prev = -1;
        for (size_t p = 0; p < m && prev < 0; ++p) {
            if (((dp[prev_mask] >> p) & 1u) && ((adj[p] >> cur) & 1u)) {
                prev = (int)p;
            }
        }
        mask = prev_mask;
        cur  = prev;
    }
    free(dp);
    return GRAPH_OK;
}

/* ------------------------------------------------------------------ */
/*  Colouring                                                           */
/* ------------------------------------------------------------------ */

static size_t finish_colors(const ListGraph *g, const int *color, size_t n, size_t *ncolors)
{
    int maxc = -1;
    for (size_t i = 0; i < n; ++i) {
        if (list_graph_is_alive(g, (int)i) && color[i] > maxc) {
            maxc = color[i];
        }
    }
    const size_t used = (size_t)(maxc + 1);
    if (ncolors != NULL) {
        *ncolors = used;
    }
    return used;
}

/** Lowest colour not used by v's coloured neighbours; scratch needs maxdeg + 2 slots. */
/* ------------------------------------------------------------------ */
/*  A small lazy-deletion heap, used to pick "the best remaining vertex" */
/*  in O(log V) instead of rescanning every vertex (which made DSatur   */
/*  and the greedy independent set O(V^2)).                             */
/* ------------------------------------------------------------------ */

typedef struct {
    long k1, k2;    /* priority keys; their meaning depends on `before` */
    int  v;
} HItem;

typedef bool (*HBefore)(const HItem *a, const HItem *b);

typedef struct {
    HItem  *a;
    size_t  n, cap;
    HBefore before;
} HHeap;

static bool hheap_init(HHeap *h, size_t cap, HBefore before)
{
    h->cap    = cap > 16 ? cap : 16;
    h->a      = malloc(h->cap * sizeof(HItem));
    h->n      = 0;
    h->before = before;
    return h->a != NULL;
}

static bool hheap_push(HHeap *h, HItem it)
{
    if (h->n == h->cap) {
        HItem *bigger = realloc(h->a, h->cap * 2 * sizeof(HItem));
        if (bigger == NULL) {
            return false;
        }
        h->a = bigger;
        h->cap *= 2;
    }
    size_t i = h->n++;
    while (i > 0 && h->before(&it, &h->a[(i - 1) / 2])) {
        h->a[i] = h->a[(i - 1) / 2];
        i = (i - 1) / 2;
    }
    h->a[i] = it;
    return true;
}

static HItem hheap_pop(HHeap *h)
{
    const HItem top  = h->a[0];
    const HItem last = h->a[--h->n];
    size_t i = 0;
    while (2 * i + 1 < h->n) {
        size_t c = 2 * i + 1;
        if (c + 1 < h->n && h->before(&h->a[c + 1], &h->a[c])) {
            ++c;
        }
        if (!h->before(&h->a[c], &last)) {
            break;
        }
        h->a[i] = h->a[c];
        i = c;
    }
    if (h->n > 0) {
        h->a[i] = last;
    }
    return top;
}

/* DSatur: most distinct neighbouring colours first, then most neighbours, then lowest ID. */
static bool dsatur_before(const HItem *a, const HItem *b)
{
    if (a->k1 != b->k1) return a->k1 > b->k1;
    if (a->k2 != b->k2) return a->k2 > b->k2;
    return a->v < b->v;
}

/* Independent set: fewest remaining neighbours first, then lowest ID. */
static bool lowdeg_before(const HItem *a, const HItem *b)
{
    if (a->k1 != b->k1) return a->k1 < b->k1;
    return a->v < b->v;
}

static int first_free_color(const MuSym *s, int v, const int *color, int *scratch, int stamp)
{
    const size_t deg = mu_sym_degree(s, v);
    for (size_t k = s->off[v]; k < s->off[v + 1]; ++k) {
        const int c = color[s->to[k]];
        if (c >= 0 && (size_t)c <= deg) {
            scratch[c] = stamp;
        }
    }
    int c = 0;
    while (scratch[c] == stamp) {
        ++c;
    }
    return c;
}

static ListStatus color_in_order(const ListGraph *g, const int *order, size_t count,
                                   int *color, size_t *ncolors)
{
    const size_t n = list_graph_id_bound(g);
    MuSym sym;
    if (!mu_sym_build_ex(g, &sym, false)) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    size_t maxdeg = 0;
    for (size_t v = 0; v < n; ++v) {
        if (mu_sym_degree(&sym, (int)v) > maxdeg) {
            maxdeg = mu_sym_degree(&sym, (int)v);
        }
    }
    int *scratch = calloc(maxdeg + 3, sizeof(int));
    if (scratch == NULL) {
        mu_sym_free(&sym);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    mu_fill_int(color, n, -1);
    int stamp = 0;
    for (size_t i = 0; i < count; ++i) {
        ++stamp;
        color[order[i]] = first_free_color(&sym, order[i], color, scratch, stamp);
    }
    free(scratch);
    mu_sym_free(&sym);
    (void)finish_colors(g, color, n, ncolors);
    return GRAPH_OK;
}

ListStatus list_color_greedy(const ListGraph *g, int *color, size_t *ncolors)
{
    if (g == NULL || color == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n = list_graph_id_bound(g);
    int *order = malloc((n > 0 ? n : 1) * sizeof(int));
    if (order == NULL) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    size_t k = 0;
    for (size_t i = 0; i < n; ++i) {
        if (list_graph_is_alive(g, (int)i)) {
            order[k++] = (int)i;
        }
    }
    const ListStatus st = color_in_order(g, order, k, color, ncolors);
    free(order);
    return st;
}

ListStatus list_color_welsh_powell(const ListGraph *g, int *color, size_t *ncolors)
{
    if (g == NULL || color == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n = list_graph_id_bound(g);
    MuSym sym;
    if (!mu_sym_build_ex(g, &sym, false)) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    int *order = malloc((n > 0 ? n : 1) * sizeof(int));
    if (order == NULL) {
        mu_sym_free(&sym);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    size_t k = 0;
    for (size_t i = 0; i < n; ++i) {
        if (list_graph_is_alive(g, (int)i)) {
            order[k++] = (int)i;
        }
    }
    /* Counting sort by degree descending, ID ascending: O(V + maxdegree) instead of the
     * O(V^2) insertion sort this replaces. */
    {
        size_t maxdeg = 0;
        for (size_t i = 0; i < k; ++i) {
            if (mu_sym_degree(&sym, order[i]) > maxdeg) {
                maxdeg = mu_sym_degree(&sym, order[i]);
            }
        }
        size_t *start = calloc(maxdeg + 2, sizeof(size_t));
        int    *sorted = malloc((k > 0 ? k : 1) * sizeof(int));
        if (start == NULL || sorted == NULL) {
            free(start); free(sorted); free(order);
            mu_sym_free(&sym);
            return GRAPH_ERR_OUT_OF_MEMORY;
        }
        for (size_t i = 0; i < k; ++i) {
            ++start[maxdeg - mu_sym_degree(&sym, order[i]) + 1];   /* bucket 0 = highest degree */
        }
        for (size_t d = 1; d <= maxdeg + 1; ++d) {
            start[d] += start[d - 1];
        }
        for (size_t i = 0; i < k; ++i) {   /* order[] is in ascending ID, so ties stay in ID order */
            sorted[start[maxdeg - mu_sym_degree(&sym, order[i])]++] = order[i];
        }
        memcpy(order, sorted, k * sizeof(int));
        free(start);
        free(sorted);
    }
    mu_fill_int(color, n, -1);
    size_t colored = 0;
    int c = 0;
    while (colored < k) {
        for (size_t i = 0; i < k; ++i) {
            const int v = order[i];
            if (color[v] >= 0) {
                continue;
            }
            bool clash = false;
            for (size_t e = sym.off[v]; e < sym.off[v + 1] && !clash; ++e) {
                clash = color[sym.to[e]] == c;
            }
            if (!clash) {
                color[v] = c;
                ++colored;
            }
        }
        ++c;
    }
    free(order);
    mu_sym_free(&sym);
    (void)finish_colors(g, color, n, ncolors);
    return GRAPH_OK;
}

ListStatus list_color_dsatur(const ListGraph *g, int *color, size_t *ncolors)
{
    if (g == NULL || color == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n = list_graph_id_bound(g);
    MuSym sym;
    if (!mu_sym_build_ex(g, &sym, false)) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    size_t maxdeg = 0;
    for (size_t v = 0; v < n; ++v) {
        if (mu_sym_degree(&sym, (int)v) > maxdeg) {
            maxdeg = mu_sym_degree(&sym, (int)v);
        }
    }
    const size_t words = mu_words(maxdeg + 2);
    uint64_t *seen = mu_calloc(n * words, sizeof(uint64_t));   /* colours seen around each vertex */
    int      *sat  = calloc(n > 0 ? n : 1, sizeof(int));
    int      *scratch = calloc(maxdeg + 3, sizeof(int));
    if (seen == NULL || sat == NULL || scratch == NULL) {
        free(seen); free(sat); free(scratch);
        mu_sym_free(&sym);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    mu_fill_int(color, n, -1);
    const size_t alive = list_graph_vertex_count(g);
    HHeap heap;
    if (!hheap_init(&heap, alive + 16, dsatur_before)) {
        free(seen); free(sat); free(scratch);
        mu_sym_free(&sym);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    bool heap_ok = true;
    for (size_t v = 0; v < n && heap_ok; ++v) {
        if (list_graph_is_alive(g, (int)v)) {
            heap_ok = hheap_push(&heap, (HItem){0, (long)mu_sym_degree(&sym, (int)v), (int)v});
        }
    }
    for (size_t step = 0; step < alive && heap_ok; ++step) {
        /* Pop until the entry is current: still uncoloured and carrying its present saturation. */
        int best = -1;
        while (heap.n > 0) {
            const HItem it = hheap_pop(&heap);
            if (color[it.v] < 0 && it.k1 == sat[it.v]) {
                best = it.v;
                break;
            }
        }
        const int c = first_free_color(&sym, best, color, scratch, (int)step + 1);
        color[best] = c;
        for (size_t e = sym.off[best]; e < sym.off[best + 1]; ++e) {
            const int u = sym.to[e];
            if (color[u] < 0 && (size_t)c < words * 64 &&
                !mu_bit_test(seen + (size_t)u * words, (size_t)c)) {
                mu_bit_set(seen + (size_t)u * words, (size_t)c);
                ++sat[u];
                heap_ok = hheap_push(&heap, (HItem){sat[u], (long)mu_sym_degree(&sym, u), u}) && heap_ok;
            }
        }
    }
    free(heap.a);
    if (!heap_ok) {
        free(seen); free(sat); free(scratch);
        mu_sym_free(&sym);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    free(seen); free(sat); free(scratch);
    mu_sym_free(&sym);
    (void)finish_colors(g, color, n, ncolors);
    return GRAPH_OK;
}

/* ------------------------------------------------------------------ */
/*  Cliques (Bron-Kerbosch with pivoting)                               */
/* ------------------------------------------------------------------ */

typedef struct {
    size_t    m;          /* local vertex count                    */
    size_t    words;
    uint64_t *nbr;        /* m rows of `words` words               */
    const int *ids;       /* local index -> vertex ID              */
    int       prefix;     /* vertex ID that belongs to every clique found, or -1 */
    size_t   *offsets;
    int      *members;
    size_t    groups;
    size_t    used;
    size_t    cap_groups;
    size_t    cap_members;
    bool      failed;
} CliqueCtx;

static void clique_report(CliqueCtx *c, const int *r, size_t len)
{
    if (c->groups + 1 >= c->cap_groups) {
        const size_t cap = c->cap_groups > 0 ? c->cap_groups * 2 : 16;
        size_t *o = realloc(c->offsets, (cap + 1) * sizeof(size_t));
        if (o == NULL) {
            c->failed = true;
            return;
        }
        c->offsets = o;
        c->cap_groups = cap;
    }
    const size_t extra = c->prefix >= 0 ? 1 : 0;
    if (c->used + len + extra > c->cap_members) {
        size_t cap = c->cap_members > 0 ? c->cap_members * 2 : 64;
        while (cap < c->used + len + extra) {
            cap *= 2;
        }
        int *mem = realloc(c->members, cap * sizeof(int));
        if (mem == NULL) {
            c->failed = true;
            return;
        }
        c->members = mem;
        c->cap_members = cap;
    }
    if (c->prefix >= 0) {
        c->members[c->used++] = c->prefix;
    }
    for (size_t i = 0; i < len; ++i) {
        c->members[c->used++] = c->ids[r[i]];
    }
    c->offsets[++c->groups] = c->used;
}

static void bron_kerbosch(CliqueCtx *c, int *r, size_t rlen, uint64_t *p, uint64_t *x)
{
    const size_t w = c->words;
    if (c->failed) {
        return;
    }
    bool p_empty = true, x_empty = true;
    for (size_t i = 0; i < w; ++i) {
        p_empty = p_empty && p[i] == 0;
        x_empty = x_empty && x[i] == 0;
    }
    if (p_empty) {
        if (x_empty) {
            clique_report(c, r, rlen);
        }
        return;
    }
    /* Pivot: the vertex of P ∪ X with the most neighbours inside P. */
    size_t pivot = 0, best = 0;
    bool   have = false;
    for (size_t i = 0; i < c->m; ++i) {
        if (!mu_bit_test(p, i) && !mu_bit_test(x, i)) {
            continue;
        }
        size_t cnt = 0;
        const uint64_t *row = c->nbr + i * w;
        for (size_t k = 0; k < w; ++k) {
            uint64_t t = p[k] & row[k];
            while (t != 0) {
                t &= t - 1;
                ++cnt;
            }
        }
        if (!have || cnt > best) {
            have  = true;
            best  = cnt;
            pivot = i;
        }
    }
    uint64_t *np = malloc(w * sizeof(uint64_t));
    uint64_t *nx = malloc(w * sizeof(uint64_t));
    if (np == NULL || nx == NULL) {
        free(np);
        free(nx);
        c->failed = true;
        return;
    }
    const uint64_t *prow = c->nbr + pivot * w;
    for (size_t v = 0; v < c->m; ++v) {
        if (!mu_bit_test(p, v) || mu_bit_test(prow, v)) {
            continue;
        }
        const uint64_t *row = c->nbr + v * w;
        for (size_t k = 0; k < w; ++k) {
            np[k] = p[k] & row[k];
            nx[k] = x[k] & row[k];
        }
        r[rlen] = (int)v;
        bron_kerbosch(c, r, rlen + 1, np, nx);
        mu_bit_clear(p, v);
        mu_bit_set(x, v);
    }
    free(np);
    free(nx);
}

/*
 * Maximal cliques, Eppstein-Loeffler-Strash style: take the vertices in
 * degeneracy order; for each vertex v, the cliques whose EARLIEST vertex (in that
 * order) is v live entirely inside N(v), so Bron-Kerbosch only has to run on the
 * small subgraph induced by v's neighbours (later neighbours = candidates P,
 * earlier ones = already-handled X).
 *
 * Why not one big BK over the whole graph: that needs a dense V x V bit matrix --
 * V^2/8 bytes, which is 125 GB for a million vertices.  Here the only dense
 * structure is one row-set of d x d bits at a time (d = the vertex's degree), so
 * memory is O(V + E + d^2/64) and the running time O(d * V * 3^(d/3)) with d the
 * degeneracy -- near-linear on sparse graphs.
 */
ListStatus list_maximal_cliques(const ListGraph *g, ListGroups *out)
{
    if (g == NULL || out == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    if (list_graph_is_directed(g)) {
        return GRAPH_ERR_UNSUPPORTED;
    }
    const size_t n = list_graph_id_bound(g);
    MuSym sym;
    if (!mu_sym_build_ex(g, &sym, false)) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }

    /* Degeneracy order by bucket peeling: repeatedly remove a vertex of least remaining degree. */
    size_t maxdeg = 0;
    int *deg   = malloc((n > 0 ? n : 1) * sizeof(int));
    int *pos   = malloc((n > 0 ? n : 1) * sizeof(int));      /* position in the order, -1 if dead */
    int *order = malloc((n > 0 ? n : 1) * sizeof(int));
    int *loc   = malloc((n > 0 ? n : 1) * sizeof(int));      /* local index while a vertex is processed */
    if (deg == NULL || pos == NULL || order == NULL || loc == NULL) {
        free(deg); free(pos); free(order); free(loc);
        mu_sym_free(&sym);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    size_t alive = 0;
    for (size_t v = 0; v < n; ++v) {
        loc[v] = -1;
        pos[v] = -1;
        deg[v] = (int)mu_sym_degree(&sym, (int)v);
        if (list_graph_is_alive(g, (int)v)) {
            ++alive;
            if ((size_t)deg[v] > maxdeg) {
                maxdeg = (size_t)deg[v];
            }
        }
    }
    {
        size_t *bstart = calloc(maxdeg + 2, sizeof(size_t));
        int    *bucketed = malloc((alive > 0 ? alive : 1) * sizeof(int));
        int    *where    = malloc((n > 0 ? n : 1) * sizeof(int));   /* index of v inside bucketed[] */
        if (bstart == NULL || bucketed == NULL || where == NULL) {
            free(bstart); free(bucketed); free(where);
            free(deg); free(pos); free(order); free(loc);
            mu_sym_free(&sym);
            return GRAPH_ERR_OUT_OF_MEMORY;
        }
        for (size_t v = 0; v < n; ++v) {
            if (list_graph_is_alive(g, (int)v)) {
                ++bstart[deg[v] + 1];
            }
        }
        for (size_t d = 1; d <= maxdeg + 1; ++d) {
            bstart[d] += bstart[d - 1];
        }
        {
            size_t *fill = malloc((maxdeg + 2) * sizeof(size_t));
            if (fill == NULL) {
                free(bstart); free(bucketed); free(where);
                free(deg); free(pos); free(order); free(loc);
                mu_sym_free(&sym);
                return GRAPH_ERR_OUT_OF_MEMORY;
            }
            memcpy(fill, bstart, (maxdeg + 2) * sizeof(size_t));
            for (size_t v = 0; v < n; ++v) {
                if (list_graph_is_alive(g, (int)v)) {
                    where[v] = (int)fill[deg[v]];
                    bucketed[fill[deg[v]]++] = (int)v;
                }
            }
            free(fill);
        }
        for (size_t i = 0; i < alive; ++i) {
            const int v = bucketed[i];      /* the least-degree vertex still present */
            order[i] = v;
            pos[v]   = (int)i;
            for (size_t e = sym.off[v]; e < sym.off[v + 1]; ++e) {
                const int u = sym.to[e];
                if (pos[u] >= 0 || deg[u] <= deg[v]) {
                    continue;                // already peeled, or can't drop below the current level
                }
                /* move u one bucket down: swap it with the first vertex of its bucket */
                const size_t du = (size_t)deg[u];
                const size_t pu = (size_t)where[u];
                size_t first = bstart[du];
                if (first <= i) {
                    first = i + 1;          // the processed prefix is not part of any bucket
                }
                const int w = bucketed[first];
                bucketed[pu]    = w;
                where[w]        = (int)pu;
                bucketed[first] = u;
                where[u]        = (int)first;
                bstart[du]      = first + 1;
                --deg[u];
            }
        }
        free(bstart); free(bucketed); free(where);
    }

    CliqueCtx c;
    memset(&c, 0, sizeof(c));
    c.prefix  = -1;
    c.offsets = malloc(17 * sizeof(size_t));
    c.cap_groups = 16;
    int *ids = malloc((maxdeg + 1) * sizeof(int));
    int *r   = malloc((maxdeg + 2) * sizeof(int));
    if (c.offsets == NULL || ids == NULL || r == NULL) {
        free(c.offsets); free(ids); free(r);
        free(deg); free(pos); free(order); free(loc);
        mu_sym_free(&sym);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    c.offsets[0] = 0;

    for (size_t i = 0; i < alive && !c.failed; ++i) {
        const int v = order[i];
        /* Local universe: v's neighbours, later ones first (they become P), earlier ones after. */
        size_t m = 0, np_count = 0;
        for (size_t e = sym.off[v]; e < sym.off[v + 1]; ++e) {
            if (pos[sym.to[e]] > (int)i) {
                loc[sym.to[e]] = (int)m;
                ids[m++]       = sym.to[e];
                ++np_count;
            }
        }
        for (size_t e = sym.off[v]; e < sym.off[v + 1]; ++e) {
            if (pos[sym.to[e]] < (int)i) {
                loc[sym.to[e]] = (int)m;
                ids[m++]       = sym.to[e];
            }
        }
        if (m == 0) {                       /* isolated vertex: the clique {v} */
            c.prefix = -1;
            c.ids    = ids;
            ids[0]   = v;
            r[0]     = 0;
            clique_report(&c, r, 1);
            continue;
        }
        c.m      = m;
        c.words  = mu_words(m);
        c.ids    = ids;
        c.prefix = v;
        c.nbr    = mu_calloc(m * c.words, sizeof(uint64_t));
        uint64_t *p = mu_calloc(c.words, sizeof(uint64_t));
        uint64_t *x = mu_calloc(c.words, sizeof(uint64_t));
        if (c.nbr == NULL || p == NULL || x == NULL) {
            free(c.nbr); free(p); free(x);
            c.failed = true;
            break;
        }
        for (size_t a = 0; a < m; ++a) {
            const int ua = ids[a];
            for (size_t e = sym.off[ua]; e < sym.off[ua + 1]; ++e) {
                const int b = loc[sym.to[e]];
                if (b >= 0 && (size_t)b != a) {
                    mu_bit_set(c.nbr + a * c.words, (size_t)b);
                }
            }
            if (a < np_count) {
                mu_bit_set(p, a);
            } else {
                mu_bit_set(x, a);
            }
        }
        bron_kerbosch(&c, r, 0, p, x);
        free(c.nbr);
        free(p);
        free(x);
        c.nbr = NULL;
        for (size_t a = 0; a < m; ++a) {
            loc[ids[a]] = -1;
        }
    }
    free(ids); free(r); free(deg); free(pos); free(order); free(loc);
    mu_sym_free(&sym);
    if (c.failed) {
        free(c.offsets);
        free(c.members);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    out->count   = c.groups;
    out->offsets = c.offsets;
    out->members = c.members != NULL ? c.members : malloc(sizeof(int));
    return GRAPH_OK;
}

ListStatus list_max_clique(const ListGraph *g, int *clique, size_t *size)
{
    if (g == NULL || clique == NULL || size == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    ListGroups groups;
    const ListStatus st = list_maximal_cliques(g, &groups);
    if (st != GRAPH_OK) {
        return st;
    }
    size_t best = 0, best_len = 0;
    for (size_t i = 0; i < groups.count; ++i) {
        const size_t len = groups.offsets[i + 1] - groups.offsets[i];
        if (len > best_len) {
            best_len = len;
            best = i;
        }
    }
    for (size_t i = 0; i < best_len; ++i) {
        clique[i] = groups.members[groups.offsets[best] + i];
    }
    *size = best_len;
    list_groups_free(&groups);
    return GRAPH_OK;
}

/* ------------------------------------------------------------------ */
/*  Vertex cover and independent set                                    */
/* ------------------------------------------------------------------ */

ListStatus list_vertex_cover_approx(const ListGraph *g, bool *in_cover, size_t *size)
{
    if (g == NULL || in_cover == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n = list_graph_id_bound(g);
    MuSym sym;
    if (!mu_sym_build_ex(g, &sym, false)) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    for (size_t i = 0; i < n; ++i) {
        in_cover[i] = false;
    }
    size_t count = 0;
    for (size_t u = 0; u < n; ++u) {
        if (!list_graph_is_alive(g, (int)u)) {
            continue;
        }
        if (list_graph_has_edge(g, (int)u, (int)u, NULL) && !in_cover[u]) {
            in_cover[u] = true;     /* a loop can only be covered by its own vertex */
            ++count;
        }
        for (size_t k = sym.off[u]; k < sym.off[u + 1]; ++k) {
            const int v = sym.to[k];
            if (!in_cover[u] && !in_cover[v]) {
                in_cover[u] = true;
                in_cover[v] = true;
                count += 2;
            }
        }
    }
    if (size != NULL) {
        *size = count;
    }
    mu_sym_free(&sym);
    return GRAPH_OK;
}

ListStatus list_independent_set_greedy(const ListGraph *g, bool *in_set, size_t *size)
{
    if (g == NULL || in_set == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n = list_graph_id_bound(g);
    MuSym sym;
    if (!mu_sym_build_ex(g, &sym, false)) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    int  *deg     = malloc((n > 0 ? n : 1) * sizeof(int));
    bool *removed = calloc(n > 0 ? n : 1, sizeof(bool));
    if (deg == NULL || removed == NULL) {
        free(deg);
        free(removed);
        mu_sym_free(&sym);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    for (size_t v = 0; v < n; ++v) {
        in_set[v] = false;
        deg[v]    = (int)mu_sym_degree(&sym, (int)v);
        /* A vertex with a self-loop is adjacent to itself and can never be chosen. */
        removed[v] = !list_graph_is_alive(g, (int)v) ||
                     list_graph_has_edge(g, (int)v, (int)v, NULL);
    }
    size_t count = 0;
    HHeap  heap;
    if (!hheap_init(&heap, n + 16, lowdeg_before)) {
        free(deg);
        free(removed);
        mu_sym_free(&sym);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    bool heap_ok = true;
    for (size_t v = 0; v < n && heap_ok; ++v) {
        if (!removed[v]) {
            heap_ok = hheap_push(&heap, (HItem){deg[v], 0, (int)v});
        }
    }
    while (heap_ok && heap.n > 0) {
        /* The cheapest remaining vertex; entries made stale by a later degree drop are skipped. */
        const HItem it = hheap_pop(&heap);
        if (removed[it.v] || it.k1 != deg[it.v]) {
            continue;
        }
        const int best = it.v;
        in_set[best] = true;
        ++count;
        removed[best] = true;
        for (size_t k = sym.off[best]; k < sym.off[best + 1]; ++k) {
            const int u = sym.to[k];
            if (removed[u]) {
                continue;
            }
            removed[u] = true;
            for (size_t j = sym.off[u]; j < sym.off[u + 1]; ++j) {
                const int w = sym.to[j];
                --deg[w];
                if (!removed[w]) {
                    heap_ok = hheap_push(&heap, (HItem){deg[w], 0, w}) && heap_ok;
                }
            }
        }
    }
    free(heap.a);
    if (!heap_ok) {
        free(deg);
        free(removed);
        mu_sym_free(&sym);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    if (size != NULL) {
        *size = count;
    }
    free(deg);
    free(removed);
    mu_sym_free(&sym);
    return GRAPH_OK;
}

/* ------------------------------------------------------------------ */
/*  DAG longest path                                                    */
/* ------------------------------------------------------------------ */

ListStatus list_dag_longest_path(const ListGraph *g, double *distance, int *parent, int *end)
{
    if (g == NULL || distance == NULL || parent == NULL || end == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    if (!list_graph_is_directed(g)) {
        return GRAPH_ERR_UNSUPPORTED;
    }
    const size_t n  = list_graph_id_bound(g);
    const size_t vc = list_graph_vertex_count(g);
    int *order = malloc((vc > 0 ? vc : 1) * sizeof(int));
    if (order == NULL) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    if (!list_topological_sort(g, order)) {
        free(order);
        return GRAPH_ERR_UNSUPPORTED;
    }
    for (size_t i = 0; i < n; ++i) {
        distance[i] = list_graph_is_alive(g, (int)i) ? 0.0 : -INFINITY;
        parent[i]   = -1;
    }
    for (size_t i = 0; i < vc; ++i) {
        const int u = order[i];
        size_t cur = 0;
        int    v;
        double w;
        while (list_graph_next_neighbor(g, u, &cur, &v, &w)) {
            if (distance[u] + w > distance[v]) {
                distance[v] = distance[u] + w;
                parent[v]   = u;
            }
        }
    }
    *end = -1;
    for (size_t i = 0; i < n; ++i) {
        if (list_graph_is_alive(g, (int)i) && (*end < 0 || distance[i] > distance[*end])) {
            *end = (int)i;
        }
    }
    free(order);
    return GRAPH_OK;
}

ListStatus list_critical_path(const ListGraph *g, int *path, size_t *length, double *cost)
{
    if (g == NULL || path == NULL || length == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n = list_graph_id_bound(g);
    double *dist   = malloc((n > 0 ? n : 1) * sizeof(double));
    int    *parent = malloc((n > 0 ? n : 1) * sizeof(int));
    if (dist == NULL || parent == NULL) {
        free(dist);
        free(parent);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    int end = -1;
    ListStatus st = list_dag_longest_path(g, dist, parent, &end);
    if (st == GRAPH_OK) {
        if (end < 0) {
            *length = 0;
            if (cost != NULL) {
                *cost = 0.0;
            }
        } else {
            size_t len = 0;
            for (int v = end; v >= 0; v = parent[v]) {
                ++len;
            }
            size_t i = len;
            for (int v = end; v >= 0; v = parent[v]) {
                path[--i] = v;
            }
            *length = len;
            if (cost != NULL) {
                *cost = dist[end];
            }
        }
    }
    free(dist);
    free(parent);
    return st;
}
