/*
 * list_util.h — small containers and helpers shared by the algorithm
 * modules.  Internal: not installed, not part of the public API.
 *
 * Everything here is `static inline` so each translation unit gets its own
 * copy and no extra object file is needed.
 */

#ifndef LIST_UTIL_H
#define LIST_UTIL_H

#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "adjacency_list.h"
#include "list_algorithms.h"

/* ------------------------------------------------------------------ */
/*  Allocation                                                          */
/* ------------------------------------------------------------------ */

/** calloc that never returns NULL for a zero-sized request. */
static inline void *mu_calloc(size_t n, size_t size)
{
    return calloc(n > 0 ? n : 1, size);
}

/**
 * n * size, or 0 if the product is too big to be a real allocation (so the caller can treat 0 as
 * TOO_LARGE).  The limit is PTRDIFF_MAX, not SIZE_MAX: no object can be larger than that, and an
 * optimising compiler knows it, so it rejects (-Walloc-size-larger-than) any malloc whose size it
 * cannot prove is below it.  Guarding here keeps every caller correct and the build warning-free
 * from -O0 to -O3.
 */
static inline size_t mu_mul(size_t a, size_t b)
{
    if (a != 0 && b > (size_t)PTRDIFF_MAX / a) {
        return 0;
    }
    return a * b;
}

static inline void mu_fill_int(int *a, size_t n, int value)
{
    for (size_t i = 0; i < n; ++i) {
        a[i] = value;
    }
}

static inline void mu_fill_double(double *a, size_t n, double value)
{
    for (size_t i = 0; i < n; ++i) {
        a[i] = value;
    }
}

/** Valid, alive vertex of a non-NULL graph. */
static inline bool mu_ok(const ListGraph *g, int v)
{
    return g != NULL && list_graph_is_alive(g, v);
}

/* ------------------------------------------------------------------ */
/*  Binary min-heap of (key, vertex); duplicates allowed (lazy deletion) */
/* ------------------------------------------------------------------ */

typedef struct {
    double *key;
    int    *val;
    size_t  size;
    size_t  cap;
} MuHeap;

static inline bool mu_heap_init(MuHeap *h, size_t cap)
{
    h->size = 0;
    h->cap  = cap > 8 ? cap : 8;
    h->key  = malloc(h->cap * sizeof(double));
    h->val  = malloc(h->cap * sizeof(int));
    if (h->key == NULL || h->val == NULL) {
        free(h->key);
        free(h->val);
        h->key = NULL;
        h->val = NULL;
        return false;
    }
    return true;
}

static inline void mu_heap_free(MuHeap *h)
{
    free(h->key);
    free(h->val);
    h->key = NULL;
    h->val = NULL;
}

static inline bool mu_heap_push(MuHeap *h, double key, int val)
{
    if (h->size == h->cap) {
        const size_t cap = h->cap * 2;
        double *nk = realloc(h->key, cap * sizeof(double));
        if (nk == NULL) {
            return false;
        }
        h->key = nk;
        int *nv = realloc(h->val, cap * sizeof(int));
        if (nv == NULL) {
            return false;
        }
        h->val = nv;
        h->cap = cap;
    }
    size_t i = h->size++;
    while (i > 0) {
        const size_t p = (i - 1) / 2;
        if (h->key[p] <= key) {
            break;
        }
        h->key[i] = h->key[p];
        h->val[i] = h->val[p];
        i = p;
    }
    h->key[i] = key;
    h->val[i] = val;
    return true;
}

/** Pops the smallest key; the heap must not be empty. */
static inline void mu_heap_pop(MuHeap *h, double *key, int *val)
{
    *key = h->key[0];
    *val = h->val[0];
    --h->size;
    if (h->size == 0) {
        return;
    }
    const double k = h->key[h->size];
    const int    v = h->val[h->size];
    size_t i = 0;
    for (;;) {
        size_t c = 2 * i + 1;
        if (c >= h->size) {
            break;
        }
        if (c + 1 < h->size && h->key[c + 1] < h->key[c]) {
            ++c;
        }
        if (h->key[c] >= k) {
            break;
        }
        h->key[i] = h->key[c];
        h->val[i] = h->val[c];
        i = c;
    }
    h->key[i] = k;
    h->val[i] = v;
}

/* ------------------------------------------------------------------ */
/*  Union-find                                                          */
/* ------------------------------------------------------------------ */

typedef struct {
    int *parent;
    int *rank;
} MuUF;

static inline bool mu_uf_init(MuUF *u, size_t n)
{
    u->parent = malloc((n > 0 ? n : 1) * sizeof(int));
    u->rank   = calloc(n > 0 ? n : 1, sizeof(int));
    if (u->parent == NULL || u->rank == NULL) {
        free(u->parent);
        free(u->rank);
        return false;
    }
    for (size_t i = 0; i < n; ++i) {
        u->parent[i] = (int)i;
    }
    return true;
}

static inline void mu_uf_free(MuUF *u)
{
    free(u->parent);
    free(u->rank);
}

static inline int mu_uf_find(MuUF *u, int x)
{
    while (u->parent[x] != x) {
        u->parent[x] = u->parent[u->parent[x]];
        x = u->parent[x];
    }
    return x;
}

/** Merges the sets of a and b; returns false if they were already joined. */
static inline bool mu_uf_union(MuUF *u, int a, int b)
{
    a = mu_uf_find(u, a);
    b = mu_uf_find(u, b);
    if (a == b) {
        return false;
    }
    if (u->rank[a] < u->rank[b]) {
        const int t = a;
        a = b;
        b = t;
    }
    u->parent[b] = a;
    if (u->rank[a] == u->rank[b]) {
        ++u->rank[a];
    }
    return true;
}

/* ------------------------------------------------------------------ */
/*  Bitsets (arrays of uint64_t)                                        */
/* ------------------------------------------------------------------ */

static inline size_t mu_words(size_t bits)
{
    return (bits + 63) / 64;
}

static inline void mu_bit_set(uint64_t *b, size_t i)
{
    b[i >> 6] |= (uint64_t)1 << (i & 63);
}

static inline void mu_bit_clear(uint64_t *b, size_t i)
{
    b[i >> 6] &= ~((uint64_t)1 << (i & 63));
}

static inline bool mu_bit_test(const uint64_t *b, size_t i)
{
    return (b[i >> 6] >> (i & 63)) & 1u;
}

static inline size_t mu_bit_count(const uint64_t *b, size_t words)
{
    size_t c = 0;
    for (size_t i = 0; i < words; ++i) {
        uint64_t x = b[i];
        while (x != 0) {
            x &= x - 1;
            ++c;
        }
    }
    return c;
}

/* ------------------------------------------------------------------ */
/*  Traversal helper                                                    */
/* ------------------------------------------------------------------ */

/**
 * mu_bfs_hops - hop distances from s along out-edges.
 * @dist:  id_bound entries; set to -1 for unreached and dead IDs.
 * @queue: id_bound entries of scratch.
 * Returns the number of vertices reached, including s.
 */
static inline size_t mu_bfs_hops(const ListGraph *g, int s, int *dist, int *queue)
{
    const size_t n = list_graph_id_bound(g);
    mu_fill_int(dist, n, -1);
    size_t head = 0, tail = 0;
    dist[s]        = 0;
    queue[tail++]  = s;
    while (head < tail) {
        const int u = queue[head++];
        size_t cur = 0;
        int    v;
        while (list_graph_next_neighbor(g, u, &cur, &v, NULL)) {
            if (dist[v] < 0) {
                dist[v]       = dist[u] + 1;
                queue[tail++] = v;
            }
        }
    }
    return tail;
}

/* ------------------------------------------------------------------ */
/*  Symmetrised adjacency (CSR)                                         */
/*                                                                      */
/*  The "underlying undirected graph": one entry per unordered pair      */
/*  {u, v}, u != v, in both rows.  For a directed graph with arcs both   */
/*  ways the weights add; self-loops are dropped.  Rows are sorted by    */
/*  neighbour ID, so everything built on it is deterministic.            */
/* ------------------------------------------------------------------ */

typedef struct {
    size_t  n;        /* id_bound                       */
    size_t *off;      /* n + 1 row offsets              */
    int    *to;
    double *w;
} MuSym;

static inline void mu_sym_free(MuSym *s)
{
    free(s->off);
    free(s->to);
    free(s->w);
    memset(s, 0, sizeof(*s));
}

static inline size_t mu_sym_degree(const MuSym *s, int v)
{
    return s->off[v + 1] - s->off[v];
}

/** One pass over the pairs; mode 0 counts into off[], mode 1 fills to[] / w[]. */
static inline void mu_sym_pass(const ListGraph *g, MuSym *s, size_t *fill, int mode)
{
    const bool   directed = list_graph_is_directed(g);
    const size_t n        = list_graph_id_bound(g);
    for (size_t ui = 0; ui < n; ++ui) {
        const int u = (int)ui;
        if (!list_graph_is_alive(g, u)) {
            continue;
        }
        size_t cur = 0;
        int    v;
        double w;
        while (list_graph_next_neighbor(g, u, &cur, &v, &w)) {
            if (v == u) {
                continue;
            }
            if (!directed) {
                /* The undirected row of u lists v, and v's row lists u. */
                if (mode == 0) {
                    ++s->off[u + 1];
                } else {
                    const size_t k = fill[u]++;
                    s->to[k] = v;
                    if (s->w != NULL) {
                        s->w[k] = w;
                    }
                }
                continue;
            }
            /* Directed: one entry per unordered pair.  With arcs both ways the weights add
             * and the smaller endpoint emits the pair. */
            double wsum = w;
            double w2;
            if (list_graph_has_edge(g, v, u, &w2)) {
                wsum = w + w2;
                if (v < u) {
                    continue;
                }
            }
            if (mode == 0) {
                ++s->off[u + 1];
                ++s->off[v + 1];
            } else {
                size_t k = fill[u]++;
                s->to[k] = v;
                if (s->w != NULL) {
                    s->w[k] = wsum;
                }
                k = fill[v]++;
                s->to[k] = u;
                if (s->w != NULL) {
                    s->w[k] = wsum;
                }
            }
        }
    }
}

/* Sorting one row of the symmetrised adjacency by neighbour ID.  Short rows (the common case) use insertion
 * sort, which is fastest on nearly sorted data; long rows (hubs) use qsort so a row of d entries costs
 * O(d log d), not O(d^2).  tmp is scratch for the weighted case (one row long) and may be NULL otherwise. */
typedef struct {
    int    to;
    double w;
} MuSymEntry;

static inline int mu_sym_cmp_int(const void *a, const void *b)
{
    const int x = *(const int *)a, y = *(const int *)b;
    return (x > y) - (x < y);
}

static inline int mu_sym_cmp_entry(const void *a, const void *b)
{
    const int x = ((const MuSymEntry *)a)->to, y = ((const MuSymEntry *)b)->to;
    return (x > y) - (x < y);
}

#define MU_SYM_INSERTION_MAX 32

static inline void mu_sym_sort_row(MuSym *s, size_t lo, size_t hi, MuSymEntry *tmp)
{
    const size_t len = hi - lo;
    if (len < 2) {
        return;
    }
    if (len > MU_SYM_INSERTION_MAX) {
        if (s->w == NULL) {
            qsort(s->to + lo, len, sizeof(int), mu_sym_cmp_int);
        } else if (tmp != NULL) {       /* always true here: mu_sym_build_ex sizes tmp for the longest row */
            for (size_t i = 0; i < len; ++i) {
                tmp[i].to = s->to[lo + i];
                tmp[i].w  = s->w[lo + i];
            }
            qsort(tmp, len, sizeof(*tmp), mu_sym_cmp_entry);
            for (size_t i = 0; i < len; ++i) {
                s->to[lo + i] = tmp[i].to;
                s->w[lo + i]  = tmp[i].w;
            }
        }
        return;
    }
    for (size_t i = lo + 1; i < hi; ++i) {
        const int    t = s->to[i];
        const double x = s->w != NULL ? s->w[i] : 0.0;
        size_t j = i;
        while (j > lo && s->to[j - 1] > t) {
            s->to[j] = s->to[j - 1];
            if (s->w != NULL) {
                s->w[j] = s->w[j - 1];
            }
            --j;
        }
        s->to[j] = t;
        if (s->w != NULL) {
            s->w[j] = x;
        }
    }
}

/**
 * Builds the symmetrised adjacency (rows sorted by neighbour ID).  Returns false on allocation failure.
 * With want_weights == false the weight array is not allocated and s->w is NULL: most algorithms
 * only walk off[] / to[], and w[] was 8 of the 12 bytes stored per entry.  Only read s->w if you asked for it.
 */
static inline bool mu_sym_build_ex(const ListGraph *g, MuSym *s, bool want_weights)
{
    memset(s, 0, sizeof(*s));
    s->n   = list_graph_id_bound(g);
    s->off = calloc(s->n + 2, sizeof(size_t));
    if (s->off == NULL) {
        return false;
    }
    mu_sym_pass(g, s, NULL, 0);
    size_t maxrow = 0;
    for (size_t i = 0; i < s->n; ++i) {
        if (s->off[i + 1] > maxrow) {
            maxrow = s->off[i + 1];
        }
        s->off[i + 1] += s->off[i];
    }
    const size_t total = s->off[s->n];
    s->to = malloc((total > 0 ? total : 1) * sizeof(int));
    if (want_weights) {
        s->w = malloc((total > 0 ? total : 1) * sizeof(double));
    }
    size_t *fill = malloc((s->n > 0 ? s->n : 1) * sizeof(size_t));
    MuSymEntry *tmp = want_weights && maxrow > MU_SYM_INSERTION_MAX ? malloc(maxrow * sizeof(*tmp)) : NULL;
    if (s->to == NULL || (want_weights && s->w == NULL) || fill == NULL ||
        (want_weights && maxrow > MU_SYM_INSERTION_MAX && tmp == NULL)) {
        free(fill);
        free(tmp);
        mu_sym_free(s);
        return false;
    }
    for (size_t i = 0; i < s->n; ++i) {
        fill[i] = s->off[i];
    }
    mu_sym_pass(g, s, fill, 1);
    free(fill);
    for (size_t u = 0; u < s->n; ++u) {
        mu_sym_sort_row(s, s->off[u], s->off[u + 1], tmp);
    }
    free(tmp);
    return true;
}

/** Builds the symmetrised adjacency with weights.  Returns false on allocation failure. */
static inline bool mu_sym_build(const ListGraph *g, MuSym *s)
{
    return mu_sym_build_ex(g, s, true);
}

#endif /* LIST_UTIL_H */
