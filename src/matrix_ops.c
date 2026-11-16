/*
 * matrix_ops.c — copy, subgraph, union, intersection, complement,
 * isomorphism, and the random / structured graph generators.
 * See matrix_algorithms.h.
 */

#include "matrix_util.h"

/* ------------------------------------------------------------------ */
/*  Building blocks                                                     */
/* ------------------------------------------------------------------ */

/** A graph over `n` IDs where exactly the IDs alive in `alive_in` are alive. */
static MatrixGraph *make_with_alive(size_t n, bool directed, const bool *alive)
{
    MatrixGraph *r = matrix_graph_create(n, directed);
    if (r == NULL) {
        return NULL;
    }
    for (size_t i = 0; i < n; ++i) {
        if (!alive[i] && !matrix_graph_remove_vertex(r, (int)i)) {
            matrix_graph_destroy(r);
            return NULL;
        }
    }
    return r;
}

/** Copies every edge of g into r (same IDs); edges already in r are kept. */
static bool copy_edges(const MatrixGraph *g, MatrixGraph *r)
{
    const size_t n = matrix_graph_id_bound(g);
    const bool directed = matrix_graph_is_directed(g);
    for (size_t u = 0; u < n; ++u) {
        if (!matrix_graph_is_alive(g, (int)u)) {
            continue;
        }
        size_t cur = 0;
        int    v;
        double w;
        while (matrix_graph_next_neighbor(g, (int)u, &cur, &v, &w)) {
            if (!directed && v < (int)u) {
                continue;       /* add_edge inserts the reverse itself */
            }
            if (!matrix_graph_has_edge(r, (int)u, v, NULL) &&
                !matrix_graph_add_edge(r, (int)u, v, w)) {
                return false;
            }
        }
    }
    return true;
}

MatrixGraph *matrix_graph_copy(const MatrixGraph *g)
{
    if (g == NULL) {
        return NULL;
    }
    const size_t n = matrix_graph_id_bound(g);
    bool *alive = calloc(n > 0 ? n : 1, sizeof(bool));
    if (alive == NULL) {
        return NULL;
    }
    for (size_t i = 0; i < n; ++i) {
        alive[i] = matrix_graph_is_alive(g, (int)i);
    }
    MatrixGraph *r = make_with_alive(n, matrix_graph_is_directed(g), alive);
    free(alive);
    if (r != NULL && !copy_edges(g, r)) {
        matrix_graph_destroy(r);
        return NULL;
    }
    return r;
}

MatrixStatus matrix_induced_subgraph(const MatrixGraph *g, const int *vertices, size_t count,
                                     MatrixGraph **out)
{
    if (g == NULL || out == NULL || (count > 0 && vertices == NULL)) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n = matrix_graph_id_bound(g);
    int *local = malloc((n > 0 ? n : 1) * sizeof(int));
    if (local == NULL) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    mu_fill_int(local, n, -1);
    for (size_t i = 0; i < count; ++i) {
        if (!mu_ok(g, vertices[i]) || local[vertices[i]] >= 0) {
            free(local);
            return GRAPH_ERR_INVALID_ARGUMENT;
        }
        local[vertices[i]] = (int)i;
    }
    const bool directed = matrix_graph_is_directed(g);
    MatrixGraph *r = matrix_graph_create(count, directed);
    if (r == NULL) {
        free(local);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    for (size_t i = 0; i < count; ++i) {
        size_t cur = 0;
        int    v;
        double w;
        while (matrix_graph_next_neighbor(g, vertices[i], &cur, &v, &w)) {
            const int j = local[v];
            if (j < 0 || (!directed && j < (int)i)) {
                continue;
            }
            if (!matrix_graph_add_edge(r, (int)i, j, w)) {
                matrix_graph_destroy(r);
                free(local);
                return GRAPH_ERR_OUT_OF_MEMORY;
            }
        }
    }
    free(local);
    *out = r;
    return GRAPH_OK;
}

MatrixStatus matrix_graph_union(const MatrixGraph *a, const MatrixGraph *b, MatrixGraph **out)
{
    if (a == NULL || b == NULL || out == NULL ||
        matrix_graph_is_directed(a) != matrix_graph_is_directed(b)) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t na = matrix_graph_id_bound(a), nb = matrix_graph_id_bound(b);
    const size_t n  = na > nb ? na : nb;
    bool *alive = calloc(n > 0 ? n : 1, sizeof(bool));
    if (alive == NULL) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    for (size_t i = 0; i < n; ++i) {
        alive[i] = (i < na && matrix_graph_is_alive(a, (int)i)) ||
                   (i < nb && matrix_graph_is_alive(b, (int)i));
    }
    MatrixGraph *r = make_with_alive(n, matrix_graph_is_directed(a), alive);
    free(alive);
    if (r == NULL || !copy_edges(a, r) || !copy_edges(b, r)) {
        matrix_graph_destroy(r);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    *out = r;
    return GRAPH_OK;
}

MatrixStatus matrix_graph_intersection(const MatrixGraph *a, const MatrixGraph *b, MatrixGraph **out)
{
    if (a == NULL || b == NULL || out == NULL ||
        matrix_graph_is_directed(a) != matrix_graph_is_directed(b)) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t na = matrix_graph_id_bound(a), nb = matrix_graph_id_bound(b);
    const size_t n  = na > nb ? na : nb;
    const bool directed = matrix_graph_is_directed(a);
    bool *alive = calloc(n > 0 ? n : 1, sizeof(bool));
    if (alive == NULL) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    for (size_t i = 0; i < n; ++i) {
        alive[i] = (i < na && matrix_graph_is_alive(a, (int)i)) &&
                   (i < nb && matrix_graph_is_alive(b, (int)i));
    }
    MatrixGraph *r = make_with_alive(n, directed, alive);
    free(alive);
    if (r == NULL) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    for (size_t u = 0; u < na; ++u) {
        if (!matrix_graph_is_alive(r, (int)u) || !matrix_graph_is_alive(a, (int)u)) {
            continue;
        }
        size_t cur = 0;
        int    v;
        double w;
        while (matrix_graph_next_neighbor(a, (int)u, &cur, &v, &w)) {
            if ((!directed && v < (int)u) || !matrix_graph_is_alive(r, v)) {
                continue;
            }
            if (matrix_graph_has_edge(b, (int)u, v, NULL) && !matrix_graph_add_edge(r, (int)u, v, w)) {
                matrix_graph_destroy(r);
                return GRAPH_ERR_OUT_OF_MEMORY;
            }
        }
    }
    *out = r;
    return GRAPH_OK;
}

MatrixStatus matrix_graph_complement(const MatrixGraph *g, MatrixGraph **out)
{
    if (g == NULL || out == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n = matrix_graph_id_bound(g);
    const bool directed = matrix_graph_is_directed(g);
    bool *alive = calloc(n > 0 ? n : 1, sizeof(bool));
    if (alive == NULL) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    for (size_t i = 0; i < n; ++i) {
        alive[i] = matrix_graph_is_alive(g, (int)i);
    }
    MatrixGraph *r = make_with_alive(n, directed, alive);
    free(alive);
    if (r == NULL) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    for (size_t u = 0; u < n; ++u) {
        if (!matrix_graph_is_alive(g, (int)u)) {
            continue;
        }
        for (size_t v = directed ? 0 : u + 1; v < n; ++v) {
            if (v == u || !matrix_graph_is_alive(g, (int)v) ||
                matrix_graph_has_edge(g, (int)u, (int)v, NULL)) {
                continue;
            }
            if (!matrix_graph_add_edge(r, (int)u, (int)v, 1.0)) {
                matrix_graph_destroy(r);
                return GRAPH_ERR_OUT_OF_MEMORY;
            }
        }
    }
    *out = r;
    return GRAPH_OK;
}

/* ------------------------------------------------------------------ */
/*  Isomorphism                                                         */
/* ------------------------------------------------------------------ */

typedef struct {
    const MatrixGraph *a, *b;
    bool   directed;
    size_t count;           /* vertices to place                        */
    int   *order;           /* a's vertices in placement order          */
    int   *map;             /* a vertex -> b vertex, -1 if not placed   */
    bool  *taken;           /* b vertex already used                    */
    size_t *sig_a, *sig_b;  /* (out-degree, in-degree, loop) packed     */
} IsoCtx;

static bool iso_place(IsoCtx *c, size_t pos)
{
    if (pos == c->count) {
        return true;
    }
    const int va = c->order[pos];
    const size_t nb = matrix_graph_id_bound(c->b);
    for (size_t wi = 0; wi < nb; ++wi) {
        const int w = (int)wi;
        if (c->taken[wi] || !matrix_graph_is_alive(c->b, w) || c->sig_b[wi] != c->sig_a[va]) {
            continue;
        }
        bool fits = true;
        for (size_t q = 0; q < pos && fits; ++q) {
            const int vq = c->order[q];
            const int wq = c->map[vq];
            fits = matrix_graph_has_edge(c->a, va, vq, NULL) == matrix_graph_has_edge(c->b, w, wq, NULL);
            if (fits && c->directed) {
                fits = matrix_graph_has_edge(c->a, vq, va, NULL) ==
                       matrix_graph_has_edge(c->b, wq, w, NULL);
            }
        }
        if (!fits) {
            continue;
        }
        c->map[va]   = w;
        c->taken[wi] = true;
        if (iso_place(c, pos + 1)) {
            return true;
        }
        c->taken[wi] = false;
        c->map[va]   = -1;
    }
    return false;
}

static void iso_signatures(const MatrixGraph *g, size_t *sig)
{
    const size_t n = matrix_graph_id_bound(g);
    const bool directed = matrix_graph_is_directed(g);
    size_t *indeg = calloc(n > 0 ? n : 1, sizeof(size_t));
    for (size_t u = 0; u < n; ++u) {
        sig[u] = 0;
        if (!matrix_graph_is_alive(g, (int)u)) {
            continue;
        }
        size_t cur = 0;
        int    v;
        while (indeg != NULL && directed && matrix_graph_next_neighbor(g, (int)u, &cur, &v, NULL)) {
            ++indeg[v];
        }
    }
    for (size_t u = 0; u < n; ++u) {
        if (!matrix_graph_is_alive(g, (int)u)) {
            continue;
        }
        const size_t out = matrix_graph_degree(g, (int)u);
        const size_t in  = directed && indeg != NULL ? indeg[u] : out;
        const size_t lp  = matrix_graph_has_edge(g, (int)u, (int)u, NULL) ? 1u : 0u;
        /* Pack three small counts into one comparable value. */
        sig[u] = (out * 1000003u + in) * 2u + lp;
    }
    free(indeg);
}

MatrixStatus matrix_is_isomorphic(const MatrixGraph *a, const MatrixGraph *b, bool *result,
                                  int *mapping)
{
    if (a == NULL || b == NULL || result == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    *result = false;
    if (matrix_graph_is_directed(a) != matrix_graph_is_directed(b) ||
        matrix_graph_vertex_count(a) != matrix_graph_vertex_count(b) ||
        matrix_graph_edge_count(a) != matrix_graph_edge_count(b)) {
        return GRAPH_OK;
    }
    const size_t na = matrix_graph_id_bound(a), nb = matrix_graph_id_bound(b);
    IsoCtx c;
    memset(&c, 0, sizeof(c));
    c.a = a;
    c.b = b;
    c.directed = matrix_graph_is_directed(a);
    c.count    = matrix_graph_vertex_count(a);
    c.order = malloc((c.count > 0 ? c.count : 1) * sizeof(int));
    c.map   = malloc((na > 0 ? na : 1) * sizeof(int));
    c.taken = calloc(nb > 0 ? nb : 1, sizeof(bool));
    c.sig_a = malloc((na > 0 ? na : 1) * sizeof(size_t));
    c.sig_b = malloc((nb > 0 ? nb : 1) * sizeof(size_t));
    if (c.order == NULL || c.map == NULL || c.taken == NULL || c.sig_a == NULL || c.sig_b == NULL) {
        free(c.order); free(c.map); free(c.taken); free(c.sig_a); free(c.sig_b);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    mu_fill_int(c.map, na, -1);
    iso_signatures(a, c.sig_a);
    iso_signatures(b, c.sig_b);

    /* Quick reject: the sorted signature multisets must agree. */
    size_t k = 0;
    for (size_t i = 0; i < na; ++i) {
        if (matrix_graph_is_alive(a, (int)i)) {
            c.order[k++] = (int)i;
        }
    }
    /* Place high-degree vertices first: they prune the search hardest. */
    for (size_t i = 1; i < k; ++i) {
        const int v = c.order[i];
        size_t j = i;
        while (j > 0 && c.sig_a[c.order[j - 1]] < c.sig_a[v]) {
            c.order[j] = c.order[j - 1];
            --j;
        }
        c.order[j] = v;
    }
    size_t *sb = malloc((c.count > 0 ? c.count : 1) * sizeof(size_t));
    bool same_multiset = sb != NULL;
    if (same_multiset) {
        size_t m = 0;
        for (size_t i = 0; i < nb; ++i) {
            if (matrix_graph_is_alive(b, (int)i)) {
                sb[m++] = c.sig_b[i];
            }
        }
        for (size_t i = 1; i < m; ++i) {            /* insertion sort, descending */
            const size_t v = sb[i];
            size_t j = i;
            while (j > 0 && sb[j - 1] < v) {
                sb[j] = sb[j - 1];
                --j;
            }
            sb[j] = v;
        }
        for (size_t i = 0; i < m && same_multiset; ++i) {
            same_multiset = sb[i] == c.sig_a[c.order[i]];
        }
    }
    free(sb);

    MatrixStatus st = GRAPH_OK;
    if (same_multiset) {
        *result = iso_place(&c, 0);
        if (*result && mapping != NULL) {
            memcpy(mapping, c.map, na * sizeof(int));
        }
    } else if (sb == NULL && c.count > 0) {
        st = GRAPH_ERR_OUT_OF_MEMORY;
    }
    free(c.order); free(c.map); free(c.taken); free(c.sig_a); free(c.sig_b);
    return st;
}

/* ------------------------------------------------------------------ */
/*  Generators                                                          */
/* ------------------------------------------------------------------ */

/** xorshift64*: tiny, fast, deterministic, good enough for graph generation. */
typedef struct {
    unsigned long long s;
} Rng;

static void rng_seed(Rng *r, unsigned long long seed)
{
    r->s = seed ^ 0x9E3779B97F4A7C15ULL;
    if (r->s == 0) {
        r->s = 0x2545F4914F6CDD1DULL;
    }
}

static unsigned long long rng_next(Rng *r)
{
    r->s ^= r->s >> 12;
    r->s ^= r->s << 25;
    r->s ^= r->s >> 27;
    return r->s * 0x2545F4914F6CDD1DULL;
}

/** Uniform double in [0, 1). */
static double rng_unit(Rng *r)
{
    return (double)(rng_next(r) >> 11) / 9007199254740992.0;
}

static bool weight_range_ok(double wmin, double wmax)
{
    return isfinite(wmin) && isfinite(wmax) && wmin <= wmax;
}

MatrixGraph *matrix_gen_erdos_renyi(size_t n, double p, unsigned long long seed, bool directed,
                                    double wmin, double wmax)
{
    if (!(p >= 0.0 && p <= 1.0) || !weight_range_ok(wmin, wmax)) {
        return NULL;
    }
    MatrixGraph *g = matrix_graph_create(n, directed);
    if (g == NULL) {
        return NULL;
    }
    Rng rng;
    rng_seed(&rng, seed);
    for (size_t u = 0; u < n; ++u) {
        for (size_t v = directed ? 0 : u + 1; v < n; ++v) {
            if (v == u) {
                continue;
            }
            if (rng_unit(&rng) < p) {
                const double w = wmin + rng_unit(&rng) * (wmax - wmin);
                if (!matrix_graph_add_edge(g, (int)u, (int)v, w)) {
                    matrix_graph_destroy(g);
                    return NULL;
                }
            }
        }
    }
    return g;
}

MatrixGraph *matrix_gen_barabasi_albert(size_t n, size_t m, unsigned long long seed)
{
    if (m == 0) {
        return NULL;
    }
    MatrixGraph *g = matrix_graph_create(n, false);
    if (g == NULL || n == 0) {
        return g;
    }
    const size_t seed_size = (m + 1 < n) ? m + 1 : n;
    size_t cap = seed_size * seed_size + 2 * m * n + 2;
    int *pool = malloc(cap * sizeof(int));        /* each vertex appears once per incident edge */
    int *pick = malloc((m + 1) * sizeof(int));
    if (pool == NULL || pick == NULL) {
        free(pool);
        free(pick);
        matrix_graph_destroy(g);
        return NULL;
    }
    size_t pn = 0;
    Rng rng;
    rng_seed(&rng, seed);
    bool ok = true;
    for (size_t u = 0; u < seed_size && ok; ++u) {
        for (size_t v = u + 1; v < seed_size && ok; ++v) {
            ok = matrix_graph_add_edge(g, (int)u, (int)v, 1.0);
            pool[pn++] = (int)u;
            pool[pn++] = (int)v;
        }
    }
    for (size_t v = seed_size; v < n && ok; ++v) {
        size_t chosen = 0;
        while (chosen < m) {
            const int t = pool[rng_next(&rng) % pn];
            bool dup = false;
            for (size_t i = 0; i < chosen && !dup; ++i) {
                dup = pick[i] == t;
            }
            if (!dup) {
                pick[chosen++] = t;
            }
        }
        for (size_t i = 0; i < m && ok; ++i) {
            ok = matrix_graph_add_edge(g, (int)v, pick[i], 1.0);
            pool[pn++] = (int)v;
            pool[pn++] = pick[i];
        }
    }
    free(pool);
    free(pick);
    if (!ok) {
        matrix_graph_destroy(g);
        return NULL;
    }
    return g;
}

MatrixGraph *matrix_gen_grid(size_t rows, size_t cols)
{
    if (rows != 0 && cols > SIZE_MAX / rows) {
        return NULL;
    }
    MatrixGraph *g = matrix_graph_create(rows * cols, false);
    if (g == NULL) {
        return NULL;
    }
    for (size_t r = 0; r < rows; ++r) {
        for (size_t c = 0; c < cols; ++c) {
            const int id = (int)(r * cols + c);
            if ((c + 1 < cols && !matrix_graph_add_edge(g, id, id + 1, 1.0)) ||
                (r + 1 < rows && !matrix_graph_add_edge(g, id, id + (int)cols, 1.0))) {
                matrix_graph_destroy(g);
                return NULL;
            }
        }
    }
    return g;
}

MatrixGraph *matrix_gen_path(size_t n, bool directed)
{
    MatrixGraph *g = matrix_graph_create(n, directed);
    for (size_t i = 0; g != NULL && i + 1 < n; ++i) {
        if (!matrix_graph_add_edge(g, (int)i, (int)i + 1, 1.0)) {
            matrix_graph_destroy(g);
            return NULL;
        }
    }
    return g;
}

MatrixGraph *matrix_gen_cycle(size_t n, bool directed)
{
    MatrixGraph *g = matrix_gen_path(n, directed);
    if (g == NULL) {
        return NULL;
    }
    /* Close the ring.  An undirected ring needs 3+ vertices (2 would repeat
     * the edge); a directed ring of 1 is a self-loop. */
    const bool close = directed ? n >= 1 : n >= 3;
    if (close && !matrix_graph_add_edge(g, (int)(n - 1), 0, 1.0)) {
        matrix_graph_destroy(g);
        return NULL;
    }
    return g;
}

MatrixGraph *matrix_gen_star(size_t n)
{
    MatrixGraph *g = matrix_graph_create(n, false);
    for (size_t i = 1; g != NULL && i < n; ++i) {
        if (!matrix_graph_add_edge(g, 0, (int)i, 1.0)) {
            matrix_graph_destroy(g);
            return NULL;
        }
    }
    return g;
}

MatrixGraph *matrix_gen_complete(size_t n, bool directed)
{
    MatrixGraph *g = matrix_graph_create(n, directed);
    for (size_t u = 0; g != NULL && u < n; ++u) {
        for (size_t v = directed ? 0 : u + 1; v < n; ++v) {
            if (v != u && !matrix_graph_add_edge(g, (int)u, (int)v, 1.0)) {
                matrix_graph_destroy(g);
                return NULL;
            }
        }
    }
    return g;
}

MatrixGraph *matrix_gen_random_dag(size_t n, double p, unsigned long long seed,
                                   double wmin, double wmax)
{
    if (!(p >= 0.0 && p <= 1.0) || !weight_range_ok(wmin, wmax)) {
        return NULL;
    }
    MatrixGraph *g = matrix_graph_create(n, true);
    Rng rng;
    rng_seed(&rng, seed);
    for (size_t u = 0; g != NULL && u < n; ++u) {
        for (size_t v = u + 1; v < n; ++v) {
            if (rng_unit(&rng) < p) {
                const double w = wmin + rng_unit(&rng) * (wmax - wmin);
                if (!matrix_graph_add_edge(g, (int)u, (int)v, w)) {
                    matrix_graph_destroy(g);
                    return NULL;
                }
            }
        }
    }
    return g;
}

MatrixGraph *matrix_gen_random_tree(size_t n, unsigned long long seed)
{
    MatrixGraph *g = matrix_graph_create(n, false);
    Rng rng;
    rng_seed(&rng, seed);
    for (size_t i = 1; g != NULL && i < n; ++i) {
        const int parent = (int)(rng_next(&rng) % i);
        if (!matrix_graph_add_edge(g, parent, (int)i, 1.0)) {
            matrix_graph_destroy(g);
            return NULL;
        }
    }
    return g;
}
