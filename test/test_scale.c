/*
 * test_scale.c — every algorithm that makes sense on a large graph, run at
 * 1 000, 10 000, 100 000 and 1 000 000 vertices and checked against an
 * independent answer.
 *
 * How the answers are obtained (never from the library under test):
 *   - oracle  : a straightforward implementation in this file working on a CSR
 *               copy of the graph (Dijkstra, union-find, Kosaraju, DP over a
 *               topological order, ...);
 *   - closed form : graphs whose answer is known (path, cycle, grid, star,
 *               chain / ring of K5 cliques, a flow network built from its
 *               capacities, ...);
 *   - invariant : a property any correct answer has (a colouring has no
 *               monochromatic edge, a topological order respects every arc, ...).
 *
 * Which sizes run: the list representation runs every size up to max_n; the
 * matrix representation is a bit matrix whose memory grows with V^2 (100 000
 * vertices already need about 2 GB, 1 000 000 would need hundreds), so it runs
 * only up to matrix_max.  Each algorithm also has its own size cap, because
 * some are quadratic or worse (all-pairs, closeness, closure, exact search)
 * and would never finish at a million vertices; those run at the sizes where
 * they are meant to be used.
 *
 * Usage:  scale_test [max_n] [matrix_max] [filter]
 *   max_n       largest size to run            (default 1000000)
 *   matrix_max  largest size for the matrix    (default 10000)
 *   filter      only run cases whose name contains this text
 *
 * Memory: every case runs in its own process with its address space capped at SCALE_MEM_MB
 * (default 4096, i.e. a 4 GiB machine), so a case that needs more fails as "OUT OF MEMORY"
 * or a clean error instead of taking the machine down; the peak memory of each case is printed.
 *
 * A case that runs longer than SCALE_TIMEOUT seconds (default 300) aborts the
 * program naming the case, so a hang is reported instead of waiting forever.
 */

#define _POSIX_C_SOURCE 200809L
#define _DEFAULT_SOURCE   /* wait4(): per-child peak memory */

#include <math.h>
#include <signal.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>

#include "graph.h"
#include "graph_advanced.h"
#include "graph_algorithms.h"
#include "graph_social.h"

/* ------------------------------------------------------------------ */
/*  Harness                                                             */
/* ------------------------------------------------------------------ */

static size_t g_max_n        = 1000000;
static size_t g_matrix_max   = 10000;
static const char *g_filter  = NULL;
static unsigned g_timeout    = 300;
static const char *g_only_rep = NULL;   /* SCALE_REP=list|matrix: run only that representation */
static size_t   g_mem_mb     = 4096;     /* SCALE_MEM_MB: hard cap on each case's address space */
static size_t   g_min_n      = 0;       /* SCALE_MIN_N: skip sizes below this */

static int         g_failures = 0;   /* failed expectations over the whole run */
static int         g_case_failed;    /* failed expectations in the current case */
static const char *g_case_name = "";
static GraphRep    g_rep;
static size_t      g_n;              /* requested size of the current run */

#define EXPECT(cond)                                                                     \
    do {                                                                                 \
        if (!(cond)) {                                                                   \
            g_case_failed++;                                                             \
            if (g_case_failed <= 5) {                                                    \
                fprintf(stderr, "      FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond);    \
            }                                                                            \
        }                                                                                \
    } while (0)

static bool near(double a, double b)
{
    if (isinf(a) || isinf(b)) return a == b;   /* unreachable == unreachable */
    return fabs(a - b) <= 1e-9 * (1.0 + fabs(a) + fabs(b));
}

static double now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec + (double)t.tv_nsec * 1e-9;
}

#define EXIT_TIMEOUT 250   /* reserved child exit code: the case ran past SCALE_TIMEOUT */

static void on_alarm(int sig)
{
    (void)sig;
    _exit(EXIT_TIMEOUT);
}

static void *xmalloc(size_t bytes)
{
    void *p = malloc(bytes ? bytes : 1);
    if (p == NULL) {
        fprintf(stderr, "out of memory (%zu bytes)\n", bytes);
        exit(4);
    }
    return p;
}

static void *xcalloc(size_t count, size_t size)
{
    void *p = calloc(count ? count : 1, size ? size : 1);
    if (p == NULL) {
        fprintf(stderr, "out of memory (%zu x %zu bytes)\n", count, size);
        exit(4);
    }
    return p;
}

static uint64_t rng_state = 0x9E3779B97F4A7C15ULL;
static uint64_t rnd64(void)
{
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return rng_state;
}
static size_t rnd_below(size_t bound) { return (size_t)(rnd64() % bound); }

/* ------------------------------------------------------------------ */
/*  Graph kinds, built once per (representation, size) and shared        */
/* ------------------------------------------------------------------ */

typedef enum {
    K_RAND,     /* undirected, weights 1..10, connected: random tree + n extra edges  */
    K_RANDD,    /* directed, weights 1..10, 3 random arcs per vertex                  */
    K_BIN01,    /* undirected like K_RAND but weights 0 or 1                          */
    K_PATH,     /* undirected path, unit weights                                      */
    K_CYCLE,    /* undirected cycle, unit weights                                     */
    K_GRID,     /* side x side grid, unit weights (side = floor(sqrt(n)))             */
    K_TREE,     /* random tree, unit weights                                          */
    K_STAR,     /* star with centre 0                                                 */
    K_DPATH,    /* directed path 0 -> 1 -> ... , unit weights                         */
    K_DCYCLE,   /* directed cycle, unit weights                                       */
    K_DAG,      /* DAG: i -> i+1 and n random forward arcs, weights 1..9              */
    K_DAGNEG,   /* the same arcs with weights -3 .. 5                                 */
    K_MULTI,    /* n/4 disjoint paths of 4 vertices                                   */
    K_CHAIN5,   /* chain of K5 cliques joined by single edges                         */
    K_RING5,    /* ring of K5 cliques (needs >= 3 cliques)                            */
    K_FLOW,     /* source 0, sink n-1, n-2 middle vertices with random capacities     */
    K_SOCIAL,   /* social network: power-law friendships, 20% of users have none      */
    K_COUNT
} Kind;

typedef struct {
    Graph *g[K_COUNT];
    size_t actual[K_COUNT];   /* vertex count actually built (grids and chains round) */
} Pool;

static Pool g_pool;

static size_t grid_side(size_t n)
{
    size_t s = (size_t)sqrt((double)n);
    while (s * s > n) s--;
    while ((s + 1) * (s + 1) <= n) s++;
    return s;
}

static Graph *new_graph(size_t n, bool directed)
{
    Graph *g = graph_create_sized(n, directed);
    if (g == NULL) {
        fprintf(stderr, "cannot create a %zu-vertex graph\n", n);
        exit(4);
    }
    graph_set_adaptive(g, false);
    return g;
}

static void edge(Graph *g, size_t u, size_t v, double w)
{
    (void)graph_add_edge(g, (int)u, (int)v, w);   /* duplicates are simply refused */
}

static void clique5(Graph *g, size_t base)
{
    for (size_t a = 0; a < 5; a++)
        for (size_t b = a + 1; b < 5; b++) edge(g, base + a, base + b, 1.0);
}

static Graph *build(Kind k, size_t n, size_t *actual)
{
    Graph *g = NULL;
    *actual = n;
    switch (k) {
    case K_RAND:
    case K_BIN01:
        g = new_graph(n, false);
        for (size_t i = 1; i < n; i++)
            edge(g, i, rnd_below(i), k == K_RAND ? (double)(1 + rnd_below(10)) : (double)rnd_below(2));
        for (size_t i = 0; i < n; i++) {
            size_t a = rnd_below(n), b = rnd_below(n);
            if (a != b) edge(g, a, b, k == K_RAND ? (double)(1 + rnd_below(10)) : (double)rnd_below(2));
        }
        break;
    case K_RANDD:
        g = new_graph(n, true);
        for (size_t i = 0; i < n; i++)
            for (int j = 0; j < 3; j++) {
                size_t b = rnd_below(n);
                if (b != i) edge(g, i, b, (double)(1 + rnd_below(10)));
            }
        break;
    case K_PATH:
        g = new_graph(n, false);
        for (size_t i = 0; i + 1 < n; i++) edge(g, i, i + 1, 1.0);
        break;
    case K_CYCLE:
        g = new_graph(n, false);
        for (size_t i = 0; i < n; i++) edge(g, i, (i + 1) % n, 1.0);
        break;
    case K_GRID: {
        size_t s = grid_side(n);
        *actual  = s * s;
        g        = new_graph(s * s, false);
        for (size_t r = 0; r < s; r++)
            for (size_t c = 0; c < s; c++) {
                if (c + 1 < s) edge(g, r * s + c, r * s + c + 1, 1.0);
                if (r + 1 < s) edge(g, r * s + c, (r + 1) * s + c, 1.0);
            }
        break;
    }
    case K_TREE:
        g = new_graph(n, false);
        for (size_t i = 1; i < n; i++) edge(g, i, rnd_below(i), 1.0);
        break;
    case K_STAR:
        g = new_graph(n, false);
        for (size_t i = 1; i < n; i++) edge(g, 0, i, 1.0);
        break;
    case K_DPATH:
        g = new_graph(n, true);
        for (size_t i = 0; i + 1 < n; i++) edge(g, i, i + 1, 1.0);
        break;
    case K_DCYCLE:
        g = new_graph(n, true);
        for (size_t i = 0; i < n; i++) edge(g, i, (i + 1) % n, 1.0);
        break;
    case K_DAG:
    case K_DAGNEG: {
        g = new_graph(n, true);
        const uint64_t saved = rng_state;
        rng_state            = 0x1234567ULL;   /* K_DAG and K_DAGNEG share their arcs */
        for (size_t i = 0; i + 1 < n; i++) {
            const double w = (double)(1 + rnd_below(9));
            edge(g, i, i + 1, k == K_DAG ? w : w - 4.0);
        }
        for (size_t t = 0; t < n; t++) {
            size_t a = rnd_below(n);
            size_t b = a + 2 + rnd_below(50);
            const double w = (double)(1 + rnd_below(9));
            if (b < n) edge(g, a, b, k == K_DAG ? w : w - 4.0);
        }
        rng_state = saved;
        break;
    }
    case K_MULTI:
        n         = n / 4 * 4;
        *actual   = n;
        g         = new_graph(n, false);
        for (size_t i = 0; i < n; i += 4) {
            edge(g, i, i + 1, 1.0);
            edge(g, i + 1, i + 2, 1.0);
            edge(g, i + 2, i + 3, 1.0);
        }
        break;
    case K_CHAIN5:
    case K_RING5: {
        const size_t cliques = n / 5;
        n                    = cliques * 5;
        *actual              = n;
        g                    = new_graph(n, false);
        for (size_t c = 0; c < cliques; c++) {
            clique5(g, c * 5);
            if (c + 1 < cliques) edge(g, c * 5 + 4, (c + 1) * 5, 1.0);
        }
        if (k == K_RING5 && cliques >= 3) edge(g, cliques * 5 - 1, 0, 1.0);
        break;
    }
    case K_SOCIAL: {
        /* Preferential attachment: a new user befriends people with probability proportional to
         * how many friends they already have, so a few hubs collect most of the links (a
         * power-law degree distribution).  Every fifth user joins but never connects to anyone. */
        g = new_graph(n, false);
        size_t cap = 4 * n + 16, ends_n = 0;
        int *ends = xmalloc(cap * sizeof *ends);      /* one entry per edge endpoint */
        size_t active = 0;
        int   *active_ids = xmalloc((n > 0 ? n : 1) * sizeof *active_ids);
        for (size_t v = 0; v < n; v++) {
            if (v % 5 == 4) continue;                 /* the silent 20% */
            const size_t links = 1 + (size_t)(-log((double)(rnd64() % 1000003 + 1) / 1000004.0) * 4.0);
            for (size_t k = 0; k < links && active > 0; k++) {
                const int friend_id = (ends_n > 0 && rnd_below(10) < 8)
                                          ? ends[rnd_below(ends_n)]
                                          : active_ids[rnd_below(active)];
                if (friend_id != (int)v && graph_add_edge(g, (int)v, friend_id, 1.0) && ends_n + 2 <= cap) {
                    ends[ends_n++] = (int)v;
                    ends[ends_n++] = friend_id;
                }
            }
            active_ids[active++] = (int)v;
        }
        free(ends);
        free(active_ids);
        break;
    }
    case K_FLOW:
        g = new_graph(n, true);
        for (size_t m = 1; m + 1 < n; m++) {
            edge(g, 0, m, (double)(1 + rnd_below(5)));
            edge(g, m, n - 1, (double)(1 + rnd_below(5)));
        }
        break;
    default:
        break;
    }
    if (g != NULL && g_rep == GRAPH_MATRIX) {
        if (!graph_convert(g, GRAPH_MATRIX)) {
            fprintf(stderr, "matrix conversion failed (kind %d, n=%zu)\n", (int)k, n);
            exit(4);
        }
    }
    return g;
}

static Graph *get(Kind k)
{
    if (g_pool.g[k] == NULL) {
        rng_state = 0x9E3779B97F4A7C15ULL + (uint64_t)k * 7919u + g_n;   /* same graph on every run */
        g_pool.g[k] = build(k, g_n, &g_pool.actual[k]);
    }
    return g_pool.g[k];
}

static size_t size_of(Kind k)
{
    (void)get(k);
    return g_pool.actual[k];
}

/* ------------------------------------------------------------------ */
/*  CSR snapshot (from graph_neighbors only) and oracles                 */
/* ------------------------------------------------------------------ */

typedef struct {
    size_t  n;      /* id bound */
    size_t  m;      /* arcs stored (an undirected edge is two arcs) */
    size_t *off;
    int    *adj;
    double *w;
} Csr;

static Csr csr_of(Graph *g)
{
    Csr c;
    c.n   = graph_id_bound(g);
    c.off = xmalloc((c.n + 1) * sizeof *c.off);
    size_t maxdeg = 0, total = 0;
    for (size_t v = 0; v < c.n; v++) {
        c.off[v] = total;
        const size_t d = graph_degree(g, (int)v);
        total += d;
        if (d > maxdeg) maxdeg = d;
    }
    c.off[c.n] = total;
    c.m        = total;
    c.adj      = xmalloc(total * sizeof *c.adj);
    c.w        = xmalloc(total * sizeof *c.w);
    for (size_t v = 0; v < c.n; v++) {
        const size_t got = graph_neighbors(g, (int)v, c.adj + c.off[v], c.w + c.off[v],
                                           c.off[v + 1] - c.off[v]);
        if (got != c.off[v + 1] - c.off[v]) {
            fprintf(stderr, "graph_neighbors returned %zu, degree is %zu\n", got,
                    c.off[v + 1] - c.off[v]);
            exit(4);
        }
    }
    (void)maxdeg;
    return c;
}

static void csr_free(Csr *c)
{
    free(c->off);
    free(c->adj);
    free(c->w);
    *c = (Csr){0};
}

/* The CSR of a pooled graph is computed on first use. */
static Csr  g_csr[K_COUNT];
static bool g_csr_ok[K_COUNT];

static const Csr *csr(Kind k)
{
    if (!g_csr_ok[k]) {
        g_csr[k]    = csr_of(get(k));
        g_csr_ok[k] = true;
    }
    return &g_csr[k];
}

/* Binary min-heap of (key, vertex) for the Dijkstra-style oracles. */
typedef struct {
    double *key;
    int    *val;
    size_t  size, cap;
} Heap;

static void heap_init(Heap *h, size_t cap)
{
    h->key  = xmalloc((cap + 1) * sizeof *h->key);
    h->val  = xmalloc((cap + 1) * sizeof *h->val);
    h->size = 0;
    h->cap  = cap + 1;
}
static void heap_free(Heap *h) { free(h->key); free(h->val); }
static void heap_push(Heap *h, double key, int val)
{
    if (h->size + 1 >= h->cap) {
        h->cap *= 2;
        h->key = realloc(h->key, h->cap * sizeof *h->key);
        h->val = realloc(h->val, h->cap * sizeof *h->val);
        if (h->key == NULL || h->val == NULL) exit(4);
    }
    size_t i = ++h->size;
    while (i > 1 && h->key[i / 2] > key) {
        h->key[i] = h->key[i / 2];
        h->val[i] = h->val[i / 2];
        i /= 2;
    }
    h->key[i] = key;
    h->val[i] = val;
}
static void heap_pop(Heap *h, double *key, int *val)
{
    *key = h->key[1];
    *val = h->val[1];
    const double k = h->key[h->size];
    const int    v = h->val[h->size];
    h->size--;
    size_t i = 1;
    while (2 * i <= h->size) {
        size_t c = 2 * i;
        if (c + 1 <= h->size && h->key[c + 1] < h->key[c]) c++;
        if (h->key[c] >= k) break;
        h->key[i] = h->key[c];
        h->val[i] = h->val[c];
        i = c;
    }
    h->key[i] = k;
    h->val[i] = v;
}

/* Shortest distances from several sources (non-negative weights). */
static void o_dijkstra(const Csr *c, const int *src, size_t ns, double *dist)
{
    Heap h;
    heap_init(&h, c->m + ns + 16);
    for (size_t v = 0; v < c->n; v++) dist[v] = INFINITY;
    for (size_t i = 0; i < ns; i++) {
        dist[src[i]] = 0;
        heap_push(&h, 0, src[i]);
    }
    while (h.size > 0) {
        double d;
        int    u;
        heap_pop(&h, &d, &u);
        if (d > dist[u]) continue;
        for (size_t e = c->off[u]; e < c->off[u + 1]; e++) {
            const double nd = d + c->w[e];
            if (nd < dist[c->adj[e]]) {
                dist[c->adj[e]] = nd;
                heap_push(&h, nd, c->adj[e]);
            }
        }
    }
    heap_free(&h);
}

/* Widest (maximum bottleneck) path from src. */
static void o_widest(const Csr *c, int src, double *width)
{
    Heap h;
    heap_init(&h, c->m + 16);
    for (size_t v = 0; v < c->n; v++) width[v] = -INFINITY;
    width[src] = INFINITY;
    heap_push(&h, -INFINITY, src);   /* key = -width */
    while (h.size > 0) {
        double k;
        int    u;
        heap_pop(&h, &k, &u);
        if (-k < width[u]) continue;
        for (size_t e = c->off[u]; e < c->off[u + 1]; e++) {
            const double nw = fmin(width[u], c->w[e]);
            if (nw > width[c->adj[e]]) {
                width[c->adj[e]] = nw;
                heap_push(&h, -nw, c->adj[e]);
            }
        }
    }
    heap_free(&h);
}

/* Hop distances (BFS); returns the number reached. */
static size_t o_bfs(const Csr *c, int src, int *level)
{
    int *queue = xmalloc(c->n * sizeof *queue);
    for (size_t v = 0; v < c->n; v++) level[v] = -1;
    size_t head = 0, tail = 0;
    level[src]  = 0;
    queue[tail++] = src;
    while (head < tail) {
        const int u = queue[head++];
        for (size_t e = c->off[u]; e < c->off[u + 1]; e++)
            if (level[c->adj[e]] < 0) {
                level[c->adj[e]] = level[u] + 1;
                queue[tail++]    = c->adj[e];
            }
    }
    free(queue);
    return tail;
}

/* Union-find. */
typedef struct { int *p; } Uf;
static Uf uf_make(size_t n)
{
    Uf u = {xmalloc(n * sizeof(int))};
    for (size_t i = 0; i < n; i++) u.p[i] = (int)i;
    return u;
}
static int uf_find(Uf *u, int x)
{
    while (u->p[x] != x) {
        u->p[x] = u->p[u->p[x]];
        x       = u->p[x];
    }
    return x;
}
static bool uf_union(Uf *u, int a, int b)
{
    a = uf_find(u, a);
    b = uf_find(u, b);
    if (a == b) return false;
    u->p[a] = b;
    return true;
}
static void uf_free(Uf *u) { free(u->p); }

/* Kosaraju (iterative): component id per vertex; returns the count. */
static size_t o_scc(const Csr *c, int *comp)
{
    const size_t n = c->n;
    /* transpose */
    size_t *toff = xcalloc(n + 2, sizeof *toff);
    for (size_t e = 0; e < c->m; e++) toff[c->adj[e] + 2]++;
    for (size_t v = 0; v < n; v++) toff[v + 2] += toff[v + 1];
    int *tadj = xmalloc((c->m + 1) * sizeof *tadj);
    for (size_t u = 0; u < n; u++)
        for (size_t e = c->off[u]; e < c->off[u + 1]; e++) tadj[toff[c->adj[e] + 1]++] = (int)u;
    /* first pass: finish order */
    int    *order = xmalloc(n * sizeof *order);
    size_t  cnt   = 0;
    char   *seen  = xcalloc(n, 1);
    int    *stack = xmalloc(n * sizeof *stack);
    size_t *it    = xmalloc(n * sizeof *it);
    for (size_t s = 0; s < n; s++) {
        if (seen[s]) continue;
        size_t sp = 0;
        stack[sp++] = (int)s;
        seen[s]     = 1;
        it[s]       = c->off[s];
        while (sp > 0) {
            const int u = stack[sp - 1];
            if (it[u] < c->off[u + 1]) {
                const int v = c->adj[it[u]++];
                if (!seen[v]) {
                    seen[v]     = 1;
                    it[v]       = c->off[v];
                    stack[sp++] = v;
                }
            } else {
                order[cnt++] = u;
                sp--;
            }
        }
    }
    /* second pass on the transpose in reverse finish order */
    for (size_t v = 0; v < n; v++) comp[v] = -1;
    size_t ncomp = 0;
    for (size_t i = n; i-- > 0;) {
        const int s = order[i];
        if (comp[s] >= 0) continue;
        size_t sp   = 0;
        stack[sp++] = s;
        comp[s]     = (int)ncomp;
        while (sp > 0) {
            const int u = stack[--sp];
            for (size_t e = toff[u]; e < toff[u + 1]; e++)
                if (comp[tadj[e]] < 0) {
                    comp[tadj[e]] = (int)ncomp;
                    stack[sp++]   = tadj[e];
                }
        }
        ncomp++;
    }
    free(toff); free(tadj); free(order); free(seen); free(stack); free(it);
    return ncomp;
}

/* Do two labelings define the same partition? */
static bool same_partition(const int *a, const int *b, size_t n)
{
    int *map = xmalloc(n * sizeof *map), *back = xmalloc(n * sizeof *back);
    for (size_t i = 0; i < n; i++) { map[i] = -1; back[i] = -1; }
    bool ok = true;
    for (size_t v = 0; v < n && ok; v++) {
        if (a[v] < 0 || b[v] < 0 || (size_t)a[v] >= n || (size_t)b[v] >= n) { ok = false; break; }
        if (map[a[v]] == -1 && back[b[v]] == -1) { map[a[v]] = b[v]; back[b[v]] = a[v]; }
        else if (map[a[v]] != b[v] || back[b[v]] != a[v]) ok = false;
    }
    free(map); free(back);
    return ok;
}

/* Are two graphs identical (same vertices, edges and weights)? */
static bool same_graph(Graph *a, Graph *b)
{
    if (graph_vertex_count(a) != graph_vertex_count(b) || graph_edge_count(a) != graph_edge_count(b) ||
        graph_is_directed(a) != graph_is_directed(b) || graph_id_bound(a) != graph_id_bound(b)) {
        return false;
    }
    const Csr ca = csr_of(a);
    bool      ok = true;
    for (size_t u = 0; u < ca.n && ok; u++)
        for (size_t e = ca.off[u]; e < ca.off[u + 1]; e++) {
            double w = 0;
            if (!graph_has_edge(b, (int)u, ca.adj[e], &w) || w != ca.w[e]) { ok = false; break; }
        }
    Csr tmp = ca;
    csr_free(&tmp);
    return ok;
}

/* ================================================================== */
/*  Cases                                                               */
/* ================================================================== */

/* ---- traversal and connectivity ---- */

static void c_bfs(void)
{
    Graph *g = get(K_RAND);
    const Csr *c = csr(K_RAND);
    const size_t n = c->n;
    int *order = xmalloc(n * sizeof *order), *parent = xmalloc(n * sizeof *parent);
    int *level = xmalloc(n * sizeof *level);
    const size_t reach = o_bfs(c, 0, level);
    const size_t got   = graph_bfs(g, 0, order, parent);
    EXPECT(got == reach && reach == n);   /* K_RAND is connected */
    bool ok = true;
    for (size_t i = 0; i < got && ok; i++) {
        const int v = order[i];
        if (i > 0 && level[order[i - 1]] > level[v]) ok = false;   /* non-decreasing levels */
        if (v != 0 && (parent[v] < 0 || level[parent[v]] + 1 != level[v])) ok = false;
    }
    EXPECT(ok);
    free(order); free(parent); free(level);
}

static void c_dfs(void)
{
    Graph *g = get(K_MULTI);
    const size_t n = size_of(K_MULTI);
    int *order = xmalloc(n * sizeof *order), *parent = xmalloc(n * sizeof *parent);
    int *pos   = xmalloc(n * sizeof *pos);
    for (size_t i = 0; i < n; i++) pos[i] = -1;
    const size_t got = graph_dfs(g, 1, order, parent);
    EXPECT(got == 4);   /* one 4-vertex component */
    bool ok = true;
    for (size_t i = 0; i < got; i++) {
        if (pos[order[i]] != -1) ok = false;   /* each vertex once */
        pos[order[i]] = (int)i;
    }
    for (size_t i = 1; i < got; i++)
        if (parent[order[i]] < 0 || pos[parent[order[i]]] < 0 || pos[parent[order[i]]] >= (int)i) ok = false;
    EXPECT(ok);
    free(order); free(parent); free(pos);

    /* a spanning DFS of the big connected graph */
    Graph *r = get(K_RAND);
    const size_t m = size_of(K_RAND);
    order = xmalloc(m * sizeof *order); parent = xmalloc(m * sizeof *parent);
    EXPECT(graph_dfs(r, 0, order, parent) == m);
    char *seen = xcalloc(m, 1);
    for (size_t i = 0; i < m; i++) { EXPECT(!seen[order[i]]); seen[order[i]] = 1; }
    free(order); free(parent); free(seen);
}

static void c_components(void)
{
    const Kind kinds[] = {K_MULTI, K_RAND, K_RANDD};
    for (size_t t = 0; t < 3; t++) {
        if (kinds[t] == K_RANDD) continue;   /* components of a digraph are not defined here */
        Graph *g = get(kinds[t]);
        const Csr *c = csr(kinds[t]);
        int *comp = xmalloc(c->n * sizeof *comp);
        const size_t got = graph_connected_components(g, comp);
        Uf uf = uf_make(c->n);
        size_t expect = c->n;
        for (size_t u = 0; u < c->n; u++)
            for (size_t e = c->off[u]; e < c->off[u + 1]; e++)
                if (uf_union(&uf, (int)u, c->adj[e])) expect--;
        EXPECT(got == expect);
        bool ok = true;
        for (size_t u = 0; u < c->n && ok; u++)
            for (size_t e = c->off[u]; e < c->off[u + 1]; e++)
                if (comp[u] != comp[c->adj[e]]) { ok = false; break; }
        EXPECT(ok);
        EXPECT(kinds[t] != K_MULTI || got == c->n / 4);
        uf_free(&uf);
        free(comp);
    }
}

static void c_predicates(void)
{
    EXPECT(graph_is_tree(get(K_TREE)) && graph_is_forest(get(K_TREE)) && graph_is_connected(get(K_TREE)));
    EXPECT(!graph_has_cycle(get(K_TREE)));
    EXPECT(graph_has_cycle(get(K_CYCLE)) && !graph_is_tree(get(K_CYCLE)) && !graph_is_forest(get(K_CYCLE)));
    EXPECT(graph_is_connected(get(K_CYCLE)));
    EXPECT(!graph_is_connected(get(K_MULTI)) && graph_is_forest(get(K_MULTI)) && !graph_is_tree(get(K_MULTI)));
    EXPECT(graph_has_cycle(get(K_RAND)));
    EXPECT(graph_is_dag(get(K_DPATH)) && !graph_has_cycle(get(K_DPATH)) && !graph_is_strongly_connected(get(K_DPATH)));
    EXPECT(!graph_is_dag(get(K_DCYCLE)) && graph_has_cycle(get(K_DCYCLE)) && graph_is_strongly_connected(get(K_DCYCLE)));
    EXPECT(graph_is_dag(get(K_DAG)));

    bool bip = false;
    int *color = xmalloc(size_of(K_GRID) * sizeof *color);
    EXPECT(graph_is_bipartite(get(K_GRID), &bip, color) == GRAPH_OK && bip);
    const Csr *c = csr(K_GRID);
    bool ok = true;
    for (size_t u = 0; u < c->n && ok; u++)
        for (size_t e = c->off[u]; e < c->off[u + 1]; e++) if (color[u] == color[c->adj[e]]) { ok = false; break; }
    EXPECT(ok);
    free(color);
    EXPECT(graph_is_bipartite(get(K_CYCLE), &bip, NULL) == GRAPH_OK && bip == (size_of(K_CYCLE) % 2 == 0));
    EXPECT(graph_is_bipartite(get(K_CHAIN5), &bip, NULL) == GRAPH_OK && !bip);   /* K5 has odd cycles */

    /* the cycle itself */
    const size_t n = size_of(K_DCYCLE);
    int *cyc = xmalloc(n * sizeof *cyc);
    const size_t len = graph_find_cycle(get(K_DCYCLE), cyc);
    EXPECT(len == n);
    ok = true;
    for (size_t i = 0; i < len; i++) if ((cyc[i] + 1) % (int)n != cyc[(i + 1) % len]) ok = false;
    EXPECT(ok);
    free(cyc);
}

static void c_topo(void)
{
    Graph *g = get(K_DAG);
    const Csr *c = csr(K_DAG);
    const size_t n = c->n;
    int *order = xmalloc(n * sizeof *order), *pos = xmalloc(n * sizeof *pos);
    EXPECT(graph_topological_sort(g, order));
    for (size_t i = 0; i < n; i++) pos[order[i]] = (int)i;
    bool ok = true;
    for (size_t u = 0; u < n && ok; u++)
        for (size_t e = c->off[u]; e < c->off[u + 1]; e++) if (pos[u] >= pos[c->adj[e]]) { ok = false; break; }
    EXPECT(ok);

    /* longest path by DP in id order (all arcs go from a smaller to a larger id) */
    double *best = xmalloc(n * sizeof *best);
    for (size_t v = 0; v < n; v++) best[v] = 0;
    double top = 0;
    for (size_t u = 0; u < n; u++) {
        for (size_t e = c->off[u]; e < c->off[u + 1]; e++)
            if (best[u] + c->w[e] > best[c->adj[e]]) best[c->adj[e]] = best[u] + c->w[e];
        if (best[u] > top) top = best[u];
    }
    double *dist = xmalloc(n * sizeof *dist);
    int *parent = xmalloc(n * sizeof *parent);
    int end = -1;
    EXPECT(graph_dag_longest_path(g, dist, parent, &end) == GRAPH_OK);
    ok = true;
    for (size_t v = 0; v < n; v++) if (dist[v] != best[v]) { ok = false; break; }
    EXPECT(ok && end >= 0 && dist[end] == top);

    int *path = xmalloc(n * sizeof *path);
    size_t plen = 0;
    double pcost = 0;
    EXPECT(graph_critical_path(g, path, &plen, &pcost) == GRAPH_OK);
    EXPECT(pcost == top && plen >= 2);
    double sum = 0;
    ok = true;
    for (size_t i = 0; i + 1 < plen; i++) {
        double w;
        if (!graph_has_edge(g, path[i], path[i + 1], &w)) { ok = false; break; }
        sum += w;
    }
    EXPECT(ok && sum == top);

    /* the advanced single-source version, from vertex 0 (reaches everything via i -> i+1) */
    double *adist = xmalloc(n * sizeof *adist);
    EXPECT(graph_adv_longest_path(g, 0, adist, NULL) == GRAPH_OK);
    ok = true;
    for (size_t v = 0; v < n; v++) if (adist[v] != best[v]) { ok = false; break; }
    EXPECT(ok);
    free(order); free(pos); free(best); free(dist); free(parent); free(path); free(adist);
}

static void c_scc(void)
{
    Graph *g = get(K_RANDD);
    const Csr *c = csr(K_RANDD);
    int *mine = xmalloc(c->n * sizeof *mine), *theirs = xmalloc(c->n * sizeof *theirs);
    const size_t want = o_scc(c, theirs);
    size_t got = 0;
    EXPECT(graph_strongly_connected_components(g, mine, &got) == GRAPH_OK);
    EXPECT(got == want);
    EXPECT(same_partition(mine, theirs, c->n));
    free(mine); free(theirs);
    EXPECT(graph_is_strongly_connected(get(K_DCYCLE)));
}

static void c_blocks(void)
{
    /* chain of k K5 cliques: k-1 bridges, 2(k-1) cut vertices, 2k-1 blocks */
    {
        Graph *g = get(K_CHAIN5);
        const size_t n = size_of(K_CHAIN5), k = n / 5;
        bool *cut = xmalloc(n * sizeof *cut);
        size_t ncut = 0;
        EXPECT(graph_articulation_points(g, cut, &ncut) == GRAPH_OK);
        EXPECT(ncut == 2 * (k - 1));
        GraphEdge *br = xmalloc(graph_edge_count(g) * sizeof *br);
        size_t nbr = 0;
        EXPECT(graph_bridges(g, br, &nbr) == GRAPH_OK && nbr == k - 1);
        GraphGroups blocks = {0};
        EXPECT(graph_biconnected_components(g, &blocks) == GRAPH_OK && blocks.count == 2 * k - 1);
        graph_groups_free(&blocks);
        free(cut); free(br);
    }
    /* path: every edge a bridge, every inner vertex a cut vertex */
    {
        Graph *g = get(K_PATH);
        const size_t n = size_of(K_PATH);
        bool *cut = xmalloc(n * sizeof *cut);
        size_t ncut = 0;
        EXPECT(graph_articulation_points(g, cut, &ncut) == GRAPH_OK && ncut == n - 2);
        GraphEdge *br = xmalloc(graph_edge_count(g) * sizeof *br);
        size_t nbr = 0;
        EXPECT(graph_bridges(g, br, &nbr) == GRAPH_OK && nbr == n - 1);
        free(cut); free(br);
    }
    /* cycle: neither */
    {
        Graph *g = get(K_CYCLE);
        const size_t n = size_of(K_CYCLE);
        bool *cut = xmalloc(n * sizeof *cut);
        size_t ncut = 9;
        EXPECT(graph_articulation_points(g, cut, &ncut) == GRAPH_OK && ncut == 0);
        GraphEdge *br = xmalloc(graph_edge_count(g) * sizeof *br);
        size_t nbr = 9;
        EXPECT(graph_bridges(g, br, &nbr) == GRAPH_OK && nbr == 0);
        GraphGroups blocks = {0};
        EXPECT(graph_biconnected_components(g, &blocks) == GRAPH_OK && blocks.count == 1);
        graph_groups_free(&blocks);
        free(cut); free(br);
    }
    /* a random tree: all bridges */
    {
        Graph *g = get(K_TREE);
        GraphEdge *br = xmalloc(graph_edge_count(g) * sizeof *br);
        size_t nbr = 0;
        EXPECT(graph_bridges(g, br, &nbr) == GRAPH_OK && nbr == size_of(K_TREE) - 1);
        free(br);
    }
}

/* ---- shortest paths ---- */

static void check_dist(Kind k, const double *dist, const int *parent, int src)
{
    Graph *g = get(k);
    const Csr *c = csr(k);
    double *want = xmalloc(c->n * sizeof *want);
    o_dijkstra(c, &src, 1, want);
    bool ok = true;
    for (size_t v = 0; v < c->n; v++) if (!near(dist[v], want[v])) { ok = false; break; }
    EXPECT(ok);
    ok = true;
    for (size_t v = 0; v < c->n && ok; v++) {
        if ((int)v == src || isinf(want[v])) continue;
        double w;
        if (parent[v] < 0 || !graph_has_edge(g, parent[v], (int)v, &w) || !near(dist[parent[v]] + w, dist[v])) ok = false;
    }
    EXPECT(ok);   /* every parent edge is tight */
    free(want);
}

static void c_dijkstra(void)
{
    const Kind kinds[] = {K_RAND, K_RANDD};
    for (int t = 0; t < 2; t++) {
        Graph *g = get(kinds[t]);
        const size_t n = size_of(kinds[t]);
        double *dist = xmalloc(n * sizeof *dist);
        int *parent  = xmalloc(n * sizeof *parent);
        EXPECT(graph_dijkstra(g, 0, dist, parent));
        check_dist(kinds[t], dist, parent, 0);
        free(dist); free(parent);
    }
}

static double manhattan_h(int v, int goal, void *ctx)
{
    const int side = *(const int *)ctx;
    return (double)(abs(v / side - goal / side) + abs(v % side - goal % side));
}

static void c_astar(void)
{
    /* grid, admissible Manhattan heuristic: cost is the Manhattan distance */
    {
        Graph *g = get(K_GRID);
        const size_t n = size_of(K_GRID);
        int side = (int)grid_side(n);
        GraphPath p = {0};
        EXPECT(graph_a_star(g, 0, (int)n - 1, manhattan_h, &side, &p) == GRAPH_OK);
        EXPECT(near(p.cost, 2.0 * (side - 1)) && p.length == (size_t)(2 * (side - 1) + 1));
        graph_path_free(&p);
    }
    /* random weighted graph: same cost as Dijkstra, with and without a heuristic */
    {
        Graph *g = get(K_RAND);
        const Csr *c = csr(K_RAND);
        double *want = xmalloc(c->n * sizeof *want);
        const int src = 0, dst = (int)c->n - 1;
        o_dijkstra(c, &src, 1, want);
        GraphPath p = {0};
        EXPECT(graph_a_star(g, src, dst, NULL, NULL, &p) == GRAPH_OK && near(p.cost, want[dst]));
        double sum = 0;
        bool ok = p.length >= 2 && p.vertices[0] == src && p.vertices[p.length - 1] == dst;
        for (size_t i = 0; ok && i + 1 < p.length; i++) {
            double w;
            if (!graph_has_edge(g, p.vertices[i], p.vertices[i + 1], &w)) ok = false; else sum += w;
        }
        EXPECT(ok && near(sum, want[dst]));
        graph_path_free(&p);
        free(want);
    }
}

static void c_01bfs(void)
{
    Graph *g = get(K_BIN01);
    const size_t n = size_of(K_BIN01);
    double *dist = xmalloc(n * sizeof *dist);
    int *parent = xmalloc(n * sizeof *parent);
    EXPECT(graph_zero_one_bfs(g, 0, dist, parent) == GRAPH_OK);
    check_dist(K_BIN01, dist, parent, 0);
    free(dist); free(parent);
}

static void c_multi(void)
{
    Graph *g = get(K_RAND);
    const Csr *c = csr(K_RAND);
    int src[8];
    for (int i = 0; i < 8; i++) src[i] = (int)rnd_below(c->n);
    double *want = xmalloc(c->n * sizeof *want), *dist = xmalloc(c->n * sizeof *dist);
    int *parent = xmalloc(c->n * sizeof *parent), *origin = xmalloc(c->n * sizeof *origin);
    o_dijkstra(c, src, 8, want);
    EXPECT(graph_dijkstra_multi(g, src, 8, dist, parent, origin) == GRAPH_OK);
    bool ok = true;
    for (size_t v = 0; v < c->n; v++) if (!near(dist[v], want[v])) { ok = false; break; }
    EXPECT(ok);
    ok = true;
    for (size_t v = 0; v < c->n; v++) {
        bool is_src = false;
        for (int i = 0; i < 8; i++) if (origin[v] == src[i]) is_src = true;
        if (!is_src) { ok = false; break; }
    }
    EXPECT(ok);
    free(want); free(dist); free(parent); free(origin);
}

static void c_bellman(void)
{
    /* positive weights: must equal Dijkstra */
    {
        Graph *g = get(K_RANDD);
        const size_t n = size_of(K_RANDD);
        double *dist = xmalloc(n * sizeof *dist);
        int *parent = xmalloc(n * sizeof *parent);
        EXPECT(graph_bellman_ford(g, 0, dist, parent) == GRAPH_OK);
        check_dist(K_RANDD, dist, parent, 0);
        free(dist); free(parent);
    }
    /* negative arcs in a DAG: DP over id order is the oracle */
    {
        Graph *g = get(K_DAGNEG);
        const Csr *c = csr(K_DAGNEG);
        double *want = xmalloc(c->n * sizeof *want), *dist = xmalloc(c->n * sizeof *dist);
        int *parent = xmalloc(c->n * sizeof *parent);
        for (size_t v = 0; v < c->n; v++) want[v] = INFINITY;
        want[0] = 0;
        for (size_t u = 0; u < c->n; u++)
            for (size_t e = c->off[u]; e < c->off[u + 1]; e++)
                if (want[u] + c->w[e] < want[c->adj[e]]) want[c->adj[e]] = want[u] + c->w[e];
        EXPECT(graph_bellman_ford(g, 0, dist, parent) == GRAPH_OK);
        bool ok = true;
        for (size_t v = 0; v < c->n; v++) if (!near(dist[v], want[v])) { ok = false; break; }
        EXPECT(ok);
        free(want); free(dist); free(parent);
    }
}

/* ---- trees, flow, matching ---- */

static double o_kruskal(const Csr *c)
{
    typedef struct { double w; int a, b; } E;
    E *es = xmalloc((c->m / 2 + 1) * sizeof *es);
    size_t m = 0;
    for (size_t u = 0; u < c->n; u++)
        for (size_t e = c->off[u]; e < c->off[u + 1]; e++)
            if ((int)u < c->adj[e]) es[m++] = (E){c->w[e], (int)u, c->adj[e]};
    /* sort by weight with a counting pass: weights are small integers */
    size_t cnt[12] = {0};
    for (size_t i = 0; i < m; i++) cnt[(int)es[i].w]++;
    E *sorted = xmalloc((m + 1) * sizeof *sorted);
    size_t start[12] = {0};
    for (int w = 1; w < 12; w++) start[w] = start[w - 1] + cnt[w - 1];
    for (size_t i = 0; i < m; i++) sorted[start[(int)es[i].w]++] = es[i];
    Uf uf = uf_make(c->n);
    double total = 0;
    for (size_t i = 0; i < m; i++) if (uf_union(&uf, sorted[i].a, sorted[i].b)) total += sorted[i].w;
    uf_free(&uf); free(es); free(sorted);
    return total;
}

static void c_mst(void)
{
    Graph *g = get(K_RAND);
    const Csr *c = csr(K_RAND);
    const size_t n = c->n;
    const double want = o_kruskal(c);
    GraphEdge *edges = xmalloc(n * sizeof *edges);
    size_t cnt = 0;
    const double kr = graph_mst_kruskal(g, edges, &cnt);
    EXPECT(near(kr, want) && cnt == n - 1);
    Uf uf = uf_make(n);
    bool ok = true;
    double sum = 0;
    for (size_t i = 0; i < cnt; i++) {
        double w;
        if (!graph_has_edge(g, edges[i].u, edges[i].v, &w) ||
            !uf_union(&uf, edges[i].u, edges[i].v)) ok = false;
        sum += w;
    }
    EXPECT(ok && near(sum, want));   /* edges exist, form a tree, and weigh the optimum */
    uf_free(&uf);
    cnt = 0;
    EXPECT(near(graph_mst_prim(g, edges, &cnt), want) && cnt == n - 1);
    free(edges);
}

static void c_flow(void)
{
    /* source -> middle -> sink with capacities a, b: the flow is sum(min(a, b)) */
    Graph *g = get(K_FLOW);
    const Csr *c = csr(K_FLOW);
    const size_t n = c->n;
    double want = 0;
    for (size_t m = 1; m + 1 < n; m++) {
        double a = 0, b = 0;
        (void)graph_has_edge(g, 0, (int)m, &a);
        (void)graph_has_edge(g, (int)m, (int)n - 1, &b);
        want += fmin(a, b);
    }
    double value = 0;
    EXPECT(graph_max_flow(g, 0, (int)n - 1, GRAPH_FLOW_DINIC, &value, NULL) == GRAPH_OK && near(value, want));
    if (n <= 2000) {
        EXPECT(graph_max_flow(g, 0, (int)n - 1, GRAPH_FLOW_EDMONDS_KARP, &value, NULL) == GRAPH_OK && near(value, want));
    }
    GraphEdge *cut = xmalloc((graph_edge_count(g) + 1) * sizeof *cut);
    size_t ncut = 0;
    EXPECT(graph_min_cut(g, 0, (int)n - 1, &value, cut, &ncut) == GRAPH_OK && near(value, want));
    double capsum = 0;
    bool ok = true;
    for (size_t i = 0; i < ncut; i++) {
        double w;
        if (!graph_has_edge(g, cut[i].u, cut[i].v, &w)) ok = false; else capsum += w;
    }
    EXPECT(ok && near(capsum, want));   /* the cut's edges weigh exactly the flow */
    free(cut);

    /* undirected unit grid, opposite corners: both corners have degree 2 */
    Graph *grid = get(K_GRID);
    const size_t gn = size_of(K_GRID);
    EXPECT(graph_max_flow(grid, 0, (int)gn - 1, GRAPH_FLOW_DINIC, &value, NULL) == GRAPH_OK && near(value, 2.0));
}

static void c_matching(void)
{
    Graph *g = get(K_GRID);
    const size_t n = size_of(K_GRID);
    int *match = xmalloc(n * sizeof *match);
    size_t size = 0;
    EXPECT(graph_bipartite_matching(g, match, &size) == GRAPH_OK);
    EXPECT(size == n / 2);   /* a grid always has a maximum matching of floor(n/2) */
    bool ok = true;
    size_t pairs = 0;
    for (size_t v = 0; v < n; v++) {
        if (match[v] < 0) continue;
        pairs++;
        if (match[match[v]] != (int)v || !graph_has_edge(g, (int)v, match[v], NULL)) ok = false;
    }
    EXPECT(ok && pairs == 2 * size);
    free(match);
    /* not bipartite -> refused */
    match = xmalloc(size_of(K_CHAIN5) * sizeof *match);
    EXPECT(graph_bipartite_matching(get(K_CHAIN5), match, &size) == GRAPH_ERR_UNSUPPORTED);
    free(match);
}

/* ---- combinatorial ---- */

static void check_coloring(Kind k, const int *color, size_t ncolors, bool tight)
{
    const Csr *c = csr(k);
    size_t maxdeg = 0;
    bool ok = true;
    for (size_t u = 0; u < c->n; u++) {
        if (c->off[u + 1] - c->off[u] > maxdeg) maxdeg = c->off[u + 1] - c->off[u];
        if (color[u] < 0 || (size_t)color[u] >= ncolors) ok = false;
        for (size_t e = c->off[u]; ok && e < c->off[u + 1]; e++) if (color[u] == color[c->adj[e]]) ok = false;
    }
    EXPECT(ok);                       /* proper: no edge joins two vertices of one colour */
    EXPECT(ncolors <= maxdeg + 1);    /* the greedy guarantee */
    (void)tight;
}

typedef GraphStatus (*ColorFn)(Graph *, int *, size_t *);

static void run_coloring(ColorFn fn, bool exact_on_bipartite)
{
    const Kind kinds[] = {K_RAND, K_GRID, K_CHAIN5};
    for (int t = 0; t < 3; t++) {
        Graph *g = get(kinds[t]);
        const size_t n = size_of(kinds[t]);
        int *color = xmalloc(n * sizeof *color);
        size_t nc = 0;
        EXPECT(fn(g, color, &nc) == GRAPH_OK);
        check_coloring(kinds[t], color, nc, false);
        if (exact_on_bipartite && kinds[t] == K_GRID) EXPECT(nc == 2);   /* DSatur is exact on bipartite graphs */
        if (kinds[t] == K_CHAIN5) EXPECT(nc >= 5);                      /* K5 needs 5 colours */
        free(color);
    }
}
static void c_color_greedy(void)  { run_coloring(graph_color_greedy, false); }
static void c_color_wp(void)      { run_coloring(graph_color_welsh_powell, false); }
static void c_color_dsatur(void)  { run_coloring(graph_color_dsatur, true); }

static void c_cover_independent(void)
{
    const Kind kinds[] = {K_RAND, K_GRID};
    for (int t = 0; t < 2; t++) {
        Graph *g = get(kinds[t]);
        const Csr *c = csr(kinds[t]);
        bool *in = xmalloc(c->n * sizeof *in);
        size_t size = 0;
        EXPECT(graph_vertex_cover_approx(g, in, &size) == GRAPH_OK);
        bool ok = true;
        size_t counted = 0;
        for (size_t u = 0; u < c->n; u++) {
            if (in[u]) counted++;
            for (size_t e = c->off[u]; e < c->off[u + 1]; e++) if (!in[u] && !in[c->adj[e]]) ok = false;
        }
        EXPECT(ok && counted == size);   /* every edge has an end in the cover */
        if (kinds[t] == K_GRID) {
            /* optimum is floor(n/2) (König), a 2-approximation is at most n */
            EXPECT(size >= c->n / 2 && size <= c->n);
        }
        EXPECT(graph_independent_set_greedy(g, in, &size) == GRAPH_OK);
        ok = true;
        counted = 0;
        for (size_t u = 0; u < c->n; u++) {
            if (!in[u]) continue;
            counted++;
            for (size_t e = c->off[u]; e < c->off[u + 1]; e++) if (in[c->adj[e]]) ok = false;
        }
        EXPECT(ok && counted == size);   /* no edge inside the set */
        /* maximal: every outside vertex has a neighbour inside */
        ok = true;
        for (size_t u = 0; u < c->n && ok; u++) {
            if (in[u]) continue;
            bool has = false;
            for (size_t e = c->off[u]; e < c->off[u + 1]; e++) if (in[c->adj[e]]) has = true;
            if (!has && c->off[u + 1] > c->off[u]) ok = false;
        }
        EXPECT(ok);
        free(in);
    }
    /* star: the greedy picks all leaves */
    Graph *s = get(K_STAR);
    const size_t n = size_of(K_STAR);
    bool *in = xmalloc(n * sizeof *in);
    size_t size = 0;
    EXPECT(graph_independent_set_greedy(s, in, &size) == GRAPH_OK && size == n - 1);
    free(in);
}

static void c_maximal_cliques(void)
{
    Graph *g = get(K_CHAIN5);
    const size_t n = size_of(K_CHAIN5), k = n / 5;
    GraphGroups cl = {0};
    EXPECT(graph_maximal_cliques(g, &cl) == GRAPH_OK);
    EXPECT(cl.count == 2 * k - 1);   /* k K5 cliques and k-1 connecting edges */
    size_t fives = 0, twos = 0;
    for (size_t i = 0; i < cl.count; i++) {
        const size_t sz = cl.offsets[i + 1] - cl.offsets[i];
        if (sz == 5) fives++; else if (sz == 2) twos++;
    }
    EXPECT(fives == k && twos == k - 1);
    graph_groups_free(&cl);
}

static void c_max_clique(void)
{
    int clique[8];
    size_t size = 0;
    EXPECT(graph_max_clique(get(K_CHAIN5), clique, &size) == GRAPH_OK && size == 5);
}

static void c_triangles_stats(void)
{
    Graph *g = get(K_CHAIN5);
    const size_t n = size_of(K_CHAIN5), k = n / 5;
    long long tri = 0;
    EXPECT(graph_triangle_count(g, &tri) == GRAPH_OK && tri == (long long)(10 * k));
    double *local = xmalloc(n * sizeof *local), avg = 0;
    EXPECT(graph_clustering_coefficient(g, local, &avg) == GRAPH_OK);
    const double connectors = 2.0 * (double)(k - 1);
    EXPECT(near(avg, ((double)n - connectors + 0.6 * connectors) / (double)n));
    free(local);

    int common[8];
    size_t ncommon = 0;
    EXPECT(graph_common_neighbors(g, 1, 2, common, &ncommon) == GRAPH_OK && ncommon == 3);
    double jac = 0;
    EXPECT(graph_jaccard(g, 1, 2, &jac) == GRAPH_OK && near(jac, 3.0 / 5.0));

    size_t hist[8] = {0};
    EXPECT(graph_degree_histogram(g, hist, 8) == GRAPH_OK);
    EXPECT(hist[4] == n - (size_t)connectors && hist[5] == (size_t)connectors);

    GraphNetworkStats st;
    EXPECT(graph_network_stats(g, &st) == GRAPH_OK && st.max_degree == 5);
    EXPECT(near(st.density, (double)graph_edge_count(g) / ((double)n * ((double)n - 1) / 2.0)));
}

static void c_euler(void)
{
    const size_t n = size_of(K_CYCLE);
    int *path = xmalloc((n + 2) * sizeof *path);
    size_t len = 0;
    EXPECT(graph_euler_trail(get(K_CYCLE), true, path, &len) == GRAPH_OK && len == n + 1);
    bool ok = path[0] == path[n];
    for (size_t i = 0; ok && i < n; i++) if (!graph_has_edge(get(K_CYCLE), path[i], path[i + 1], NULL)) ok = false;
    EXPECT(ok);
    free(path);

    const size_t pn = size_of(K_PATH);
    path = xmalloc((pn + 2) * sizeof *path);
    EXPECT(graph_euler_trail(get(K_PATH), false, path, &len) == GRAPH_OK && len == pn);
    EXPECT(graph_euler_trail(get(K_PATH), true, path, &len) == GRAPH_ERR_NOT_FOUND);   /* odd degrees */
    free(path);
}

/* ---- ranking and communities ---- */

static void c_pagerank(void)
{
    Graph *g = get(K_RANDD);
    const Csr *c = csr(K_RANDD);
    const size_t n = c->n;
    double *rank = xmalloc(n * sizeof *rank);
    size_t iters = 0;
    EXPECT(graph_pagerank(g, 0.85, 1e-12, 500, rank, &iters) == GRAPH_OK);
    double sum = 0;
    for (size_t v = 0; v < n; v++) sum += rank[v];
    EXPECT(fabs(sum - 1.0) < 1e-6);

    /* oracle: power iteration, transition probability proportional to weight */
    double *r = xmalloc(n * sizeof *r), *nr = xmalloc(n * sizeof *nr), *strength = xcalloc(n, sizeof *strength);
    for (size_t u = 0; u < n; u++) for (size_t e = c->off[u]; e < c->off[u + 1]; e++) strength[u] += c->w[e];
    for (size_t v = 0; v < n; v++) r[v] = 1.0 / (double)n;
    for (int it = 0; it < 500; it++) {
        double dangling = 0;
        for (size_t v = 0; v < n; v++) nr[v] = 0;
        for (size_t u = 0; u < n; u++) {
            if (strength[u] == 0) { dangling += r[u]; continue; }
            for (size_t e = c->off[u]; e < c->off[u + 1]; e++) nr[c->adj[e]] += r[u] * c->w[e] / strength[u];
        }
        double diff = 0;
        for (size_t v = 0; v < n; v++) {
            nr[v] = (1 - 0.85) / (double)n + 0.85 * (nr[v] + dangling / (double)n);
            diff += fabs(nr[v] - r[v]);
        }
        double *t = r; r = nr; nr = t;
        if (diff < 1e-13) break;
    }
    double worst = 0;
    for (size_t v = 0; v < n; v++) worst = fmax(worst, fabs(r[v] - rank[v]));
    EXPECT(worst < 1e-9);
    free(rank); free(r); free(nr); free(strength);

    /* a directed cycle is symmetric: every vertex gets 1/n */
    const size_t cn = size_of(K_DCYCLE);
    rank = xmalloc(cn * sizeof *rank);
    EXPECT(graph_pagerank(get(K_DCYCLE), 0.85, 1e-12, 200, rank, &iters) == GRAPH_OK);
    bool ok = true;
    for (size_t v = 0; v < cn; v++) if (fabs(rank[v] - 1.0 / (double)cn) > 1e-9) ok = false;
    EXPECT(ok);
    free(rank);
}

static void c_eigenvector_compare(void)
{
    /* eigenvector centrality: the advanced version (identity shift) gives the same scores as
     * the existing one on a graph where plain power iteration converges, and is uniform on
     * the (bipartite for even n) cycle where it would oscillate */
    {
        const size_t n = size_of(K_RAND);
        double *a = xmalloc(n * sizeof *a), *b = xmalloc(n * sizeof *b);
        EXPECT(graph_eigenvector_centrality(get(K_RAND), 1e-10, 2000, a) == GRAPH_OK);
        EXPECT(graph_adv_eigenvector_centrality(get(K_RAND), 2000, 1e-10, b, NULL) == GRAPH_OK);
        double worst = 0;
        for (size_t v = 0; v < n; v++) worst = fmax(worst, fabs(a[v] - b[v]));
        EXPECT(worst < 1e-5);
        free(a); free(b);
    }
}

static void c_eigenvector_symmetry(void)
{
    {
        const size_t n = size_of(K_CYCLE);
        double *b = xmalloc(n * sizeof *b);
        EXPECT(graph_adv_eigenvector_centrality(get(K_CYCLE), 200, 1e-9, b, NULL) == GRAPH_OK);
        bool ok = true;
        for (size_t v = 0; v < n; v++) if (fabs(b[v] - b[0]) > 1e-6) ok = false;
        EXPECT(ok);
        free(b);
    }
    {
        const size_t n = size_of(K_STAR);
        double *b = xmalloc(n * sizeof *b);
        EXPECT(graph_adv_eigenvector_centrality(get(K_STAR), 1000, 1e-12, b, NULL) == GRAPH_OK);
        bool ok = true;
        for (size_t v = 1; v < n; v++) if (b[0] <= b[v] || fabs(b[v] - b[1]) > 1e-9) ok = false;
        EXPECT(ok);
        double *deg = xmalloc(n * sizeof *deg);
        EXPECT(graph_degree_centrality(get(K_STAR), deg) == GRAPH_OK && near(deg[0], 1.0));
        free(deg); free(b);
    }
}

static void c_hits(void)
{
    /* HITS on a directed cycle: all hubs equal, all authorities equal */
    {
        const size_t n = size_of(K_DCYCLE);
        double *hub = xmalloc(n * sizeof *hub), *auth = xmalloc(n * sizeof *auth);
        EXPECT(graph_hits(get(K_DCYCLE), 1e-8, 200, hub, auth) == GRAPH_OK);
        bool ok = true;
        for (size_t v = 0; v < n; v++) if (fabs(hub[v] - hub[0]) > 1e-9 || fabs(auth[v] - auth[0]) > 1e-9) ok = false;
        EXPECT(ok);
        free(hub); free(auth);
    }
}

static void c_cores(void)
{
    const struct { Kind k; int core; } want[] = {{K_GRID, 2}, {K_TREE, 1}, {K_CHAIN5, 4}, {K_CYCLE, 2}};
    for (int t = 0; t < 4; t++) {
        const size_t n = size_of(want[t].k);
        int *core = xmalloc(n * sizeof *core);
        EXPECT(graph_core_numbers(get(want[t].k), core) == GRAPH_OK);
        bool ok = true;
        for (size_t v = 0; v < n; v++) if (core[v] != want[t].core) { ok = false; break; }
        EXPECT(ok);
        free(core);
    }
    double r = 0;
    EXPECT(graph_degree_assortativity(get(K_RAND), &r) == GRAPH_OK && r >= -1.0 - 1e-9 && r <= 1.0 + 1e-9);
}

static void c_communities(void)
{
    Graph *g = get(K_RING5);
    const size_t n = size_of(K_RING5), k = n / 5;
    if (k < 3) return;
    int *planted = xmalloc(n * sizeof *planted), *comm = xmalloc(n * sizeof *comm);
    for (size_t v = 0; v < n; v++) planted[v] = (int)(v / 5);
    double q_planted = 0;
    EXPECT(graph_modularity(g, planted, &q_planted) == GRAPH_OK && q_planted > 0.5);

    double q = 0;
    size_t count = 0;
    EXPECT(graph_louvain(g, comm, &q, &count) == GRAPH_OK);
    double q_check = 0;
    EXPECT(graph_modularity(g, comm, &q_check) == GRAPH_OK && near(q, q_check));   /* self-consistent */
    EXPECT(q >= q_planted - 0.02);        /* finds (almost) the planted structure */
    EXPECT(count >= 1 && count <= k);

    EXPECT(graph_label_propagation(g, comm, &count) == GRAPH_OK && count >= 1 && count <= k);
    bool ok = true;
    for (size_t c = 0; c < k && ok; c++)       /* a K5 never gets split by label propagation */
        for (size_t j = 1; j < 5; j++) if (comm[c * 5 + j] != comm[c * 5]) { ok = false; break; }
    EXPECT(ok);
    free(planted); free(comm);
}


static FILE *tmp(void)
{
    FILE *f = tmpfile();
    if (f == NULL) { fprintf(stderr, "tmpfile failed\n"); exit(4); }
    return f;
}

/* ---- a social network: hubs, long tail, users with no connections ---- */

static void c_social(void)
{
    Graph *g = get(K_SOCIAL);
    const Csr *c = csr(K_SOCIAL);
    const size_t n = c->n;

    /* degree structure: every fifth user is silent, the degrees are heavy-tailed */
    size_t isolated = 0, maxdeg = 0, maxdeg_v = 0, total = 0;
    for (size_t v = 0; v < n; v++) {
        const size_t d = c->off[v + 1] - c->off[v];
        total += d;
        if (d == 0) isolated++;
        if (d > maxdeg) { maxdeg = d; maxdeg_v = v; }
    }
    EXPECT(isolated >= n / 5);                           /* at least the silent 20% */
    EXPECT(total == 2 * graph_edge_count(g));
    EXPECT(maxdeg > 10 * (total / n + 1));               /* a hub far above the average: power law */
    size_t *hist = xcalloc(maxdeg + 2, sizeof *hist);
    EXPECT(graph_degree_histogram(g, hist, maxdeg + 2) == GRAPH_OK);
    size_t sum = 0;
    for (size_t d = 0; d <= maxdeg + 1; d++) sum += hist[d];
    EXPECT(sum == n && hist[0] == isolated);
    free(hist);

    /* components: each silent user is its own component, the rest form the network */
    int *comp = xmalloc(n * sizeof *comp);
    Uf uf = uf_make(n);
    size_t expect = n;
    for (size_t u = 0; u < n; u++)
        for (size_t e = c->off[u]; e < c->off[u + 1]; e++) if (uf_union(&uf, (int)u, c->adj[e])) expect--;
    const size_t got = graph_connected_components(g, comp);
    EXPECT(got == expect && got >= isolated + 1);
    uf_free(&uf);
    free(comp);

    /* reach from the biggest hub = size of its component; a silent user reaches only itself */
    int *level = xmalloc(n * sizeof *level), *order = xmalloc(n * sizeof *order), *parent = xmalloc(n * sizeof *parent);
    const size_t reach = o_bfs(c, (int)maxdeg_v, level);
    EXPECT(graph_bfs(g, (int)maxdeg_v, order, parent) == reach && reach > n / 2);
    size_t silent = 0;
    while (silent < n && c->off[silent + 1] != c->off[silent]) silent++;
    if (silent < n) EXPECT(graph_bfs(g, (int)silent, order, parent) == 1);
    free(level); free(order); free(parent);

    /* PageRank copes with the silent users (dangling vertices): a distribution that sums to 1,
     * identical for all silent users, with the hub far above them */
    double *rank = xmalloc(n * sizeof *rank);
    size_t iters = 0;
    EXPECT(graph_pagerank(g, 0.85, 1e-9, 200, rank, &iters) == GRAPH_OK);
    double s = 0;
    for (size_t v = 0; v < n; v++) s += rank[v];
    EXPECT(fabs(s - 1.0) < 1e-6);
    bool same = true;
    for (size_t v = 0; v < n; v++)
        if (c->off[v + 1] == c->off[v] && fabs(rank[v] - rank[silent]) > 1e-12) { same = false; break; }
    EXPECT(same && rank[maxdeg_v] > 50 * rank[silent]);
    free(rank);

    /* communities and cores never fail on silent users */
    int *label = xmalloc(n * sizeof *label), *core = xmalloc(n * sizeof *core);
    size_t ncomm = 0;
    EXPECT(graph_label_propagation(g, label, &ncomm) == GRAPH_OK && ncomm >= isolated);
    EXPECT(graph_core_numbers(g, core) == GRAPH_OK);
    bool zero_core = true;
    for (size_t v = 0; v < n; v++) if (c->off[v + 1] == c->off[v] && core[v] != 0) zero_core = false;
    EXPECT(zero_core);
    free(label); free(core);

    long long tri = 0;
    EXPECT(graph_triangle_count(g, &tri) == GRAPH_OK && tri >= 0);
}


/* The per-user queries a social app runs: they must stay fast however many users exist. */
static void c_social_queries(void)
{
    Graph *g = get(K_SOCIAL);
    const Csr *c = csr(K_SOCIAL);
    const size_t n = c->n;
    int *dist = xmalloc(n * sizeof *dist);
    int *users = xmalloc(n * sizeof *users), *hops = xmalloc(n * sizeof *hops);
    const double t0 = now();
    size_t queries = 0;
    for (int q = 0; q < 40; q++) {
        int u = (int)rnd_below(n);
        if (c->off[u + 1] - c->off[u] > 2000) continue;     /* skip the biggest hubs: their 2-hop ball is the graph */
        (void)o_bfs(c, u, dist);

        /* everyone within 2 hops, nearest first, with the right distances */
        size_t count = 0, want = 0;
        EXPECT(graph_social_within_hops(g, u, 2, n, users, hops, &count, NULL) == GRAPH_OK);
        for (size_t v = 0; v < n; v++) if ((int)v != u && dist[v] >= 1 && dist[v] <= 2) want++;
        bool ok = count == want;
        for (size_t i = 0; ok && i < count; i++) if (hops[i] != dist[users[i]] || (i && hops[i - 1] > hops[i])) ok = false;
        EXPECT(ok);

        /* suggestions: every one is a real friend-of-friend with the stated number of mutual friends */
        GraphSuggestion sug[5];
        size_t sc = 0;
        EXPECT(graph_social_suggest_friends(g, u, 5, sug, &sc) == GRAPH_OK);
        ok = true;
        for (size_t i = 0; i < sc; i++) {
            size_t mutual = 0;
            for (size_t e = c->off[u]; e < c->off[u + 1]; e++)
                if (graph_has_edge(g, c->adj[e], sug[i].user, NULL)) mutual++;
            if (mutual != sug[i].mutual || dist[sug[i].user] != 2 || (i && sug[i - 1].mutual < sug[i].mutual)) ok = false;
        }
        EXPECT(ok);

        /* how are we connected: the chain is as short as breadth-first search says */
        const int t = (int)rnd_below(n);
        GraphPath p = {0};
        const GraphStatus st = graph_social_path(g, u, t, SIZE_MAX, &p);
        if (dist[t] < 0) {
            EXPECT(st == GRAPH_ERR_NOT_FOUND);
        } else {
            EXPECT(st == GRAPH_OK && p.length == (size_t)dist[t] + 1);
            ok = st == GRAPH_OK;
            for (size_t i = 0; ok && i + 1 < p.length; i++) if (!graph_has_edge(g, p.vertices[i], p.vertices[i + 1], NULL)) ok = false;
            EXPECT(ok);
            graph_path_free(&p);
        }
        queries++;
    }
    EXPECT(queries >= 10);
    /* the oracle BFS dominates this time; the library queries alone are far faster than the
     * whole-graph algorithms (checked loosely so a slow machine does not fail it) */
    EXPECT(now() - t0 < 120.0);
    free(dist); free(users); free(hops);
}

/* The memory budget at scale: a graph this size is small as a list and impossible as a matrix.
 * Converting or loading it into the matrix must be refused with TOO_LARGE, not crash. */
static void c_budget_refusal(void)
{
    if (g_n < 100000) {
        return;     /* below this the matrix still fits in 4 GiB, so there is nothing to refuse */
    }
    Graph *g = get(K_SOCIAL);
    EXPECT(!graph_convert(g, GRAPH_MATRIX));
    EXPECT(graph_last_error() == GRAPH_ERR_TOO_LARGE);
    EXPECT(graph_representation(g) == GRAPH_LIST);       /* untouched and still usable */
    EXPECT(graph_degree(g, 0) < graph_vertex_count(g));

    FILE *f = tmp();
    EXPECT(graph_save_binary(g, f) == GRAPH_OK);
    rewind(f);
    Graph *back = NULL;
    EXPECT(graph_load_binary(f, GRAPH_MATRIX, &back) != GRAPH_OK && back == NULL);
    EXPECT(graph_last_error() == GRAPH_ERR_TOO_LARGE);
    fclose(f);
}

/* ---- graph operations, serialisation, generators, mutation ---- */

static void c_builders(void)
{
    /* transpose of a directed path is the reversed path */
    Graph *p = get(K_DPATH);
    const size_t n = size_of(K_DPATH);
    Graph *t = NULL;
    EXPECT(graph_transpose(p, &t) == GRAPH_OK && t != NULL);
    bool ok = graph_edge_count(t) == n - 1;
    for (size_t i = 0; ok && i + 1 < n; i++) if (!graph_has_edge(t, (int)i + 1, (int)i, NULL) || graph_has_edge(t, (int)i, (int)i + 1, NULL)) ok = false;
    EXPECT(ok);
    graph_destroy(t);

    /* copy, union and intersection of a graph with itself */
    Graph *r = get(K_RAND);
    Graph *cp = graph_copy(r);
    EXPECT(cp != NULL && same_graph(r, cp));
    Graph *u = NULL, *i = NULL;
    EXPECT(graph_union(r, cp, &u) == GRAPH_OK && same_graph(r, u));
    EXPECT(graph_intersection(r, cp, &i) == GRAPH_OK && same_graph(r, i));
    graph_destroy(u); graph_destroy(i); graph_destroy(cp);

    /* induced subgraph on the first half of a path is a shorter path */
    Graph *pa = get(K_PATH);
    const size_t pn = size_of(K_PATH), half = pn / 2;
    int *keep = xmalloc(half * sizeof *keep);
    for (size_t v = 0; v < half; v++) keep[v] = (int)v;
    Graph *sub = NULL;
    EXPECT(graph_induced_subgraph(pa, keep, half, &sub) == GRAPH_OK && sub != NULL);
    EXPECT(graph_vertex_count(sub) == half && graph_edge_count(sub) == half - 1);
    graph_destroy(sub);
    free(keep);
}

static void c_quadratic_ops(void)
{
    /* closure of a path DAG has n(n-1)/2 arcs and its reduction is the path again */
    Graph *p = get(K_DPATH);
    const size_t n = size_of(K_DPATH);
    Graph *cl = NULL, *red = NULL, *co = NULL;
    EXPECT(graph_transitive_closure(p, &cl) == GRAPH_OK && graph_edge_count(cl) == n * (n - 1) / 2);
    EXPECT(graph_transitive_reduction(cl, &red) == GRAPH_OK && same_graph(red, p));
    graph_destroy(cl); graph_destroy(red);
    /* complement of a path: all other pairs */
    Graph *u = get(K_PATH);
    EXPECT(graph_complement(u, &co) == GRAPH_OK && graph_edge_count(co) == n * (n - 1) / 2 - (n - 1));
    graph_destroy(co);
}

static void c_io_roundtrip(void)
{
    /* Load into the OTHER representation to prove the snapshot is representation-independent --
     * but only where that is affordable: the matrix needs ~V^2/8 bytes, so for a large graph
     * the snapshot is loaded back into the same (list) representation instead. */
    const GraphRep other_rep = (g_rep == GRAPH_LIST && g_n <= 10000) ? GRAPH_MATRIX
                               : g_rep == GRAPH_MATRIX                ? GRAPH_LIST
                                                                      : GRAPH_LIST;
    const Kind kinds[] = {K_RAND, K_RANDD};
    for (int t = 0; t < 2; t++) {
        Graph *g = get(kinds[t]);
        FILE *f = tmp();
        EXPECT(graph_save_binary(g, f) == GRAPH_OK);
        rewind(f);
        Graph *back = NULL;
        EXPECT(graph_load_binary(f, other_rep, &back) == GRAPH_OK && back != NULL);
        if (back != NULL) {
            EXPECT(same_graph(g, back));   /* exact, whatever representation it was loaded into */
            graph_destroy(back);
        }
        fclose(f);

        f = tmp();
        EXPECT(graph_write_edge_list(g, f) == GRAPH_OK);
        rewind(f);
        EXPECT(graph_read_edge_list(f, graph_is_directed(g), g_rep, &back) == GRAPH_OK && back != NULL);
        if (back != NULL) {
            EXPECT(graph_edge_count(back) == graph_edge_count(g) && graph_vertex_count(back) == graph_vertex_count(g));
            graph_destroy(back);
        }
        fclose(f);
    }
}

static void c_io_export(void)
{
    Graph *g = get(K_RAND);
    GraphStatus (*const writers[])(Graph *, FILE *) = {graph_write_graphml, graph_write_gml, graph_write_json};
    for (int i = 0; i < 3; i++) {
        FILE *f = tmp();
        EXPECT(writers[i](g, f) == GRAPH_OK);
        EXPECT(ftell(f) > (long)size_of(K_RAND));   /* wrote at least a few bytes per vertex */
        fclose(f);
    }
    FILE *f = tmp();
    EXPECT(graph_write_dot(g, f, "g", NULL, true) == GRAPH_OK && ftell(f) > (long)size_of(K_RAND));
    fclose(f);
}

static void c_generators(void)
{
    const size_t n = g_n;
    const GraphRep rep = g_rep;
    Graph *g;
    g = graph_gen_path(rep, n, false);
    EXPECT(g && graph_vertex_count(g) == n && graph_edge_count(g) == n - 1 && graph_is_tree(g));
    graph_destroy(g);
    g = graph_gen_cycle(rep, n, true);
    EXPECT(g && graph_edge_count(g) == n && graph_is_strongly_connected(g));
    graph_destroy(g);
    g = graph_gen_star(rep, n);
    EXPECT(g && graph_edge_count(g) == n - 1 && graph_degree(g, 0) == n - 1);
    graph_destroy(g);
    const size_t s = grid_side(n);
    g = graph_gen_grid(rep, s, s);
    EXPECT(g && graph_vertex_count(g) == s * s && graph_edge_count(g) == 2 * s * (s - 1));
    graph_destroy(g);
    g = graph_gen_random_tree(rep, n, 7);
    EXPECT(g && graph_edge_count(g) == n - 1 && graph_is_tree(g));
    graph_destroy(g);
    g = graph_gen_barabasi_albert(rep, n, 3, 7);
    EXPECT(g && graph_vertex_count(g) == n && graph_is_connected(g) && graph_edge_count(g) >= n - 1);
    graph_destroy(g);
    if (n <= 20000) {
        const double p = 8.0 / (double)n;
        g = graph_gen_erdos_renyi(rep, n, p, 7, false, 1, 5);
        const double expect = p * (double)n * ((double)n - 1) / 2.0;
        EXPECT(g && fabs((double)graph_edge_count(g) - expect) < 6 * sqrt(expect));   /* 6 sigma */
        graph_destroy(g);
        g = graph_gen_random_dag(rep, n, p, 7, 1, 5);
        EXPECT(g && graph_is_dag(g) && graph_edge_count(g) > 0);
        graph_destroy(g);
    }
    if (n <= 1000) {
        g = graph_gen_complete(rep, n, false);
        EXPECT(g && graph_edge_count(g) == n * (n - 1) / 2);
        graph_destroy(g);
    }
}

static void c_mutation(void)
{
    Graph *g = graph_copy(get(K_RAND));
    EXPECT(g != NULL);
    if (g == NULL) return;
    graph_set_adaptive(g, false);
    const size_t n = graph_id_bound(g);
    size_t removed = 0;
    for (size_t v = 0; v < n; v += 7) if (graph_remove_vertex(g, (int)v)) removed++;
    EXPECT(graph_vertex_count(g) == n - removed);
    const size_t edges = graph_edge_count(g);

    /* a surviving edge is still there, a removed vertex has none */
    const Csr c = csr_of(g);
    bool ok = true;
    for (size_t u = 0; u < c.n; u++) {
        if (u % 7 == 0) { if (c.off[u + 1] != c.off[u] || graph_is_alive(g, (int)u)) ok = false; continue; }
        for (size_t e = c.off[u]; e < c.off[u + 1]; e++) if (c.adj[e] % 7 == 0) ok = false;
    }
    EXPECT(ok);
    Csr tmpc = c; csr_free(&tmpc);

    int *remap = xmalloc(n * sizeof *remap);
    EXPECT(graph_compact(g, remap));
    EXPECT(graph_id_bound(g) == n - removed && graph_vertex_count(g) == n - removed && graph_edge_count(g) == edges);
    size_t mapped = 0;
    for (size_t v = 0; v < n; v++) if (remap[v] >= 0) mapped++;
    EXPECT(mapped == n - removed);
    free(remap);

    /* removing and re-adding an edge */
    const int a = 1, b = 2;
    const bool had = graph_has_edge(g, a, b, NULL);
    if (!had) EXPECT(graph_add_edge(g, a, b, 3.5));
    double w = 0;
    EXPECT(graph_has_edge(g, a, b, &w));
    EXPECT(graph_remove_edge(g, a, b) && !graph_has_edge(g, a, b, NULL));
    graph_destroy(g);
}

static void c_convert(void)
{
    /* list -> matrix -> list keeps everything (only run where the matrix fits) */
    Graph *g = graph_copy(get(K_RANDD));
    EXPECT(g != NULL);
    if (g == NULL) return;
    Graph *orig = graph_copy(g);
    EXPECT(graph_convert(g, GRAPH_MATRIX) && graph_representation(g) == GRAPH_MATRIX);
    EXPECT(same_graph(orig, g));
    EXPECT(graph_convert(g, GRAPH_LIST) && graph_representation(g) == GRAPH_LIST);
    EXPECT(same_graph(orig, g));
    graph_destroy(g); graph_destroy(orig);
}

/* ---- advanced algorithms ---- */

static void c_adv_tree(void)
{
    Graph *g = get(K_TREE);
    const Csr *c = csr(K_TREE);
    const size_t n = c->n;
    int *level = xmalloc(n * sizeof *level);
    (void)o_bfs(c, 0, level);
    int *out = xmalloc(n * sizeof *out), *pos = xmalloc(n * sizeof *pos);
    size_t cnt = 0;

    EXPECT(graph_adv_tree_traversal(g, 0, GRAPH_ADV_TREE_PREORDER, out, &cnt) == GRAPH_OK && cnt == n);
    for (size_t i = 0; i < n; i++) pos[out[i]] = (int)i;
    bool ok = out[0] == 0;
    for (size_t u = 0; ok && u < n; u++)           /* a parent comes before each child */
        for (size_t e = c->off[u]; e < c->off[u + 1]; e++)
            if (level[c->adj[e]] == level[u] + 1 && pos[c->adj[e]] < pos[u]) { ok = false; break; }
    EXPECT(ok);

    EXPECT(graph_adv_tree_traversal(g, 0, GRAPH_ADV_TREE_POSTORDER, out, &cnt) == GRAPH_OK && cnt == n);
    for (size_t i = 0; i < n; i++) pos[out[i]] = (int)i;
    ok = out[n - 1] == 0;
    for (size_t u = 0; ok && u < n; u++)           /* a child comes before its parent */
        for (size_t e = c->off[u]; e < c->off[u + 1]; e++)
            if (level[c->adj[e]] == level[u] + 1 && pos[c->adj[e]] > pos[u]) { ok = false; break; }
    EXPECT(ok);

    EXPECT(graph_adv_tree_traversal(g, 0, GRAPH_ADV_TREE_LEVELORDER, out, &cnt) == GRAPH_OK && cnt == n);
    ok = true;
    for (size_t i = 1; i < n; i++) if (level[out[i - 1]] > level[out[i]]) { ok = false; break; }
    EXPECT(ok);
    free(level); free(out); free(pos);
}

static void c_adv_refine_cycle(void)
{
    size_t nc = 0, rounds = 0;
    {
        const size_t n = size_of(K_CYCLE);
        int *color = xmalloc(n * sizeof *color);
        EXPECT(graph_adv_color_refinement(get(K_CYCLE), color, &nc, &rounds) == GRAPH_OK && nc == 1);
        free(color);
    }
    {
        const size_t n = size_of(K_STAR);
        int *color = xmalloc(n * sizeof *color);
        EXPECT(graph_adv_color_refinement(get(K_STAR), color, &nc, &rounds) == GRAPH_OK && nc == 2);
        free(color);
    }
}

static void c_adv_refine_path(void)
{
    size_t nc = 0, rounds = 0;
    {   /* a path has ceil(n/2) classes: v and n-1-v are interchangeable */
        const size_t n = size_of(K_PATH);
        int *color = xmalloc(n * sizeof *color);
        EXPECT(graph_adv_color_refinement(get(K_PATH), color, &nc, &rounds) == GRAPH_OK && nc == (n + 1) / 2);
        bool ok = true;
        for (size_t v = 0; v < n; v++) if (color[v] != color[n - 1 - v]) { ok = false; break; }
        EXPECT(ok);
        free(color);
    }
}

static void run_iso(Kind kind)
{
    const size_t n = size_of(kind);
    Graph *a = get(kind);
    const bool directed = graph_is_directed(a);
    const Csr *c = csr(kind);
    int *perm = xmalloc(n * sizeof *perm);
    for (size_t i = 0; i < n; i++) perm[i] = (int)i;
    for (size_t i = n - 1; i > 0; i--) { size_t j = rnd_below(i + 1); int t = perm[i]; perm[i] = perm[j]; perm[j] = t; }
    Graph *b = new_graph(n, directed);
    for (size_t u = 0; u < n; u++)
        for (size_t e = c->off[u]; e < c->off[u + 1]; e++)
            if (directed || (int)u < c->adj[e]) (void)graph_add_edge(b, perm[u], perm[c->adj[e]], c->w[e]);
    bool iso = false;
    int *map = xmalloc(n * sizeof *map);
    EXPECT(graph_adv_graph_isomorphism(a, b, true, &iso, map) == GRAPH_OK && iso);
    bool ok = iso;
    for (size_t u = 0; ok && u < n; u++)           /* the mapping really is an isomorphism */
        for (size_t e = c->off[u]; e < c->off[u + 1]; e++) {
            double w;
            if (!graph_has_edge(b, map[u], map[c->adj[e]], &w) || w != c->w[e]) { ok = false; break; }
        }
    EXPECT(ok);

    /* change one weight: structure the same, weights differ */
    int u0 = perm[0], v0 = -1;
    int nb[64];
    double w0[64];
    const size_t d = graph_neighbors(b, u0, nb, w0, 64);
    if (d > 0) {
        v0 = nb[0];
        graph_remove_edge(b, u0, v0);
        graph_add_edge(b, u0, v0, w0[0] + 1.0);
        EXPECT(graph_adv_graph_isomorphism(a, b, true, &iso, NULL) == GRAPH_OK && !iso);
        EXPECT(graph_adv_graph_isomorphism(a, b, false, &iso, NULL) == GRAPH_OK && iso);
    }
    /* move an edge: no longer isomorphic */
    Graph *e = graph_copy(a);
    for (size_t t = 0; t < 1000 && e != NULL; t++) {   /* find a non-adjacent pair and rewire */
        int x = (int)rnd_below(n), y = (int)rnd_below(n);
        if (x != y && !graph_has_edge(e, x, y, NULL)) {
            int nb2[64]; double w2[64];
            size_t d2 = graph_neighbors(e, x, nb2, w2, 64);
            if (d2 > 1 && d2 < 64) { graph_remove_edge(e, x, nb2[0]); graph_add_edge(e, x, y, 1.0); break; }
        }
    }
    EXPECT(graph_adv_graph_isomorphism(a, e, false, &iso, NULL) == GRAPH_OK && !iso);
    graph_destroy(e); graph_destroy(b);
    free(perm); free(map);
}

/* A random digraph has almost no symmetry: colour refinement separates nearly every vertex,
 * so the search is short and this scales.  The tree-plus-edges graph is full of interchangeable
 * leaves, which is the hard case for exact search, so it only runs on small sizes. */
static void c_adv_iso_asym(void) { run_iso(K_RANDD); }
static void c_adv_iso_sym(void)  { run_iso(K_RAND); }

static void c_adv_subiso(void)
{
    Graph *tri = new_graph(3, false);
    edge(tri, 0, 1, 1); edge(tri, 1, 2, 1); edge(tri, 0, 2, 1);
    bool found = false;
    int map[3];
    EXPECT(graph_adv_subgraph_isomorphism(tri, get(K_CHAIN5), true, false, map, &found) == GRAPH_OK && found);
    EXPECT(graph_has_edge(get(K_CHAIN5), map[0], map[1], NULL) && graph_has_edge(get(K_CHAIN5), map[1], map[2], NULL) &&
           graph_has_edge(get(K_CHAIN5), map[0], map[2], NULL));
    /* proving that a pattern is ABSENT scans every candidate, which is quadratic in the target
     * for this implementation (measured: 0.8 s at 1 000, over a minute at 10 000) */
    if (g_n <= 2000) {
        EXPECT(graph_adv_subgraph_isomorphism(tri, get(K_TREE), false, false, NULL, &found) == GRAPH_OK && !found);
        EXPECT(graph_adv_subgraph_isomorphism(tri, get(K_GRID), false, false, NULL, &found) == GRAPH_OK && !found);
    }
    graph_destroy(tri);
}

static void c_adv_paths(void)
{
    /* widest path = maximum bottleneck, checked against a modified Dijkstra */
    {
        Graph *g = get(K_RAND);
        const Csr *c = csr(K_RAND);
        double *want = xmalloc(c->n * sizeof *want), *got = xmalloc(c->n * sizeof *got);
        o_widest(c, 0, want);
        EXPECT(graph_adv_widest_path(g, 0, got, NULL) == GRAPH_OK);
        bool ok = true;
        for (size_t v = 1; v < c->n; v++) if (!near(got[v], want[v])) { ok = false; break; }
        EXPECT(ok);
        free(want); free(got);
    }
    /* flood on a grid: the corner message reaches the far corner after 2(side-1) rounds */
    {
        const size_t n = size_of(K_GRID);
        const int side = (int)grid_side(n);
        int *round = xmalloc(n * sizeof *round), *parent = xmalloc(n * sizeof *parent);
        size_t msgs = 0, reached = 0;
        EXPECT(graph_adv_flood(get(K_GRID), 0, (size_t)n, round, parent, &msgs, &reached) == GRAPH_OK);
        EXPECT(reached == n && round[n - 1] == 2 * (side - 1));
        free(round); free(parent);
    }
    /* jump point search on an open grid (the graph's vertices are the cells): octile distance */
    {
        const size_t side = grid_side(g_n);
        Graph *cells = new_graph(side * side, false);
        if (g_rep == GRAPH_MATRIX) EXPECT(graph_convert(cells, GRAPH_MATRIX));
        GraphPath p = {0};
        EXPECT(graph_adv_jump_point_search(cells, side, 0, (int)(side * side) - 1, &p) == GRAPH_OK);
        EXPECT(near(p.cost, (double)(side - 1) * sqrt(2.0)));
        graph_path_free(&p);
        graph_destroy(cells);
    }
}

static void c_adv_gnn(void)
{
    /* GNN on a regular graph with constant features and identity weights: output == input */
    {
        const size_t n = size_of(K_CYCLE);
        double *feat = xmalloc(n * 2 * sizeof *feat), *out = xmalloc(n * 2 * sizeof *out);
        for (size_t i = 0; i < n * 2; i++) feat[i] = (i % 2) ? 3.0 : 1.5;
        const double w[4] = {1, 0, 0, 1};
        const GraphAdvGnnLayer layer = {2, 2, w, NULL, GRAPH_ADV_ACT_NONE};
        EXPECT(graph_adv_gnn_forward(get(K_CYCLE), feat, 2, &layer, 1, out) == GRAPH_OK);
        bool ok = true;
        for (size_t i = 0; i < n * 2; i++) if (!near(out[i], feat[i])) { ok = false; break; }
        EXPECT(ok);
        free(feat); free(out);
    }
}

static void c_adv_circle(void)
{
    /* circle layout: every vertex at distance r from the centre */
    {
        const size_t n = size_of(K_CYCLE);
        double *x = xmalloc(n * sizeof *x), *y = xmalloc(n * sizeof *y);
        EXPECT(graph_adv_circle_layout(get(K_CYCLE), 5.0, x, y) == GRAPH_OK);
        bool ok = true;
        for (size_t i = 0; i < n; i++) if (fabs(hypot(x[i], y[i]) - 5.0) > 1e-6) { ok = false; break; }
        EXPECT(ok);
        free(x); free(y);
    }
}

static void c_adv_lexbfs(void)
{
    /* lex-BFS gives a permutation */
    {
        const size_t n = size_of(K_RAND);
        int *order = xmalloc(n * sizeof *order);
        char *seen = xcalloc(n, 1);
        EXPECT(graph_adv_lex_bfs(get(K_RAND), 0, order) == GRAPH_OK);
        bool ok = order[0] == 0;
        for (size_t i = 0; i < n; i++) { if (order[i] < 0 || (size_t)order[i] >= n || seen[order[i]]) { ok = false; break; } seen[order[i]] = 1; }
        EXPECT(ok);
        free(order); free(seen);
    }
}

static void c_adv_colorcoding(void)
{
    /* colour coding finds a 4-vertex path in a path graph */
    {
        bool found = false;
        int path[8];
        EXPECT(graph_adv_color_coding_path(get(K_PATH), 4, 400, 11, &found, path) == GRAPH_OK && found);
        bool ok = found;
        for (int i = 0; ok && i + 1 < 4; i++) if (!graph_has_edge(get(K_PATH), path[i], path[i + 1], NULL)) ok = false;
        EXPECT(ok);
    }
}

static void c_adv_dynconn(void)
{
    Graph *g = graph_copy(get(K_PATH));
    EXPECT(g != NULL);
    if (g == NULL) return;
    graph_set_adaptive(g, false);
    const int n = (int)graph_vertex_count(g), mid = n / 2;
    GraphAdvDynConn *dc = graph_adv_dynconn_create(g);
    EXPECT(dc != NULL);
    if (dc == NULL) { graph_destroy(g); return; }
    bool con = false;
    size_t comps = 0;
    EXPECT(graph_adv_dynconn_connected(dc, 0, n - 1, &con) == GRAPH_OK && con);
    EXPECT(graph_adv_dynconn_components(dc, &comps) == GRAPH_OK && comps == 1);
    EXPECT(graph_remove_edge(g, mid, mid + 1));
    EXPECT(graph_adv_dynconn_connected(dc, 0, n - 1, &con) == GRAPH_OK && !con);
    EXPECT(graph_adv_dynconn_connected(dc, 0, mid, &con) == GRAPH_OK && con);
    EXPECT(graph_adv_dynconn_components(dc, &comps) == GRAPH_OK && comps == 2);
    EXPECT(graph_add_edge(g, mid, mid + 1, 1.0));
    EXPECT(graph_adv_dynconn_connected(dc, 0, n - 1, &con) == GRAPH_OK && con);
    graph_adv_dynconn_destroy(dc);
    graph_destroy(g);
}

/* ---- algorithms that are meant for small and medium graphs only ---- */

static void c_distance_metrics(void)
{
    Graph *g = get(K_PATH);
    const size_t n = size_of(K_PATH);
    int *ecc = xmalloc(n * sizeof *ecc);
    EXPECT(graph_eccentricity(g, ecc) == GRAPH_OK);
    bool ok = true;
    for (size_t v = 0; v < n; v++) if ((size_t)ecc[v] != (v > n - 1 - v ? v : n - 1 - v)) ok = false;
    EXPECT(ok);
    int diam = 0, rad = 0;
    EXPECT(graph_diameter_radius(g, &diam, &rad) == GRAPH_OK && (size_t)diam == n - 1 && (size_t)rad == n / 2);

    double *close = xmalloc(n * sizeof *close), *btw = xmalloc(n * sizeof *btw), *harm = xmalloc(n * sizeof *harm);
    EXPECT(graph_closeness_centrality(g, close) == GRAPH_OK);
    EXPECT(graph_betweenness_centrality(g, btw) == GRAPH_OK);
    EXPECT(graph_adv_harmonic_centrality(g, harm) == GRAPH_OK);
    double *hn = xmalloc((n + 1) * sizeof *hn);   /* harmonic numbers */
    hn[0] = 0;
    for (size_t k = 1; k <= n; k++) hn[k] = hn[k - 1] + 1.0 / (double)k;
    bool okc = true, okb = true, okh = true;
    for (size_t v = 0; v < n; v++) {
        const double sum = (double)v * (double)(v + 1) / 2.0 + (double)(n - 1 - v) * (double)(n - v) / 2.0;
        if (!near(close[v], (double)(n - 1) / sum)) okc = false;
        if (!near(btw[v], (double)v * (double)(n - 1 - v))) okb = false;
        if (!near(harm[v], (hn[v] + hn[n - 1 - v]) / (double)(n - 1))) okh = false;
    }
    EXPECT(okc); EXPECT(okb); EXPECT(okh);
    free(ecc); free(close); free(btw); free(harm); free(hn);
}

static void c_all_pairs(void)
{
    Graph *g = get(K_RANDD);
    const Csr *c = csr(K_RANDD);
    const size_t n = c->n;
    double *fw = xmalloc(n * n * sizeof *fw), *jo = xmalloc(n * n * sizeof *jo), *aj = xmalloc(n * n * sizeof *aj);
    int *fnext = xmalloc(n * n * sizeof *fnext), *jnext = xmalloc(n * n * sizeof *jnext), *apred = xmalloc(n * n * sizeof *apred);
    EXPECT(graph_floyd_warshall(g, fw, fnext) == GRAPH_OK);
    EXPECT(graph_johnson(g, jo, jnext) == GRAPH_OK);
    EXPECT(graph_adv_johnson(g, aj, apred) == GRAPH_OK);
    double *want = xmalloc(n * sizeof *want);
    bool okf = true, okj = true, oka = true;
    for (size_t s = 0; s < n; s++) {
        const int src = (int)s;
        o_dijkstra(c, &src, 1, want);
        for (size_t t = 0; t < n; t++) {
            if (!near(fw[s * n + t], want[t])) okf = false;
            if (!near(jo[s * n + t], want[t])) okj = false;
            if (!near(aj[s * n + t], want[t])) oka = false;
        }
    }
    EXPECT(okf); EXPECT(okj); EXPECT(oka);
    /* paths rebuilt from `next` and from `pred` have the right cost */
    bool okp = true;
    for (int trial = 0; trial < 200; trial++) {
        const int s = (int)rnd_below(n), t = (int)rnd_below(n);
        if (s == t || isinf(fw[(size_t)s * n + t])) continue;
        int path[600];
        const size_t len = graph_path_from_next(fnext, n, s, t, path, 600);
        double cost = 0;
        for (size_t i = 0; i + 1 < len; i++) {
            double w;
            if (!graph_has_edge(g, path[i], path[i + 1], &w)) { okp = false; break; }
            cost += w;
        }
        if (len < 2 || !near(cost, fw[(size_t)s * n + t])) okp = false;
        GraphPath ap = {0};
        if (graph_adv_apsp_path(aj, apred, n, s, t, &ap) != GRAPH_OK || !near(ap.cost, fw[(size_t)s * n + t])) okp = false;
        graph_path_free(&ap);
    }
    EXPECT(okp);
    free(fw); free(jo); free(aj); free(fnext); free(jnext); free(apred); free(want);
}

static void c_k_shortest(void)
{
    Graph *g = get(K_RAND);
    const Csr *c = csr(K_RAND);
    const int src = 0, dst = (int)c->n - 1;
    double *want = xmalloc(c->n * sizeof *want);
    o_dijkstra(c, &src, 1, want);
    GraphPath paths[4] = {{0}};
    size_t found = 0;
    EXPECT(graph_k_shortest_paths(g, src, dst, 4, paths, &found) == GRAPH_OK && found >= 1);
    EXPECT(found > 0 && near(paths[0].cost, want[dst]));
    bool ok = true;
    for (size_t i = 0; i < found; i++) {
        if (i > 0 && paths[i].cost < paths[i - 1].cost - 1e-9) ok = false;   /* non-decreasing */
        double sum = 0;
        for (size_t j = 0; j + 1 < paths[i].length; j++) {
            double w;
            if (!graph_has_edge(g, paths[i].vertices[j], paths[i].vertices[j + 1], &w)) ok = false; else sum += w;
        }
        if (!near(sum, paths[i].cost)) ok = false;
        for (size_t j = 0; j < paths[i].length && ok; j++)      /* loopless */
            for (size_t l = j + 1; l < paths[i].length; l++) if (paths[i].vertices[j] == paths[i].vertices[l]) ok = false;
        if (i > 0 && paths[i].length == paths[i - 1].length &&
            memcmp(paths[i].vertices, paths[i - 1].vertices, paths[i].length * sizeof(int)) == 0) ok = false;   /* distinct */
        graph_path_free(&paths[i]);
    }
    EXPECT(ok);
    free(want);
}

static void c_exact_small(void)
{
    /* Hamiltonian search, exact cliques, the Hungarian method: sizes fixed by the algorithm */
    Graph *c16 = new_graph(16, false), *p16 = new_graph(16, false), *s16 = new_graph(16, false);
    for (int i = 0; i < 16; i++) { edge(c16, i, (i + 1) % 16, 1); if (i < 15) edge(p16, i, i + 1, 1); if (i > 0) edge(s16, 0, i, 1); }
    int path[16];
    EXPECT(graph_hamiltonian(c16, true, path) == GRAPH_OK);
    bool ok = true;
    for (int i = 0; i < 16; i++) if (!graph_has_edge(c16, path[i], path[(i + 1) % 16], NULL)) ok = false;
    EXPECT(ok);
    EXPECT(graph_hamiltonian(p16, false, path) == GRAPH_OK);
    EXPECT(graph_hamiltonian(p16, true, path) == GRAPH_ERR_NOT_FOUND);
    EXPECT(graph_hamiltonian(s16, false, path) == GRAPH_ERR_NOT_FOUND);
    graph_destroy(c16); graph_destroy(p16); graph_destroy(s16);

    /* cost[i][j] = i * j is minimised by pairing small with large: sum i*(n-1-i) */
    enum { N = 300 };
    double *cost = xmalloc((size_t)N * N * sizeof *cost);
    int *assign = xmalloc(N * sizeof *assign);
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) cost[(size_t)i * N + j] = (double)i * j;
    double total = 0, want = 0;
    EXPECT(graph_hungarian(cost, N, N, assign, &total) == GRAPH_OK);
    for (int i = 0; i < N; i++) want += (double)i * (N - 1 - i);
    EXPECT(near(total, want));
    free(cost); free(assign);
}

static void c_adv_small(void)
{
    /* Coffman-Graham: a valid schedule on 3 processors */
    {
        Graph *g = get(K_DAG);
        const Csr *c = csr(K_DAG);
        const size_t n = c->n;
        int *label = xmalloc(n * sizeof *label), *slot = xmalloc(n * sizeof *slot);
        size_t span = 0;
        EXPECT(graph_adv_coffman_graham(g, 3, label, slot, &span) == GRAPH_OK);
        EXPECT(span >= (n + 2) / 3);
        bool ok = true;
        size_t *load = xcalloc(span + 1, sizeof *load);
        for (size_t u = 0; u < n; u++) {
            if (slot[u] < 0 || (size_t)slot[u] >= span || ++load[slot[u]] > 3) ok = false;
            for (size_t e = c->off[u]; ok && e < c->off[u + 1]; e++) if (slot[u] >= slot[c->adj[e]]) ok = false;
        }
        EXPECT(ok);
        free(label); free(slot); free(load);
    }
    /* crossing numbers: a cycle in order has none; the heuristic's count matches a recount */
    {
        const size_t n = size_of(K_CYCLE);
        int *order = xmalloc(n * sizeof *order);
        for (size_t i = 0; i < n; i++) order[i] = (int)i;
        long long cr = -1;
        EXPECT(graph_adv_crossing_count(get(K_CYCLE), order, n, &cr) == GRAPH_OK && cr == 0);
        long long heur = -1, recount = -2;
        EXPECT(graph_adv_crossing_number_heuristic(get(K_RING5), 0, &heur, order) == GRAPH_OK);
        const size_t rn = size_of(K_RING5);
        EXPECT(graph_adv_crossing_count(get(K_RING5), order, rn, &recount) == GRAPH_OK && heur == recount);
        free(order);
    }
    /* reverse-delete builds a tree as light as Kruskal's */
    {
        Graph *g = get(K_RAND);
        const size_t n = size_of(K_RAND);
        GraphEdge *e1 = xmalloc(n * sizeof *e1), *e2 = xmalloc(n * sizeof *e2);
        size_t c1 = 0, c2 = 0, comps = 0;
        double total = 0;
        EXPECT(graph_adv_reverse_delete(g, e1, &c1, &total, &comps) == GRAPH_OK && comps == 1 && c1 == n - 1);
        EXPECT(near(total, graph_mst_kruskal(g, e2, &c2)));
        free(e1); free(e2);
    }
    /* cycles: the directed cycle is the only simple cycle; with total weight 0 it is a zero cycle */
    {
        const size_t n = size_of(K_DCYCLE);
        GraphAdvCycles cyc = {0};
        EXPECT(graph_adv_rocha_thatte(get(K_DCYCLE), 0, 0, &cyc) == GRAPH_OK && cyc.count == 1);
        EXPECT(cyc.count == 1 && cyc.offsets[1] == n);
        graph_adv_cycles_free(&cyc);

        Graph *z = new_graph(n, true);
        for (size_t i = 0; i + 1 < n; i++) edge(z, i, i + 1, 1.0);
        edge(z, n - 1, 0, -(double)(n - 1));
        GraphAdvCycleKind kind;
        int *cycle = xmalloc(n * sizeof *cycle);
        size_t len = 0;
        EXPECT(graph_adv_zero_weight_cycle(z, 1e-9, &kind, cycle, &len) == GRAPH_OK && kind == GRAPH_ADV_CYCLE_ZERO && len == n);
        free(cycle);
        graph_destroy(z);
    }
    /* force layout stays in the unit square */
    {
        const size_t n = size_of(K_GRID);
        double *x = xmalloc(n * sizeof *x), *y = xmalloc(n * sizeof *y);
        EXPECT(graph_adv_force_layout(get(K_GRID), 30, x, y) == GRAPH_OK);
        bool ok = true;
        for (size_t i = 0; i < n; i++) if (!(x[i] >= 0 && x[i] <= 1 && y[i] >= 0 && y[i] <= 1)) ok = false;
        EXPECT(ok);
        free(x); free(y);
    }
    /* uninformed searches on a path: depth n - 1 */
    {
        const size_t n = size_of(K_PATH);
        GraphPath p = {0};
        EXPECT(graph_adv_iddfs(get(K_PATH), 0, (int)n - 1, SIZE_MAX, &p) == GRAPH_OK && p.length == n && near(p.cost, (double)(n - 1)));
        graph_path_free(&p);
        EXPECT(graph_adv_ida_star(get(K_PATH), 0, (int)n - 1, NULL, NULL, 0, &p) == GRAPH_OK && near(p.cost, (double)(n - 1)));
        graph_path_free(&p);
    }
}

/* ================================================================== */
/*  Registry and driver                                                 */
/* ================================================================== */

typedef struct {
    const char *name;
    void      (*fn)(void);
    size_t      cap_list;     /* largest n the case runs at on the list representation   */
    size_t      cap_matrix;   /* ... on the matrix representation                        */
} Case;

#define M 1000000u
static const Case CASES[] = {
    /* traversal and connectivity */
    {"bfs",                    c_bfs,                M,       100000},
    {"dfs",                    c_dfs,                M,       100000},
    {"connected components",   c_components,         M,       100000},
    {"graph predicates",       c_predicates,         M,       100000},
    {"topological/DAG paths",  c_topo,               M,       100000},
    {"strongly connected",     c_scc,                M,       100000},
    {"bridges/cut vertices",   c_blocks,             M,       100000},
    /* shortest paths */
    {"dijkstra",               c_dijkstra,           M,       100000},
    {"a-star",                 c_astar,              M,       100000},
    {"0-1 bfs",                c_01bfs,              M,       100000},
    {"dijkstra multi-source",  c_multi,              M,       100000},
    {"bellman-ford",           c_bellman,            M,       10000},
    {"k shortest paths",       c_k_shortest,         M,       10000},
    {"all-pairs (fw/johnson)", c_all_pairs,          500,     500},
    /* trees, flow, matching */
    {"mst kruskal/prim",       c_mst,                M,       100000},
    {"max flow / min cut",     c_flow,               M,       100000},
    {"bipartite matching",     c_matching,           M,       100000},
    /* combinatorial */
    {"colour: greedy",         c_color_greedy,       M,       100000},
    {"colour: welsh-powell",   c_color_wp,           M,       100000},
    {"colour: dsatur",         c_color_dsatur,       M,       100000},
    {"cover / independent",    c_cover_independent,  M,       100000},
    {"maximal cliques",        c_maximal_cliques,    M,       100000},
    {"max clique",             c_max_clique,         M,       100000},
    {"triangles/stats",        c_triangles_stats,    M,       100000},
    {"euler trail",            c_euler,              M,       100000},
    {"hamiltonian/hungarian",  c_exact_small,        1000,    1000},
    /* ranking and communities */
    {"pagerank",               c_pagerank,           M,       100000},
    {"eigenvector vs adv",     c_eigenvector_compare, 10000,  10000},
    {"eigenvector symmetry",   c_eigenvector_symmetry, 100000, 100000},
    {"hits",                   c_hits,               M,       100000},
    {"core numbers",           c_cores,              M,       100000},
    {"louvain/label prop.",    c_communities,        M,       100000},
    {"distance metrics",       c_distance_metrics,   2000,    2000},
    /* social network */
    {"social network",         c_social,             M,       10000},
    {"social queries",         c_social_queries,     M,       10000},
    {"memory budget refusal",  c_budget_refusal,     M,       0},
    /* graph operations and I/O */
    {"copy/union/transpose",   c_builders,           M,       100000},
    {"closure/complement",     c_quadratic_ops,      2000,    2000},
    {"binary/edge-list io",    c_io_roundtrip,       M,       100000},
    {"graphml/gml/json/dot",   c_io_export,          M,       10000},
    {"generators",             c_generators,         M,       100000},
    {"mutation + compaction",  c_mutation,           M,       100000},
    {"list<->matrix convert",  c_convert,            0,       100000},
    /* advanced */
    {"adv tree traversal",     c_adv_tree,           M,       100000},
    {"adv refine: cycle/star", c_adv_refine_cycle,   M,       100000},
    {"adv refine: path",       c_adv_refine_path,    3000,    3000},
    {"adv iso (asymmetric)",   c_adv_iso_asym,       3000,    3000},
    {"adv iso (symmetric)",    c_adv_iso_sym,        2000,    2000},
    {"adv subgraph iso",       c_adv_subiso,         M,       100000},
    {"adv widest/flood/jps",   c_adv_paths,          M,       100000},
    {"adv gnn",                c_adv_gnn,            M,       100000},
    {"adv circle layout",      c_adv_circle,         M,       100000},
    {"adv lex-bfs",            c_adv_lexbfs,         10000,   10000},
    {"adv colour coding",      c_adv_colorcoding,    M,       100000},
    {"adv dynamic connectivity", c_adv_dynconn,      M,       100000},
    {"adv small-graph algos",  c_adv_small,          2000,    2000},
};
#undef M

int main(int argc, char **argv)
{
    if (argc > 1) g_max_n = (size_t)strtoull(argv[1], NULL, 10);
    if (argc > 2) g_matrix_max = (size_t)strtoull(argv[2], NULL, 10);
    if (argc > 3) g_filter = argv[3];
    g_only_rep = getenv("SCALE_REP");
    if (getenv("SCALE_MEM_MB") != NULL) g_mem_mb = (size_t)strtoull(getenv("SCALE_MEM_MB"), NULL, 10);
    if (getenv("SCALE_MIN_N") != NULL) g_min_n = (size_t)strtoull(getenv("SCALE_MIN_N"), NULL, 10);
    if (getenv("SCALE_TIMEOUT") != NULL) g_timeout = (unsigned)atoi(getenv("SCALE_TIMEOUT"));
    const size_t sizes[] = {1000, 10000, 100000, 500000, 1000000};
    const GraphRep reps[] = {GRAPH_LIST, GRAPH_MATRIX};
    const double t_start = now();
    int ran = 0, passed = 0;

    printf("%-26s %-7s %9s %10s %9s  %s\n", "case", "rep", "vertices", "time", "peak mem", "result");
    for (int r = 0; r < 2; r++) {
        for (size_t s = 0; s < sizeof sizes / sizeof *sizes; s++) {
            const size_t n = sizes[s];
            if ((g_only_rep != NULL && strcmp(g_only_rep, r == 0 ? "list" : "matrix") != 0) ||
                n < g_min_n || n > g_max_n || (reps[r] == GRAPH_MATRIX && n > g_matrix_max)) continue;
            g_rep = reps[r];
            g_n   = n;
            for (size_t i = 0; i < sizeof CASES / sizeof *CASES; i++) {
                const Case *cs = &CASES[i];
                const size_t cap = reps[r] == GRAPH_LIST ? cs->cap_list : cs->cap_matrix;
                if (n > cap || (g_filter != NULL && strstr(cs->name, g_filter) == NULL)) continue;
                g_case_name   = cs->name;
                g_case_failed = 0;
                fflush(stdout);
                fflush(stderr);
                const double t0 = now();
                const pid_t pid = fork();
                if (pid < 0) { perror("fork"); return 4; }
                if (pid == 0) {                       /* child: build what the case needs, run it */
                    signal(SIGALRM, on_alarm);
                    alarm(g_timeout);
                    if (g_mem_mb != 0) {              /* behave as on a machine with this much memory */
                        struct rlimit lim = {(rlim_t)g_mem_mb << 20, (rlim_t)g_mem_mb << 20};
                        setrlimit(RLIMIT_AS, &lim);
                    }
                    cs->fn();
                    _exit(g_case_failed > 200 ? 200 : g_case_failed);
                }
                int status = 0;
                struct rusage ru;
                wait4(pid, &status, 0, &ru);
                const long peak_mb = ru.ru_maxrss / 1024;   /* Linux reports KiB */
                const double dt = now() - t0;
                char verdict[48] = "ok";
                int  bad = 0;
                if (WIFEXITED(status) && WEXITSTATUS(status) == EXIT_TIMEOUT) {
                    snprintf(verdict, sizeof verdict, "TIMEOUT (> %u s)", g_timeout);
                    bad = 1;
                } else if (WIFSIGNALED(status)) {
                    snprintf(verdict, sizeof verdict, "CRASHED (signal %d)", WTERMSIG(status));
                    bad = 1;
                } else if (WIFEXITED(status) && WEXITSTATUS(status) == 4) {
                    snprintf(verdict, sizeof verdict, "OUT OF MEMORY");
                    bad = 1;
                } else if (WIFEXITED(status) && WEXITSTATUS(status) != 0) {
                    snprintf(verdict, sizeof verdict, "FAILED (%d checks)", WEXITSTATUS(status));
                    bad = WEXITSTATUS(status);
                }
                ran++;
                if (bad == 0) passed++; else g_failures += bad;
                printf("%-26s %-7s %9zu %8.0f ms %6ld MB  %s\n", cs->name, reps[r] == GRAPH_LIST ? "list" : "matrix",
                       n, dt * 1000.0, peak_mb, verdict);
                fflush(stdout);
            }
        }
    }
    printf("\n%d of %d case runs passed, %d failed expectations, %.1f s total\n", passed, ran, g_failures,
           now() - t_start);
    return g_failures == 0 ? 0 : 1;
}
