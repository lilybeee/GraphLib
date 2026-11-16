/*
 * test_adaptive.c — tests for the list backend adapter, migration, the
 * versioned cache and the observer.
 *
 * The matrix side is played by the stub backend (stub_backend.c) until the
 * real adjacency-matrix backend is merged in.
 *
 * Compile & run:  make adaptive_test && ./adaptive_test
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cache.h"
#include "list_backend.h"
#include "migrate.h"
#include "policy.h"
#include "stub_backend.h"

/* ------------------------------------------------------------------ */
/*  Helpers                                                             */
/* ------------------------------------------------------------------ */

static unsigned long long rng_state = 0x9E3779B97F4A7C15ULL;

static unsigned rng(unsigned bound)
{
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return (unsigned)(rng_state % bound);
}

/** Two graphs, possibly in different representations, hold the same thing. */
static void assert_same_graph(const GraphBackend *a, const void *ia,
                              const GraphBackend *b, const void *ib)
{
    assert(a->is_directed(ia) == b->is_directed(ib));
    assert(a->id_bound(ia) == b->id_bound(ib));
    assert(a->vertex_count(ia) == b->vertex_count(ib));
    assert(a->edge_count(ia) == b->edge_count(ib));

    const size_t n = a->id_bound(ia);
    for (size_t u = 0; u < n; ++u) {
        assert(a->is_alive(ia, (int)u) == b->is_alive(ib, (int)u));
        assert(a->degree(ia, (int)u) == b->degree(ib, (int)u));
        for (size_t v = 0; v < n; ++v) {
            double wa = 0.0, wb = 0.0;
            const bool ha = a->has_edge(ia, (int)u, (int)v, &wa);
            const bool hb = b->has_edge(ib, (int)u, (int)v, &wb);
            assert(ha == hb);
            if (ha) {
                assert(wa == wb);
            }
        }
    }
}

/** Random mutations through the vtable; returns the number applied. */
static void mutate_randomly(const GraphBackend *b, void *g, int steps, int max_ids)
{
    for (int i = 0; i < steps; ++i) {
        const unsigned op = rng(100);
        const int n = (int)b->id_bound(g);
        const int u = n == 0 ? 0 : (int)rng((unsigned)n);
        const int v = n == 0 ? 0 : (int)rng((unsigned)n);
        if (op < 10) {
            if (n < max_ids || b->vertex_count(g) < (size_t)n) {
                (void)b->add_vertex(g);
            }
        } else if (op < 18) {
            (void)b->remove_vertex(g, u);
        } else if (op < 75) {
            (void)b->add_edge(g, u, v, (double)rng(40) / 4.0);
        } else {
            (void)b->remove_edge(g, u, v);
        }
    }
}

/* ------------------------------------------------------------------ */
/*  1. The list backend through its vtable                              */
/* ------------------------------------------------------------------ */

static void test_list_backend_vtable(void)
{
    const GraphBackend *b = list_backend();
    assert(strcmp(b->name, "list") == 0);

    /* Everything the observer and migration rely on is present. */
    assert(b->create && b->destroy && b->iter_begin && b->iter_next);
    assert(b->add_new_edge && b->estimate_cost && b->estimate_memory &&
           b->estimate_build_cost);
    assert(b->bfs && b->dfs && b->dijkstra && b->mst_kruskal);

    void *g = b->create(4, false);
    assert(g != NULL);
    assert(b->add_edge(g, 0, 1, 1.5));
    assert(b->add_edge(g, 0, 2, 2.5));
    assert(b->add_new_edge(g, 2, 3, 0.5));
    assert(b->edge_count(g) == 3);
    assert(b->degree(g, 0) == 2);

    /* Iterate vertex 0's neighbourhood: {1: 1.5, 2: 2.5}. */
    GraphIter it;
    int v;
    double w, total = 0.0;
    int count = 0, seen_sum = 0;
    b->iter_begin(g, 0, &it);
    while (b->iter_next(&it, &v, &w)) {
        ++count;
        seen_sum += v;
        total += w;
    }
    assert(count == 2 && seen_sum == 3 && fabs(total - 4.0) < 1e-12);

    /* Dead and out-of-range vertices iterate as empty. */
    assert(b->remove_vertex(g, 1));
    b->iter_begin(g, 1, &it);
    assert(!b->iter_next(&it, &v, &w));
    b->iter_begin(g, 99, &it);
    assert(!b->iter_next(&it, &v, &w));

    /* Kernels run through the vtable. */
    int order[4], parent[4];
    assert(b->bfs(g, 0, order, parent) == 3); /* 0, 2, 3 */

    /* The cost model grows with the work asked of it. */
    GraphShape small = { 100, 100, 200, false, false};
    GraphShape big   = { 100, 100, 2000, false, false };
    Workload algo;
    workload_reset(&algo);
    algo.n[OP_ALGO_RUN] = 10;
    assert(b->estimate_cost(&big, &algo) > b->estimate_cost(&small, &algo));
    assert(b->estimate_memory(&big) > b->estimate_memory(&small));
    assert(b->estimate_build_cost(&big) > b->estimate_build_cost(&small));

    b->destroy(g);
}

/* ------------------------------------------------------------------ */
/*  2. List and stub stay in lock-step under random operations          */
/* ------------------------------------------------------------------ */

static void lockstep(bool directed)
{
    const GraphBackend *l = list_backend();
    const GraphBackend *s = stub_backend();
    void *gl = l->create(6, directed);
    void *gs = s->create(6, directed);
    assert(gl != NULL && gs != NULL);

    for (int step = 0; step < 1500; ++step) {
        const unsigned op = rng(100);
        const int n = (int)l->id_bound(gl);
        const int u = (int)rng((unsigned)n + 1);   /* sometimes out of range */
        const int v = (int)rng((unsigned)n + 1);

        if (op < 10) {
            if (n < 40 || l->vertex_count(gl) < (size_t)n) {
                assert(l->add_vertex(gl) == s->add_vertex(gs));
            }
        } else if (op < 18) {
            assert(l->remove_vertex(gl, u) == s->remove_vertex(gs, u));
        } else if (op < 75) {
            const double w = (double)rng(40) / 4.0;
            assert(l->add_edge(gl, u, v, w) == s->add_edge(gs, u, v, w));
        } else {
            assert(l->remove_edge(gl, u, v) == s->remove_edge(gs, u, v));
        }
        assert(l->edge_count(gl) == s->edge_count(gs));
        assert(l->vertex_count(gl) == s->vertex_count(gs));
        if (step % 50 == 0) {
            assert_same_graph(l, gl, s, gs);
        }
    }
    assert_same_graph(l, gl, s, gs);
    l->destroy(gl);
    s->destroy(gs);
}

static void test_lockstep(void)
{
    stub_set_cost_mode(STUB_COST_MATRIX, 1.0);
    lockstep(false);
    lockstep(true);
}

/* ------------------------------------------------------------------ */
/*  3. Migration preserves IDs and content                              */
/* ------------------------------------------------------------------ */

static void clone_round_trip(bool directed)
{
    const GraphBackend *l = list_backend();
    const GraphBackend *s = stub_backend();

    void *src = l->create(12, directed);
    assert(src != NULL);
    mutate_randomly(l, src, 400, 30);
    assert(l->vertex_count(src) < l->id_bound(src) || l->id_bound(src) > 0);

    const uint64_t version = l->version(src);
    void *mid = graph_clone(l, src, s);       /* list -> matrix */
    assert(mid != NULL);
    assert_same_graph(l, src, s, mid);
    assert(l->version(src) == version);       /* the source was not touched */

    void *back = graph_clone(s, mid, l);      /* matrix -> list */
    assert(back != NULL);
    assert_same_graph(s, mid, l, back);
    assert_same_graph(l, src, l, back);

    /* Dead IDs come back dead and reusable: the next vertex reuses one if any. */
    const size_t dead = l->id_bound(back) - l->vertex_count(back);
    const int fresh = l->add_vertex(back);
    assert(fresh >= 0);
    if (dead != 0) {
        assert((size_t)fresh < l->id_bound(back) - 1 + 1 &&
               l->id_bound(back) == l->id_bound(src));
    }

    l->destroy(src);
    s->destroy(mid);
    l->destroy(back);
}

static void test_clone(void)
{
    stub_set_cost_mode(STUB_COST_MATRIX, 1.0);
    clone_round_trip(false);
    clone_round_trip(true);

    /* Cloning an empty graph works, and failure leaves the caller with NULL. */
    const GraphBackend *l = list_backend();
    void *empty = l->create(0, true);
    void *copy  = graph_clone(l, empty, stub_backend());
    assert(copy != NULL && stub_backend()->id_bound(copy) == 0);
    stub_backend()->destroy(copy);
    l->destroy(empty);

    void *one = l->create(3, false);
    stub_set_fail_create(true);
    assert(graph_clone(l, one, stub_backend()) == NULL);
    stub_set_fail_create(false);
    assert(graph_clone(NULL, one, stub_backend()) == NULL);
    l->destroy(one);
}

/* ------------------------------------------------------------------ */
/*  4. The cache                                                        */
/* ------------------------------------------------------------------ */

static CacheKey key(uint32_t algo, int a, int b)
{
    CacheKey k = { algo, a, b };
    return k;
}

static void test_cache_basics(void)
{
    Cache *c = cache_create(8, 1024);
    assert(c != NULL);

    int result[4] = { 1, 2, 3, 4 };
    size_t size = 0;

    assert(cache_lookup(c, key(CACHE_BFS, 0, 0), 1, &size) == NULL);   /* miss */
    assert(cache_store(c, key(CACHE_BFS, 0, 0), 1, result, sizeof(result)));
    const int *hit = cache_lookup(c, key(CACHE_BFS, 0, 0), 1, &size);
    assert(hit != NULL && size == sizeof(result) && memcmp(hit, result, size) == 0);

    /* The cache stores a copy. */
    result[0] = 99;
    hit = cache_lookup(c, key(CACHE_BFS, 0, 0), 1, NULL);
    assert(hit != NULL && hit[0] == 1);

    /* Different algo / arguments are different keys. */
    assert(cache_lookup(c, key(CACHE_DFS, 0, 0), 1, NULL) == NULL);
    assert(cache_lookup(c, key(CACHE_BFS, 1, 0), 1, NULL) == NULL);
    assert(cache_lookup(c, key(CACHE_BFS, 0, 1), 1, NULL) == NULL);

    /* Storing the same key again replaces the value. */
    int updated[2] = { 7, 8 };
    assert(cache_store(c, key(CACHE_BFS, 0, 0), 1, updated, sizeof(updated)));
    hit = cache_lookup(c, key(CACHE_BFS, 0, 0), 1, &size);
    assert(hit != NULL && size == sizeof(updated) && hit[1] == 8);
    assert(cache_stats(c).entries == 1);

    /* A zero-byte result (e.g. "no cycle") is a valid entry. */
    assert(cache_store(c, key(CACHE_CYCLE, 0, 0), 1, NULL, 0));
    assert(cache_lookup(c, key(CACHE_CYCLE, 0, 0), 1, &size) != NULL && size == 0);

    /* A new graph version invalidates everything on the next access. */
    assert(cache_lookup(c, key(CACHE_BFS, 0, 0), 2, NULL) == NULL);
    CacheStats st = cache_stats(c);
    assert(st.entries == 0 && st.bytes == 0 && st.invalidations == 2);
    assert(st.hits >= 3 && st.misses >= 5);

    /* ...and an old entry cannot come back under the old version. */
    assert(cache_store(c, key(CACHE_BFS, 0, 0), 2, result, sizeof(result)));
    assert(cache_lookup(c, key(CACHE_BFS, 0, 0), 1, NULL) == NULL);

    /* NULL-safety. */
    assert(cache_lookup(NULL, key(CACHE_BFS, 0, 0), 1, NULL) == NULL);
    assert(!cache_store(NULL, key(CACHE_BFS, 0, 0), 1, result, 4));
    cache_destroy(NULL);
    cache_destroy(c);
}

static void test_cache_eviction(void)
{
    /* Entry limit: the least recently used entry goes first. */
    Cache *c = cache_create(3, 1 << 20);
    int x = 0;
    for (int i = 0; i < 3; ++i) {
        assert(cache_store(c, key(CACHE_BFS, i, 0), 1, &x, sizeof(x)));
    }
    assert(cache_lookup(c, key(CACHE_BFS, 0, 0), 1, NULL) != NULL); /* touch 0 */
    assert(cache_store(c, key(CACHE_BFS, 3, 0), 1, &x, sizeof(x)));  /* evicts 1 */
    assert(cache_lookup(c, key(CACHE_BFS, 1, 0), 1, NULL) == NULL);
    assert(cache_lookup(c, key(CACHE_BFS, 0, 0), 1, NULL) != NULL);
    assert(cache_lookup(c, key(CACHE_BFS, 2, 0), 1, NULL) != NULL);
    assert(cache_lookup(c, key(CACHE_BFS, 3, 0), 1, NULL) != NULL);
    assert(cache_stats(c).evictions == 1 && cache_stats(c).entries == 3);
    cache_destroy(c);

    /* Byte limit: big results push out old ones; an oversized one is refused. */
    c = cache_create(100, 100);
    char blob[60] = { 0 };
    assert(cache_store(c, key(CACHE_DIJKSTRA, 1, 0), 1, blob, 60));
    assert(cache_store(c, key(CACHE_DIJKSTRA, 2, 0), 1, blob, 60)); /* evicts 1 */
    assert(cache_lookup(c, key(CACHE_DIJKSTRA, 1, 0), 1, NULL) == NULL);
    assert(cache_stats(c).bytes == 60);
    char huge[101] = { 0 };
    assert(!cache_store(c, key(CACHE_DIJKSTRA, 3, 0), 1, huge, sizeof(huge)));
    assert(cache_lookup(c, key(CACHE_DIJKSTRA, 2, 0), 1, NULL) != NULL); /* untouched */
    cache_destroy(c);

    /* A disabled cache stores nothing. */
    c = cache_create(0, 100);
    assert(!cache_store(c, key(CACHE_BFS, 0, 0), 1, &x, sizeof(x)));
    cache_destroy(c);
}

/* Many keys force the hash table to grow; every value must survive. */
static void test_cache_many_keys(void)
{
    const int n = 5000;
    Cache *c = cache_create((size_t)n, (size_t)n * sizeof(int));
    for (int i = 0; i < n; ++i) {
        assert(cache_store(c, key(CACHE_USER + (uint32_t)(i % 3), i, -i), 7,
                           &i, sizeof(i)));
    }
    for (int i = 0; i < n; ++i) {
        size_t size;
        const int *v = cache_lookup(c, key(CACHE_USER + (uint32_t)(i % 3), i, -i), 7, &size);
        assert(v != NULL && size == sizeof(int) && *v == i);
    }
    assert(cache_stats(c).entries == (size_t)n && cache_stats(c).evictions == 0);
    cache_clear(c);
    assert(cache_stats(c).entries == 0 && cache_stats(c).bytes == 0);
    cache_destroy(c);
}

/* A cache fed from a changing graph never serves a stale answer. */
static void test_cache_tracks_graph_version(void)
{
    const GraphBackend *b = list_backend();
    void *g = b->create(5, false);
    Cache *c = cache_create(16, 4096);
    const CacheKey k = key(CACHE_COMPONENTS, 0, 0);

    for (int round = 0; round < 40; ++round) {
        (void)b->add_edge(g, (int)rng(5), (int)rng(5), 1.0);
        int comp[5];
        const size_t fresh = b->connected_components(g, comp);

        size_t size;
        const int *cached = cache_lookup(c, k, b->version(g), &size);
        if (cached != NULL) {
            assert(size == sizeof(comp) + sizeof(fresh));
            size_t stored;
            memcpy(&stored, cached + 5, sizeof(stored));
            assert(stored == fresh && memcmp(cached, comp, sizeof(comp)) == 0);
        } else {
            unsigned char buf[sizeof(comp) + sizeof(fresh)];
            memcpy(buf, comp, sizeof(comp));
            memcpy(buf + sizeof(comp), &fresh, sizeof(fresh));
            assert(cache_store(c, k, b->version(g), buf, sizeof(buf)));
        }
    }
    assert(cache_stats(c).invalidations > 0);
    cache_destroy(c);
    b->destroy(g);
}

/* ------------------------------------------------------------------ */
/*  5. The observer                                                     */
/* ------------------------------------------------------------------ */

typedef struct {
    const GraphBackend *backend;
    void               *impl;
} Held;

static const GraphBackend *candidates[2];

static void setup_candidates(void)
{
    candidates[0] = list_backend();
    candidates[1] = stub_backend();
}

/** Undirected list graph with `n` vertices and each pair joined with probability p%. */
static Held make_graph(size_t n, unsigned percent)
{
    Held h = { list_backend(), NULL };
    h.impl = h.backend->create(n, false);
    assert(h.impl != NULL);
    for (size_t u = 0; u < n; ++u) {
        for (size_t v = u + 1; v < n; ++v) {
            if (rng(100) < percent) {
                assert(h.backend->add_edge(h.impl, (int)u, (int)v, 1.0 + rng(9)));
            }
        }
    }
    return h;
}

static void feed(Observer *obs, OpKind op, int times)
{
    for (int i = 0; i < times; ++i) {
        observer_note(obs, op);
    }
}

/** feed + one evaluation, the way a facade would after each batch. */
static bool step(Observer *obs, Held *h, OpKind op, int times, Decision *d)
{
    feed(obs, op, times);
    return observer_maybe_migrate(obs, &h->backend, &h->impl, candidates, 2, d);
}

static void test_observer_clock(void)
{
    Observer *obs = observer_create(NULL);
    assert(obs != NULL);
    assert(!observer_due(obs));
    feed(obs, OP_CACHE_HIT, 5000);            /* hits do no graph work */
    feed(obs, OP_CACHE_MISS, 5000);
    assert(!observer_due(obs));
    feed(obs, OP_HAS_EDGE, 999);
    assert(!observer_due(obs));
    feed(obs, OP_HAS_EDGE, 1);
    assert(observer_due(obs));
    assert(observer_workload(obs)->n[OP_CACHE_HIT] == 5000);
    assert(workload_hit_rate(observer_workload(obs)) == 0.5);

    /* Before the clock fires, maybe_migrate does nothing at all. */
    Observer *quiet = observer_create(NULL);
    Held h = make_graph(10, 50);
    Decision d;
    assert(!step(quiet, &h, OP_HAS_EDGE, 10, &d));
    assert(d.kind == DECISION_STAY && d.target == NULL);
    h.backend->destroy(h.impl);
    observer_destroy(quiet);
    observer_destroy(obs);
    observer_destroy(NULL);
}

/** A dense graph under point-lookup load moves to the matrix, and stays. */
static void test_observer_migrates_and_settles(void)
{
    stub_set_cost_mode(STUB_COST_MATRIX, 1.0);
    setup_candidates();

    Held h = make_graph(200, 50);                       /* avg degree ~100 */
    void *reference = graph_clone(h.backend, h.impl, h.backend);
    const GraphBackend *ref_backend = h.backend;

    FILE *log = tmpfile();
    assert(log != NULL);
    ObserverConfig cfg = observer_default_config();
    cfg.log = log;
    Observer *obs = observer_create(&cfg);

    Decision d;
    assert(step(obs, &h, OP_HAS_EDGE, 1000, &d));       /* first check migrates */
    assert(d.kind == DECISION_MIGRATE && d.target == stub_backend());
    assert(d.cost_current > cfg.hysteresis * d.cost_target);
    assert(h.backend == stub_backend());
    assert(observer_migrations(obs) == 1);
    assert_same_graph(ref_backend, reference, h.backend, h.impl);  /* content kept */

    /* The decision was logged. */
    rewind(log);
    char line[256] = { 0 };
    assert(fgets(line, sizeof(line), log) != NULL);
    assert(strstr(line, "list -> matrix-stub") != NULL);

    /* The same workload for many more windows: no ping-pong. */
    for (int i = 0; i < 25; ++i) {
        assert(!step(obs, &h, OP_HAS_EDGE, 1000, &d));
        assert(h.backend == stub_backend());
    }
    assert(observer_migrations(obs) == 1);

    /* The graph thins out and the workload turns to whole-graph algorithms:
     * the list is now the better home, and the observer notices. */
    const int n = (int)h.backend->id_bound(h.impl);
    for (int u = 0; u < n; ++u) {
        for (int v = u + 2; v < n; ++v) {
            (void)h.backend->remove_edge(h.impl, u, v);
        }
        if (u + 1 < n && !h.backend->has_edge(h.impl, u, u + 1, NULL)) {
            assert(h.backend->add_edge(h.impl, u, u + 1, 1.0));
        }
    }
    assert(h.backend->edge_count(h.impl) == (size_t)n - 1);   /* exactly a path */
    bool back = false;
    for (int i = 0; i < 5 && !back; ++i) {
        back = step(obs, &h, OP_ALGO_RUN, 1000, &d);
    }
    assert(back && h.backend == list_backend());
    assert(observer_migrations(obs) == 2);
    assert(h.backend->edge_count(h.impl) == (size_t)n - 1);
    assert(h.backend->has_edge(h.impl, 5, 6, NULL));
    assert(!h.backend->has_edge(h.impl, 5, 7, NULL));

    fclose(log);
    observer_destroy(obs);
    h.backend->destroy(h.impl);
    ref_backend->destroy(reference);
}

/** Hysteresis and payback each block a migration at the right boundary. */
static void test_observer_thresholds(void)
{
    setup_candidates();
    ObserverConfig cfg = observer_default_config();

    /* Stub only 1.25x cheaper: better, but inside the 1.5x hysteresis band. */
    stub_set_cost_mode(STUB_COST_SCALED, 0.8);
    Held h = make_graph(200, 50);
    Observer *obs = observer_create(&cfg);
    Decision d;
    assert(!step(obs, &h, OP_HAS_EDGE, 5000, &d));
    assert(d.kind == DECISION_NOT_ENOUGH_GAIN && h.backend == list_backend());
    observer_destroy(obs);

    /* 2x cheaper, but too little work in the window to repay the rebuild. */
    stub_set_cost_mode(STUB_COST_SCALED, 0.5);
    obs = observer_create(&cfg);
    assert(!step(obs, &h, OP_HAS_EDGE, 1000, &d));
    assert(d.kind == DECISION_NOT_WORTH_IT && h.backend == list_backend());
    observer_destroy(obs);

    /* 2x cheaper with plenty of work: now it pays. */
    obs = observer_create(&cfg);
    assert(step(obs, &h, OP_HAS_EDGE, 5000, &d));
    assert(d.kind == DECISION_MIGRATE && h.backend == stub_backend());
    observer_destroy(obs);
    h.backend->destroy(h.impl);

    /* Stub more expensive than the list: never moves. */
    stub_set_cost_mode(STUB_COST_SCALED, 2.0);
    h = make_graph(200, 50);
    obs = observer_create(&cfg);
    assert(!step(obs, &h, OP_HAS_EDGE, 5000, &d));
    assert(d.kind == DECISION_STAY);
    observer_destroy(obs);
    h.backend->destroy(h.impl);
}

/** After a migration the observer holds still until the cooldown has passed. */
static void test_observer_cooldown(void)
{
    setup_candidates();
    ObserverConfig cfg = observer_default_config();
    cfg.cooldown_ops = 20000;

    stub_set_cost_mode(STUB_COST_SCALED, 0.5);
    Held h = make_graph(200, 50);
    Observer *obs = observer_create(&cfg);
    Decision d;
    assert(step(obs, &h, OP_HAS_EDGE, 5000, &d));        /* list -> stub */

    /* The market flips: the stub is now 2x dearer than the list. */
    stub_set_cost_mode(STUB_COST_SCALED, 2.0);
    assert(!step(obs, &h, OP_HAS_EDGE, 5000, &d));
    assert(d.kind == DECISION_COOLDOWN && h.backend == stub_backend());

    /* Once enough operations have passed, it is free to move back. */
    assert(step(obs, &h, OP_HAS_EDGE, 20000, &d));
    assert(d.kind == DECISION_MIGRATE && h.backend == list_backend());
    assert(observer_migrations(obs) == 2);

    observer_destroy(obs);
    h.backend->destroy(h.impl);
    stub_set_cost_mode(STUB_COST_MATRIX, 1.0);
}

/** Memory is a hard constraint, in both directions. */
static void test_observer_memory_budget(void)
{
    stub_set_cost_mode(STUB_COST_MATRIX, 1.0);
    setup_candidates();

    /* 1) The matrix would win on speed but does not fit: stay on the list. */
    Held h = make_graph(1000, 20);                       /* avg degree ~200 */
    GraphShape shape = graph_shape_of(h.backend, h.impl);
    const size_t list_bytes   = h.backend->memory_bytes(h.impl);
    const size_t matrix_bytes = stub_backend()->estimate_memory(&shape);
    assert(list_bytes < matrix_bytes);

    ObserverConfig cfg = observer_default_config();
    Observer *control = observer_create(&cfg);           /* no budget: it moves */
    Held ctl = make_graph(1000, 20);
    Decision d;
    assert(step(control, &ctl, OP_HAS_EDGE, 10000, &d));
    assert(ctl.backend == stub_backend());
    observer_destroy(control);
    ctl.backend->destroy(ctl.impl);

    cfg.memory_budget = (list_bytes + matrix_bytes) / 2;  /* list fits, matrix not */
    Observer *obs = observer_create(&cfg);
    assert(!step(obs, &h, OP_HAS_EDGE, 10000, &d));
    assert(d.kind == DECISION_STAY && h.backend == list_backend());
    observer_destroy(obs);
    h.backend->destroy(h.impl);

    /* 2) A matrix that has outgrown the budget is forced back to the list,
     *    even with no useful workload and no hysteresis margin. */
    Held big = { stub_backend(), stub_backend()->create(300, false) };
    assert(big.impl != NULL);
    for (int i = 0; i + 1 < 300; ++i) {
        assert(big.backend->add_edge(big.impl, i, i + 1, 1.0));
    }
    cfg = observer_default_config();
    cfg.memory_budget = 100 * 1024;
    assert(big.backend->memory_bytes(big.impl) > cfg.memory_budget);
    obs = observer_create(&cfg);
    assert(step(obs, &big, OP_ADD_VERTEX, 1000, &d));
    assert(d.kind == DECISION_FORCED_MEMORY && big.backend == list_backend());
    assert(big.backend->edge_count(big.impl) == 299);
    assert(big.backend->has_edge(big.impl, 41, 42, NULL));
    observer_destroy(obs);
    big.backend->destroy(big.impl);

    /* 3) Over budget with nowhere to go: stay put rather than fail. */
    h = make_graph(100, 50);
    cfg = observer_default_config();
    cfg.memory_budget = 1;
    obs = observer_create(&cfg);
    assert(!step(obs, &h, OP_HAS_EDGE, 1000, &d));
    assert(d.kind == DECISION_STAY && h.backend == list_backend());
    observer_destroy(obs);
    h.backend->destroy(h.impl);
}

/** A rebuild that fails must leave the live graph exactly as it was. */
static void test_observer_failed_migration(void)
{
    stub_set_cost_mode(STUB_COST_MATRIX, 1.0);
    setup_candidates();

    Held h = make_graph(200, 50);
    void *reference = graph_clone(h.backend, h.impl, h.backend);
    Observer *obs = observer_create(NULL);
    Decision d;

    stub_set_fail_create(true);
    assert(!step(obs, &h, OP_HAS_EDGE, 1000, &d));
    stub_set_fail_create(false);
    assert(d.kind == DECISION_MIGRATE);                  /* it tried... */
    assert(h.backend == list_backend());                 /* ...and kept the list */
    assert(observer_migrations(obs) == 0);
    assert_same_graph(list_backend(), reference, h.backend, h.impl);

    /* The next check succeeds. */
    assert(step(obs, &h, OP_HAS_EDGE, 1000, &d));
    assert(h.backend == stub_backend());

    observer_destroy(obs);
    h.backend->destroy(h.impl);
    list_backend()->destroy(reference);
}

/* ------------------------------------------------------------------ */
/*  main                                                                */
/* ------------------------------------------------------------------ */

int main(void)
{
    test_list_backend_vtable();
    test_lockstep();
    test_clone();
    test_cache_basics();
    test_cache_eviction();
    test_cache_many_keys();
    test_cache_tracks_graph_version();
    test_observer_clock();
    test_observer_migrates_and_settles();
    test_observer_thresholds();
    test_observer_cooldown();
    test_observer_memory_budget();
    test_observer_failed_migration();
    puts("All adaptive-layer tests passed.");
    return 0;
}
