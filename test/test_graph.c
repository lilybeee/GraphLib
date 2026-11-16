/*
 * test_graph.c — tests for the unified Graph API (graph.h) and NamedGraph.
 *
 * Run: make build/graph_test && ./build/graph_test
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "graph.h"
#include "mapping.h"

/* ---- helpers ---- */

static GraphConfig config(bool directed, GraphRep start, bool adaptive)
{
    GraphConfig c = graph_default_config(directed);
    c.start       = start;
    c.adaptive    = adaptive;
    return c;
}

static Graph *make(bool directed, GraphRep start, bool adaptive, size_t n)
{
    const GraphConfig c = config(directed, start, adaptive);
    Graph *g = graph_create_with(&c, n);
    assert(g != NULL);
    return g;
}

/** A small weighted undirected graph: a 6-cycle with a chord, plus an island. */
static void fill_sample(Graph *g)
{
    for (int i = 0; i < 6; ++i) {
        assert(graph_add_edge(g, i, (i + 1) % 6, 1.0 + i));
    }
    assert(graph_add_edge(g, 0, 3, 10.0));
}

/** Everything observable about a graph, via the public API only. */
typedef struct {
    size_t vertices, bound, edges;
    size_t bfs_count, comps;
    bool   cyclic;
    double dist[128];
    double mst;
} Snapshot;

static Snapshot snapshot(Graph *g)
{
    Snapshot s;
    memset(&s, 0, sizeof(s));
    s.vertices = graph_vertex_count(g);
    s.bound    = graph_id_bound(g);
    s.edges    = graph_edge_count(g);
    int order[128], parent[128], comp[128];
    s.bfs_count = graph_bfs(g, 0, order, parent);
    s.comps     = graph_connected_components(g, comp);
    s.cyclic    = graph_has_cycle(g);
    assert(graph_dijkstra(g, 0, s.dist, parent));
    GraphEdge edges[128];
    size_t    ne = 0;
    s.mst = graph_mst_kruskal(g, edges, &ne);
    return s;
}

static void assert_same(const Snapshot *a, const Snapshot *b)
{
    assert(a->vertices == b->vertices && a->bound == b->bound && a->edges == b->edges);
    assert(a->bfs_count == b->bfs_count && a->comps == b->comps && a->cyclic == b->cyclic);
    for (size_t i = 0; i < a->bound; ++i) {
        assert(a->dist[i] == b->dist[i] || (isinf(a->dist[i]) && isinf(b->dist[i])));
    }
    assert(a->mst == b->mst || (isinf(a->mst) && isinf(b->mst)));
}

/* ---- tests ---- */

static void test_start_on_either_representation(void)
{
    for (int rep = GRAPH_LIST; rep <= GRAPH_MATRIX; ++rep) {
        Graph *g = make(false, (GraphRep)rep, false, 7);
        assert(graph_representation(g) == (GraphRep)rep);
        fill_sample(g);
        assert(graph_vertex_count(g) == 7);
        assert(graph_edge_count(g) == 7);
        double w = 0;
        assert(graph_has_edge(g, 3, 0, &w) && w == 10.0);   /* undirected */
        assert(!graph_has_edge(g, 1, 4, NULL));
        assert(!graph_add_edge(g, 0, 1, 1.0));              /* duplicate  */
        assert(graph_degree(g, 0) == 3);
        graph_destroy(g);
    }
}

static void test_same_answers_on_both(void)
{
    Graph *l = make(false, GRAPH_LIST, false, 7);
    Graph *m = make(false, GRAPH_MATRIX, false, 7);
    fill_sample(l);
    fill_sample(m);
    const Snapshot a = snapshot(l);
    const Snapshot b = snapshot(m);
    assert_same(&a, &b);
    assert(a.comps == 2);            /* the cycle, and the island vertex 6 */
    assert(a.cyclic);
    assert(a.bfs_count == 6);
    assert(isinf(a.mst));            /* disconnected: no spanning tree      */
    graph_destroy(l);
    graph_destroy(m);
}

static void test_manual_convert_preserves_everything(void)
{
    Graph *g = make(false, GRAPH_LIST, false, 7);
    fill_sample(g);
    const Snapshot before = snapshot(g);

    assert(graph_convert(g, GRAPH_MATRIX));
    assert(graph_representation(g) == GRAPH_MATRIX);
    const Snapshot mid = snapshot(g);
    assert_same(&before, &mid);

    assert(graph_convert(g, GRAPH_MATRIX));          /* already there: no-op */
    assert(graph_convert(g, GRAPH_LIST));
    const Snapshot after = snapshot(g);
    assert_same(&before, &after);
    assert(graph_migration_count(g) == 0);           /* manual != automatic  */
    graph_destroy(g);
}

static void test_remove_vertex_and_reuse_and_compact(void)
{
    for (int rep = GRAPH_LIST; rep <= GRAPH_MATRIX; ++rep) {
        Graph *g = make(false, (GraphRep)rep, false, 6);
        fill_sample(g);
        assert(graph_remove_vertex(g, 3));
        assert(!graph_remove_vertex(g, 3));              /* already dead      */
        assert(!graph_is_alive(g, 3));
        assert(graph_vertex_count(g) == 5);
        assert(graph_id_bound(g) == 6);
        assert(graph_edge_count(g) == 4);                /* 2-3, 3-4, 0-3 gone */
        assert(!graph_has_edge(g, 2, 3, NULL));
        assert(!graph_add_edge(g, 3, 4, 1.0));           /* dead endpoint     */
        assert(graph_degree(g, 3) == 0);

        int order[6], dist_parent[6];
        assert(graph_bfs(g, 0, order, dist_parent) == 5);
        assert(dist_parent[3] == -1);
        assert(graph_bfs(g, 3, order, NULL) == 0);       /* dead start        */

        /* The dead ID is handed out again. */
        assert(graph_add_vertex(g) == 3);
        assert(graph_is_alive(g, 3));
        assert(graph_vertex_count(g) == 6);

        /* Compaction renumbers, and reports how. */
        assert(graph_remove_vertex(g, 1));
        assert(graph_remove_vertex(g, 2));
        assert(graph_needs_compaction(g));
        int remap[6];
        assert(graph_compact(g, remap));
        assert(graph_id_bound(g) == 4 && graph_vertex_count(g) == 4);
        assert(remap[0] == 0 && remap[1] == -1 && remap[2] == -1 && remap[3] == 1);
        graph_destroy(g);
    }
}

static void test_neighbors(void)
{
    for (int rep = GRAPH_LIST; rep <= GRAPH_MATRIX; ++rep) {
        Graph *g = make(true, (GraphRep)rep, false, 4);
        assert(graph_add_edge(g, 0, 1, 1.5));
        assert(graph_add_edge(g, 0, 3, 2.5));
        int    nb[4];
        double w[4];
        assert(graph_neighbors(g, 0, nb, w, 4) == 2);
        double sum = 0;
        int    mask = 0;
        for (int i = 0; i < 2; ++i) {
            sum += w[i];
            mask |= 1 << nb[i];
        }
        assert(sum == 4.0 && mask == ((1 << 1) | (1 << 3)));
        assert(graph_neighbors(g, 0, nb, NULL, 1) == 2);  /* reports the true degree */
        assert(graph_neighbors(g, 2, nb, w, 4) == 0);
        graph_destroy(g);
    }
}

static void test_directed_algorithms(void)
{
    Graph *g = make(true, GRAPH_MATRIX, false, 5);
    assert(graph_add_edge(g, 0, 1, 1));
    assert(graph_add_edge(g, 1, 2, 1));
    assert(graph_add_edge(g, 0, 3, 5));
    assert(graph_add_edge(g, 3, 2, 1));
    int order[5];
    assert(!graph_has_cycle(g));
    assert(graph_topological_sort(g, order));
    assert(order[0] == 0);
    double dist[5];
    assert(graph_dijkstra(g, 0, dist, NULL));
    assert(dist[2] == 2.0 && isinf(dist[4]));
    assert(graph_add_edge(g, 2, 0, 1));
    assert(graph_has_cycle(g));                       /* cache must not go stale */
    assert(!graph_topological_sort(g, order));
    graph_destroy(g);
}

static void test_cache_never_goes_stale(void)
{
    Graph *g = make(false, GRAPH_LIST, false, 4);
    assert(graph_add_edge(g, 0, 1, 1));
    int order[4], parent[4];
    assert(graph_bfs(g, 0, order, parent) == 2);
    assert(graph_bfs(g, 0, order, parent) == 2);      /* served from cache  */
    assert(graph_add_edge(g, 1, 2, 1));
    assert(graph_bfs(g, 0, order, parent) == 3);      /* version changed    */
    assert(parent[2] == 1);
    assert(graph_remove_vertex(g, 1));
    assert(graph_bfs(g, 0, order, parent) == 1);
    graph_destroy(g);
}

static void test_cache_stats(void)
{
    CacheStats st;
    assert(!graph_cache_stats(NULL, &st));
    Graph *g = make(false, GRAPH_LIST, false, 4);
    assert(!graph_cache_stats(g, NULL));
    assert(graph_cache_stats(g, &st));
    const uint64_t h0 = st.hits, m0 = st.misses, i0 = st.invalidations;
    assert(graph_add_edge(g, 0, 1, 1));
    int order[4], parent[4];
    assert(graph_bfs(g, 0, order, parent) == 2);      /* miss, then stored  */
    assert(graph_bfs(g, 0, order, parent) == 2);      /* hit                */
    assert(graph_cache_stats(g, &st));
    assert(st.misses == m0 + 1 && st.hits == h0 + 1 && st.entries == 1);
    assert(graph_add_edge(g, 1, 2, 1));
    assert(graph_bfs(g, 0, order, parent) == 3);      /* stale -> miss      */
    assert(graph_cache_stats(g, &st));
    assert(st.misses == m0 + 2 && st.invalidations > i0);
    graph_destroy(g);
    GraphConfig c = config(false, GRAPH_LIST, false);
    c.cache_entries = 0;
    g = graph_create_with(&c, 3);
    assert(graph_add_edge(g, 0, 1, 1));
    assert(graph_bfs(g, 0, order, parent) == 2);
    assert(graph_bfs(g, 0, order, parent) == 2);
    if (graph_cache_stats(g, &st)) {
        assert(st.hits == 0 && st.entries == 0);      /* disabled: never serves */
    }
    graph_destroy(g);
}

static void test_adaptive_off_never_migrates(void)
{
    GraphConfig c = config(false, GRAPH_LIST, false);
    c.observer.check_interval = 50;
    c.observer.cooldown_ops   = 0;
    Graph *g = graph_create_with(&c, 60);
    for (int u = 0; u < 60; ++u) {
        for (int v = u + 1; v < 60; ++v) {
            assert(graph_add_edge(g, u, v, 1.0));     /* complete graph     */
        }
    }
    for (int i = 0; i < 20000; ++i) {
        (void)graph_has_edge(g, i % 60, (i * 7 + 1) % 60, NULL);
    }
    assert(graph_representation(g) == GRAPH_LIST);
    assert(graph_migration_count(g) == 0);
    graph_destroy(g);
}

static void test_adaptive_on_migrates_to_matrix_for_point_queries(void)
{
    GraphConfig c = config(false, GRAPH_LIST, true);
    c.observer.check_interval = 200;
    c.observer.cooldown_ops   = 0;
    Graph *g = graph_create_with(&c, 120);
    for (int u = 0; u < 120; ++u) {
        for (int v = u + 1; v < 120; ++v) {
            assert(graph_add_edge(g, u, v, 1.0));     /* dense              */
        }
    }
    const Snapshot before = snapshot(g);
    for (int i = 0; i < 20000 && graph_representation(g) == GRAPH_LIST; ++i) {
        (void)graph_has_edge(g, i % 120, (i * 7 + 1) % 120, NULL);
    }
    assert(graph_representation(g) == GRAPH_MATRIX);
    assert(graph_migration_count(g) >= 1);
    const Snapshot after = snapshot(g);
    assert(before.vertices == after.vertices && before.edges == after.edges);
    assert(before.bfs_count == after.bfs_count && before.comps == after.comps);

    /* Opting out mid-flight freezes it where it is. */
    graph_set_adaptive(g, false);
    const size_t migrations = graph_migration_count(g);
    for (int i = 0; i < 20000; ++i) {
        (void)graph_neighbors(g, i % 120, NULL, NULL, 0);
    }
    assert(graph_migration_count(g) == migrations);
    assert(!graph_is_adaptive(g));
    graph_destroy(g);
}

static void test_compare_does_not_touch_the_graph(void)
{
    Graph *g = make(false, GRAPH_LIST, true, 7);
    fill_sample(g);
    const Snapshot before = snapshot(g);

    GraphCompareResult all[GRAPH_ALGO_COUNT];
    assert(graph_compare_all(g, 0, 3, all));
    for (int a = 0; a < GRAPH_ALGO_COUNT; ++a) {
        assert(all[a].agreed);
        assert(all[a].list_ns > 0 && all[a].matrix_ns > 0);
        assert(all[a].list_bytes > 0 && all[a].matrix_bytes > 0);
        assert(strlen(graph_algo_name((GraphAlgo)a)) > 1);
    }
    assert(all[GRAPH_ALGO_BFS].list_ok && all[GRAPH_ALGO_BFS].matrix_ok);
    assert(!all[GRAPH_ALGO_TOPO].list_ok && !all[GRAPH_ALGO_TOPO].matrix_ok); /* undirected */
    assert(!all[GRAPH_ALGO_MST].list_ok);                                    /* disconnected */

    assert(graph_representation(g) == GRAPH_LIST);
    assert(graph_migration_count(g) == 0);
    const Snapshot after = snapshot(g);
    assert_same(&before, &after);

    /* A dead source falls back to the first alive vertex. */
    assert(graph_remove_vertex(g, 0));
    GraphCompareResult r;
    assert(graph_compare(g, GRAPH_ALGO_BFS, 0, 1, &r) && r.source == 1 && r.agreed);

    FILE *sink = fopen("/dev/null", "w");
    if (sink != NULL) {
        graph_compare_print(&r, sink);
        fclose(sink);
    }

    Graph *empty = make(false, GRAPH_LIST, false, 0);
    assert(!graph_compare(empty, GRAPH_ALGO_BFS, 0, 1, &r));
    graph_destroy(empty);
    graph_destroy(g);
}

static void test_compare_on_a_dag(void)
{
    Graph *g = make(true, GRAPH_MATRIX, false, 5);
    assert(graph_add_edge(g, 0, 1, 1) && graph_add_edge(g, 1, 2, 2));
    assert(graph_add_edge(g, 0, 3, 7) && graph_add_edge(g, 3, 4, 1));
    GraphCompareResult r;
    assert(graph_compare(g, GRAPH_ALGO_TOPO, 0, 2, &r) && r.list_ok && r.matrix_ok && r.agreed);
    assert(graph_compare(g, GRAPH_ALGO_DIJKSTRA, 0, 2, &r) && r.agreed);
    assert(graph_compare(g, GRAPH_ALGO_CYCLE, 0, 2, &r) && r.agreed);
    graph_destroy(g);
}

static void test_null_and_bad_arguments(void)
{
    graph_destroy(NULL);
    assert(graph_vertex_count(NULL) == 0);
    assert(!graph_add_edge(NULL, 0, 1, 1));
    assert(graph_add_vertex(NULL) == -1);
    assert(!graph_remove_vertex(NULL, 0));
    assert(!graph_convert(NULL, GRAPH_LIST));
    Graph *g = graph_create(false);
    assert(g != NULL && graph_vertex_count(g) == 0);
    assert(!graph_remove_vertex(g, 0));
    assert(!graph_add_edge(g, 0, 1, 1));
    assert(graph_add_vertex(g) == 0 && graph_add_vertex(g) == 1);
    assert(!graph_add_edge(g, 0, 1, INFINITY));
    assert(graph_add_edge(g, 0, 1, 1.0));
    graph_destroy(g);
}

static void test_named_graph(void)
{
    NamedGraph *ng = named_graph_create(false, KEYOPS_STRING);
    assert(ng != NULL);
    assert(named_add_vertex(ng, "Paris") == 0);
    assert(named_add_vertex(ng, "Rome") == 1);
    assert(named_add_edge(ng, "Paris", "Rome", 1100.0));
    assert(named_add_edge(ng, "Paris", "Oslo", 1800.0));     /* creates Oslo */
    assert(graph_edge_count(ng->g) == 2 && graph_vertex_count(ng->g) == 3);
    assert(strcmp(named_vertex_data(ng, 2), "Oslo") == 0);
    named_graph_destroy(ng);
}

int main(void)
{
    test_start_on_either_representation();
    test_same_answers_on_both();
    test_manual_convert_preserves_everything();
    test_remove_vertex_and_reuse_and_compact();
    test_neighbors();
    test_directed_algorithms();
    test_cache_never_goes_stale();
    test_cache_stats();
    test_adaptive_off_never_migrates();
    test_adaptive_on_migrates_to_matrix_for_point_queries();
    test_compare_does_not_touch_the_graph();
    test_compare_on_a_dag();
    test_null_and_bad_arguments();
    test_named_graph();
    printf("All unified-graph tests passed.\n");
    return 0;
}
