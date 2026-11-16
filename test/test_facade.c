/*
 * test_facade.c — the unified API beyond the basics:
 *   1. error codes recorded by the facade,
 *   2. list-vs-matrix cross-check of the algorithm wrappers on random graphs,
 *   3. observer pricing of the new workload kinds,
 *   4. frozen (read-only, thread-shareable) handles,
 *   5. NamedGraph vertex removal and compaction,
 *   6. graph-building wrappers (generators, transpose, union, readers),
 *   7. the list-vs-matrix crossover report.
 */

#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "backend.h"
#include "graph.h"
#include "graph_algorithms.h"
#include "graph_internal.h"      /* to read the observer counters */
#include "list_backend.h"
#include "mapping.h"
#include "matrix_backend.h"

#define CHECK(cond)                                                                  \
    do {                                                                             \
        if (!(cond)) {                                                               \
            fprintf(stderr, "%s:%d: CHECK failed: %s\n", __FILE__, __LINE__, #cond); \
            exit(1);                                                                 \
        }                                                                            \
    } while (0)

static unsigned long long rng_state = 88172645463325252ULL;
static unsigned rnd(unsigned bound)
{
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return (unsigned)(rng_state % bound);
}

static bool near(double a, double b)
{
    if (isinf(a) || isinf(b)) {
        return isinf(a) && isinf(b) && (a > 0) == (b > 0);
    }
    return fabs(a - b) <= 1e-9 * fmax(1.0, fmax(fabs(a), fabs(b)));
}

/* ------------------------------------------------------------------ */
/*  1. Error codes                                                      */
/* ------------------------------------------------------------------ */

#define EXPECT_ERR(expr, code) \
    do { graph_clear_error(); (void)(expr); CHECK(graph_last_error() == (code)); } while (0)

static void test_error_codes(void)
{
    Graph *g = graph_create_sized(4, false);
    Graph *d = graph_create_sized(3, true);
    CHECK(g != NULL && d != NULL);
    CHECK(graph_add_edge(g, 0, 1, 1.0));

    EXPECT_ERR(graph_add_edge(g, 0, 1, 1.0), GRAPH_ERR_DUPLICATE_EDGE);
    EXPECT_ERR(graph_add_edge(g, 0, 99, 1.0), GRAPH_ERR_INVALID_VERTEX);
    EXPECT_ERR(graph_add_edge(g, 1, 2, NAN), GRAPH_ERR_INVALID_WEIGHT);
    EXPECT_ERR(graph_remove_edge(g, 2, 3), GRAPH_ERR_EDGE_NOT_FOUND);
    EXPECT_ERR(graph_remove_vertex(g, 77), GRAPH_ERR_INVALID_VERTEX);
    EXPECT_ERR(graph_add_edge(NULL, 0, 1, 1.0), GRAPH_ERR_INVALID_ARGUMENT);
    EXPECT_ERR(graph_add_vertex(NULL), GRAPH_ERR_INVALID_ARGUMENT);

    int order[8], parent[8];
    double dist[8];
    EXPECT_ERR(graph_bfs(g, 50, order, parent), GRAPH_ERR_INVALID_VERTEX);
    EXPECT_ERR(graph_bfs(NULL, 0, order, parent), GRAPH_ERR_INVALID_ARGUMENT);
    EXPECT_ERR(graph_topological_sort(g, order), GRAPH_ERR_UNSUPPORTED);
    CHECK(graph_add_edge(d, 0, 1, 1.0) && graph_add_edge(d, 1, 2, 1.0) && graph_add_edge(d, 2, 0, 1.0));
    EXPECT_ERR(graph_topological_sort(d, order), GRAPH_ERR_CYCLE);
    EXPECT_ERR(graph_topological_sort(d, order), GRAPH_ERR_CYCLE);   /* served from the cache */
    EXPECT_ERR(graph_dijkstra(g, 0, NULL, parent), GRAPH_ERR_INVALID_ARGUMENT);

    CHECK(graph_add_edge(g, 2, 3, -1.0));
    EXPECT_ERR(graph_dijkstra(g, 0, dist, parent), GRAPH_ERR_NEGATIVE_WEIGHT);

    GraphEdge edges[8];
    size_t ne;
    CHECK(isinf(graph_mst_kruskal(g, edges, &ne)));            /* {0,1} and {2,3} are separate */
    CHECK(graph_last_error() == GRAPH_ERR_NOT_FOUND);
    graph_clear_error();
    CHECK(isinf(graph_mst_kruskal(g, edges, &ne)));            /* cached failure still reports */
    CHECK(graph_last_error() == GRAPH_ERR_NOT_FOUND);

    /* Library wrappers return the status directly. */
    long long tri;
    CHECK(graph_triangle_count(NULL, &tri) == GRAPH_ERR_INVALID_ARGUMENT);
    CHECK(graph_triangle_count(d, &tri) == GRAPH_ERR_UNSUPPORTED);
    double bf[8];
    CHECK(graph_bellman_ford(g, 0, bf, parent) == GRAPH_OK);

    graph_destroy(g);
    graph_destroy(d);
}

/* ------------------------------------------------------------------ */
/*  2. Cross-check of both representations                              */
/* ------------------------------------------------------------------ */

#define MAXN 32

typedef struct {
    GraphStatus st[48];
    double      num[48][MAXN];
    int         ints[16][MAXN];
    double      scalar[48];
    long long   big;
    size_t      count[16];
    bool        flag[16];
    bool        bools[4][MAXN];
    int         misc[8];
    size_t      common;
} Results;

static void gather(Graph *g, Results *r)
{
    memset(r, 0, sizeof(*r));
    const size_t n = graph_id_bound(g);
    int s = -1, t = -1;
    for (size_t v = 0; v < n; ++v) {
        if (graph_is_alive(g, (int)v)) {
            if (s < 0) { s = (int)v; }
            t = (int)v;
        }
    }
    int k = 0;
    r->st[k++] = graph_degree_centrality(g, r->num[0]);
    r->st[k++] = graph_closeness_centrality(g, r->num[1]);
    r->st[k++] = graph_betweenness_centrality(g, r->num[2]);
    r->st[k++] = graph_triangle_count(g, &r->big);
    r->st[k++] = graph_clustering_coefficient(g, r->num[3], &r->scalar[0]);
    r->st[k++] = graph_is_bipartite(g, &r->flag[0], r->ints[0]);
    r->st[k++] = graph_eccentricity(g, r->ints[1]);
    r->st[k++] = graph_diameter_radius(g, &r->misc[0], &r->misc[1]);
    r->flag[1] = graph_is_connected(g);
    r->flag[2] = graph_is_strongly_connected(g);
    r->flag[3] = graph_is_dag(g);
    r->flag[4] = graph_is_tree(g);
    r->flag[5] = graph_is_forest(g);
    r->st[k++] = graph_strongly_connected_components(g, r->ints[2], &r->count[0]);
    r->st[k++] = graph_articulation_points(g, r->bools[0], &r->count[1]);
    GraphEdge br[MAXN * MAXN];
    r->st[k++] = graph_bridges(g, br, &r->count[2]);
    r->st[k++] = graph_bellman_ford(g, s, r->num[4], r->ints[4]);
    static double fw[MAXN * MAXN];
    static int    nx[MAXN * MAXN];
    r->st[k++] = graph_floyd_warshall(g, fw, nx);
    r->scalar[3] = fw[(size_t)s * n + (size_t)t];
    GraphEdge mst[MAXN];
    r->scalar[4] = graph_mst_prim(g, mst, &r->count[3]);
    r->st[k++] = graph_max_flow(g, s, t, GRAPH_FLOW_DINIC, &r->scalar[5], NULL);
    r->st[k++] = graph_max_flow(g, s, t, GRAPH_FLOW_EDMONDS_KARP, &r->scalar[6], NULL);
    r->st[k++] = graph_bipartite_matching(g, r->ints[5], &r->count[4]);
    r->st[k++] = graph_color_greedy(g, r->ints[6], &r->count[5]);
    r->st[k++] = graph_color_welsh_powell(g, r->ints[7], &r->count[6]);
    r->st[k++] = graph_color_dsatur(g, r->ints[8], &r->count[7]);
    r->st[k++] = graph_max_clique(g, r->ints[9], &r->count[8]);
    GraphGroups cl;
    r->st[k] = graph_maximal_cliques(g, &cl);
    r->count[9] = r->st[k] == GRAPH_OK ? cl.count : 0;
    if (r->st[k] == GRAPH_OK) {
        graph_groups_free(&cl);
    }
    ++k;
    r->st[k++] = graph_pagerank(g, 0.85, 1e-12, 500, r->num[5], &r->count[10]);
    r->st[k++] = graph_core_numbers(g, r->ints[10]);
    r->st[k++] = graph_label_propagation(g, r->ints[11], &r->count[11]);
    r->st[k++] = graph_louvain(g, r->ints[12], &r->scalar[7], &r->count[12]);
    r->st[k++] = graph_degree_assortativity(g, &r->scalar[8]);
    r->st[k++] = graph_vertex_cover_approx(g, r->bools[1], &r->count[13]);
    r->st[k++] = graph_independent_set_greedy(g, r->bools[2], &r->count[14]);
    GraphNetworkStats ns;
    r->st[k] = graph_network_stats(g, &ns);
    r->scalar[9] = ns.density;
    r->scalar[10] = ns.reciprocity;
    r->count[15] = ns.max_degree;
    ++k;
    r->st[k++] = graph_jaccard(g, s, t, &r->scalar[11]);
    int cn[MAXN];
    r->st[k++] = graph_common_neighbors(g, s, t, cn, &r->common);
    r->st[k++] = graph_dag_longest_path(g, r->num[6], r->ints[15], &r->misc[2]);
}

static void compare_results(const Results *a, const Results *b, size_t n)
{
    for (int i = 0; i < 40; ++i) {
        CHECK(a->st[i] == b->st[i]);
    }
    /* Only compare payloads of calls that succeeded (st is identical on both). */
    #define OK(i) (a->st[i] == GRAPH_OK)
    if (OK(0)) { for (size_t v = 0; v < n; ++v) { CHECK(near(a->num[0][v], b->num[0][v])); } }
    if (OK(1)) { for (size_t v = 0; v < n; ++v) { CHECK(near(a->num[1][v], b->num[1][v])); } }
    if (OK(2)) { for (size_t v = 0; v < n; ++v) { CHECK(near(a->num[2][v], b->num[2][v])); } }
    if (OK(3)) { CHECK(a->big == b->big); }
    if (OK(4)) { for (size_t v = 0; v < n; ++v) { CHECK(near(a->num[3][v], b->num[3][v])); }
                 CHECK(near(a->scalar[0], b->scalar[0])); }
    if (OK(5)) { CHECK(a->flag[0] == b->flag[0]); }
    if (OK(6)) { for (size_t v = 0; v < n; ++v) { CHECK(a->ints[1][v] == b->ints[1][v]); } }
    if (OK(7)) { CHECK(a->misc[0] == b->misc[0] && a->misc[1] == b->misc[1]); }
    for (int i = 1; i <= 5; ++i) { CHECK(a->flag[i] == b->flag[i]); }
    if (OK(9)) { for (size_t v = 0; v < n; ++v) { CHECK(a->ints[2][v] == b->ints[2][v]); }
                 CHECK(a->count[0] == b->count[0]); }
    if (OK(10)) { CHECK(a->count[1] == b->count[1]); }
    if (OK(11)) { CHECK(a->count[2] == b->count[2]); }
    if (OK(12)) { for (size_t v = 0; v < n; ++v) { CHECK(near(a->num[4][v], b->num[4][v])); } }
    if (OK(13)) { CHECK(near(a->scalar[3], b->scalar[3])); }
    CHECK(near(a->scalar[4], b->scalar[4]));
    CHECK(a->count[3] == b->count[3]);
    if (OK(14)) { CHECK(near(a->scalar[5], b->scalar[5])); }
    if (OK(15)) { CHECK(near(a->scalar[6], b->scalar[6])); }
    if (OK(16)) { CHECK(a->count[4] == b->count[4]); }
    if (OK(17)) { CHECK(a->count[5] == b->count[5]); }
    if (OK(18)) { CHECK(a->count[6] == b->count[6]); }
    if (OK(19)) { CHECK(a->count[7] == b->count[7]); }
    if (OK(20)) { CHECK(a->count[8] == b->count[8]); }
    if (OK(21)) { CHECK(a->count[9] == b->count[9]); }
    if (OK(22)) { for (size_t v = 0; v < n; ++v) { CHECK(near(a->num[5][v], b->num[5][v])); } }
    if (OK(23)) { for (size_t v = 0; v < n; ++v) { CHECK(a->ints[10][v] == b->ints[10][v]); } }
    if (OK(24)) { CHECK(a->count[11] == b->count[11]); }
    if (OK(25)) { CHECK(near(a->scalar[7], b->scalar[7])); CHECK(a->count[12] == b->count[12]); }
    if (OK(26)) { CHECK(near(a->scalar[8], b->scalar[8])); }
    if (OK(27)) { CHECK(a->count[13] == b->count[13]); }
    if (OK(28)) { CHECK(a->count[14] == b->count[14]); }
    if (OK(29)) { CHECK(near(a->scalar[9], b->scalar[9])); CHECK(near(a->scalar[10], b->scalar[10]));
                  CHECK(a->count[15] == b->count[15]); }
    if (OK(30)) { CHECK(near(a->scalar[11], b->scalar[11])); }
    if (OK(31)) { CHECK(a->common == b->common); }
    #undef OK
}

static void test_cross_check(void)
{
    unsigned compared = 0, ok_calls = 0;
    for (unsigned round = 0; round < 80; ++round) {
        const bool   directed = (round & 1) != 0;
        const size_t n = 6 + rnd(14);
        const double p = 0.08 + 0.1 * (double)rnd(6);
        Graph *g = graph_gen_erdos_renyi(round % 3 == 0 ? GRAPH_MATRIX : GRAPH_LIST, n, p,
                                         1000 + round, directed, 1.0, 9.0);
        CHECK(g != NULL);
        graph_set_adaptive(g, false);
        for (unsigned k = 0; k < rnd(3); ++k) {
            (void)graph_remove_vertex(g, (int)rnd((unsigned)n));
        }
        if (graph_vertex_count(g) < 3) {
            graph_destroy(g);
            continue;
        }
        Results on_list, on_matrix;
        CHECK(graph_convert(g, GRAPH_LIST));
        gather(g, &on_list);
        CHECK(graph_convert(g, GRAPH_MATRIX));
        gather(g, &on_matrix);
        compare_results(&on_list, &on_matrix, graph_id_bound(g));
        for (int i = 0; i < 40; ++i) {
            ok_calls += on_list.st[i] == GRAPH_OK;
        }
        ++compared;
        graph_destroy(g);
    }
    CHECK(compared >= 60);
    CHECK(ok_calls > 60 * 15);     /* the comparison was not vacuous */
}

/* ------------------------------------------------------------------ */
/*  3. Observer pricing                                                 */
/* ------------------------------------------------------------------ */

static double cost_of(const GraphBackend *be, size_t n, double density, OpKind op)
{
    GraphShape shape = { n, n, (size_t)(density * (double)n * (double)(n - 1) / 2.0), false, false};
    Workload w;
    workload_reset(&w);
    w.n[op] = 1.0;
    return be->estimate_cost(&shape, &w);
}

static void test_cost_model(void)
{
    const GraphBackend *l = list_backend(), *m = matrix_backend();
    /* Pair probing is cheap on a dense matrix, sweeping over edges is cheap on a list. */
    CHECK(cost_of(m, 2000, 0.2, OP_ADJACENCY_RUN) < cost_of(l, 2000, 0.2, OP_ADJACENCY_RUN));
    CHECK(cost_of(l, 2000, 0.01, OP_ITERATIVE_RUN) < cost_of(m, 2000, 0.01, OP_ITERATIVE_RUN));
    CHECK(cost_of(l, 2000, 0.2, OP_ITERATIVE_RUN) < cost_of(m, 2000, 0.2, OP_ITERATIVE_RUN));
    /* Costs grow with the graph. */
    for (int op = OP_ITERATIVE_RUN; op <= OP_ADJACENCY_RUN; ++op) {
        CHECK(cost_of(l, 1000, 0.1, (OpKind)op) > cost_of(l, 100, 0.1, (OpKind)op));
        CHECK(cost_of(m, 1000, 0.1, (OpKind)op) > cost_of(m, 100, 0.1, (OpKind)op));
        CHECK(cost_of(l, 100, 0.1, (OpKind)op) > 0.0 && cost_of(m, 100, 0.1, (OpKind)op) > 0.0);
    }

    /* The wrappers feed those kinds to the observer: a probing-heavy workload
     * on a dense list graph makes the observer move to the matrix. */
    GraphConfig cfg = graph_default_config(false);
    cfg.observer.check_interval = 8;
    cfg.observer.cooldown_ops   = 0;
    cfg.observer.hysteresis     = 1.1;
    Graph *g = graph_create_with(&cfg, 0);
    CHECK(g != NULL);
    Graph *src = graph_gen_erdos_renyi(GRAPH_LIST, 600, 0.3, 5, false, 1.0, 5.0);
    CHECK(src != NULL);
    for (int v = 0; v < 600; ++v) { CHECK(graph_add_vertex(g) == v); }
    for (int u = 0; u < 600; ++u) {
        for (int v = u + 1; v < 600; ++v) {
            if (graph_has_edge(src, u, v, NULL)) { CHECK(graph_add_edge(g, u, v, 1.0)); }
        }
    }
    graph_set_adaptive(g, false);
    CHECK(graph_representation(g) == GRAPH_LIST);
    graph_set_adaptive(g, true);
    long long tri;
    for (int i = 0; i < 40; ++i) {
        CHECK(graph_triangle_count(g, &tri) == GRAPH_OK);
    }
    CHECK(graph_representation(g) == GRAPH_MATRIX);
    CHECK(graph_migration_count(g) >= 1);
    graph_destroy(g);
    graph_destroy(src);
}

/* ------------------------------------------------------------------ */
/*  4. Frozen handles                                                   */
/* ------------------------------------------------------------------ */

typedef struct {
    Graph  *g;
    double  pagerank_sum;
    long long triangles;
    size_t  visited;
    bool    ok;
} Job;

static void *reader(void *arg)
{
    Job *job = arg;
    job->ok = true;
    const size_t n = graph_id_bound(job->g);
    double *pr = malloc(n * sizeof(double));
    int    *order = malloc(n * sizeof(int));
    for (int rep = 0; rep < 20; ++rep) {
        job->ok = job->ok && graph_pagerank(job->g, 0.85, 1e-10, 100, pr, NULL) == GRAPH_OK;
        job->ok = job->ok && graph_triangle_count(job->g, &job->triangles) == GRAPH_OK;
        job->visited = graph_bfs(job->g, 0, order, NULL);
        job->ok = job->ok && graph_dijkstra(job->g, 0, pr, NULL);
        graph_clear_error();
        CHECK(!graph_add_edge(job->g, 0, 1, 1.0));          /* mutation refused */
        job->ok = job->ok && graph_last_error() == GRAPH_ERR_FROZEN;
    }
    double s = 0.0;
    CHECK(graph_pagerank(job->g, 0.85, 1e-10, 100, pr, NULL) == GRAPH_OK);
    for (size_t v = 0; v < n; ++v) { s += pr[v]; }
    job->pagerank_sum = s;
    free(pr);
    free(order);
    return NULL;
}

static void test_frozen(void)
{
    Graph *g = graph_gen_erdos_renyi(GRAPH_LIST, 120, 0.1, 9, false, 1.0, 5.0);
    CHECK(g != NULL);
    CHECK(!graph_is_frozen(g));
    graph_freeze(g);
    CHECK(graph_is_frozen(g));
    const size_t edges = graph_edge_count(g);

    EXPECT_ERR(graph_add_vertex(g), GRAPH_ERR_FROZEN);
    EXPECT_ERR(graph_remove_vertex(g, 0), GRAPH_ERR_FROZEN);
    EXPECT_ERR(graph_remove_edge(g, 0, 1), GRAPH_ERR_FROZEN);
    EXPECT_ERR(graph_compact(g, NULL), GRAPH_ERR_FROZEN);
    EXPECT_ERR(graph_convert(g, GRAPH_MATRIX), GRAPH_ERR_FROZEN);
    CHECK(graph_representation(g) == GRAPH_LIST && graph_edge_count(g) == edges);

    const double ops_before     = workload_graph_ops(observer_workload(g->obs));
    const double lookups_before = observer_workload(g->obs)->n[OP_CACHE_HIT] +
                                  observer_workload(g->obs)->n[OP_CACHE_MISS];
    enum { THREADS = 6 };
    pthread_t th[THREADS];
    Job jobs[THREADS];
    for (int i = 0; i < THREADS; ++i) {
        jobs[i] = (Job){ .g = g };
        CHECK(pthread_create(&th[i], NULL, reader, &jobs[i]) == 0);
    }
    for (int i = 0; i < THREADS; ++i) {
        CHECK(pthread_join(th[i], NULL) == 0);
    }
    for (int i = 0; i < THREADS; ++i) {
        CHECK(jobs[i].ok);
        CHECK(near(jobs[i].pagerank_sum, jobs[0].pagerank_sum));
        CHECK(jobs[i].triangles == jobs[0].triangles && jobs[i].visited == jobs[0].visited);
    }
    CHECK(graph_migration_count(g) == 0);       /* nothing was recorded or migrated */
    CHECK(workload_graph_ops(observer_workload(g->obs)) == ops_before);   /* no counters touched */
    CHECK(observer_workload(g->obs)->n[OP_CACHE_HIT] +
          observer_workload(g->obs)->n[OP_CACHE_MISS] == lookups_before); /* cache not consulted */

    graph_thaw(g);
    CHECK(!graph_is_frozen(g));
    CHECK(graph_add_vertex(g) >= 0);
    graph_destroy(g);
}

/* ------------------------------------------------------------------ */
/*  5. NamedGraph removal                                               */
/* ------------------------------------------------------------------ */

static void test_named_removal(void)
{
    NamedGraph *ng = named_graph_create(false, KEYOPS_STRING);
    CHECK(ng != NULL);
    CHECK(named_add_edge(ng, "ann", "bob", 1.0));
    CHECK(named_add_edge(ng, "bob", "cy", 2.0));
    CHECK(named_add_edge(ng, "cy", "dee", 3.0));
    CHECK(named_add_edge(ng, "dee", "ann", 4.0));
    CHECK(graph_vertex_count(ng->g) == 4 && graph_edge_count(ng->g) == 4);
    const int bob = named_vertex_id(ng, "bob");
    CHECK(bob >= 0);

    CHECK(named_remove_vertex(ng, "bob"));
    CHECK(named_vertex_id(ng, "bob") < 0);
    CHECK(named_vertex_data(ng, bob) == NULL);
    CHECK(graph_vertex_count(ng->g) == 3 && graph_edge_count(ng->g) == 2);
    CHECK(!graph_is_alive(ng->g, bob));
    graph_clear_error();
    CHECK(!named_remove_vertex(ng, "bob"));
    CHECK(graph_last_error() == GRAPH_ERR_NOT_FOUND);
    CHECK(!named_remove_vertex(NULL, "x") && graph_last_error() == GRAPH_ERR_INVALID_ARGUMENT);

    /* The next vertex recycles bob's ID, and the name lookup follows it. */
    const int eve = named_add_vertex(ng, "eve");
    CHECK(eve == bob);
    CHECK(strcmp(named_vertex_data(ng, eve), "eve") == 0);
    CHECK(named_vertex_id(ng, "eve") == eve);
    CHECK(named_add_edge(ng, "eve", "cy", 5.0));
    double w;
    CHECK(graph_has_edge(ng->g, eve, named_vertex_id(ng, "cy"), &w) && w == 5.0);

    CHECK(named_remove_edge(ng, "eve", "cy"));
    CHECK(!named_remove_edge(ng, "eve", "cy"));
    CHECK(graph_last_error() == GRAPH_ERR_EDGE_NOT_FOUND);
    CHECK(!named_remove_edge(ng, "eve", "nobody") && graph_last_error() == GRAPH_ERR_NOT_FOUND);

    /* Remove several, compact, and check every surviving name still maps to its vertex. */
    char name[16];
    for (int i = 0; i < 40; ++i) {
        snprintf(name, sizeof(name), "v%d", i);
        CHECK(named_add_edge(ng, name, "ann", (double)(i + 1)));
    }
    for (int i = 0; i < 40; i += 3) {
        snprintf(name, sizeof(name), "v%d", i);
        CHECK(named_remove_vertex(ng, name));
    }
    CHECK(named_remove_vertex(ng, "ann"));
    CHECK(named_compact(ng));
    CHECK(graph_vertex_count(ng->g) == graph_id_bound(ng->g));
    CHECK(mapping_count(ng->map) == graph_vertex_count(ng->g));
    CHECK(mapping_id_bound(ng->map) == mapping_count(ng->map));
    for (int i = 0; i < 40; ++i) {
        snprintf(name, sizeof(name), "v%d", i);
        const int id = named_vertex_id(ng, name);
        if (i % 3 == 0) {
            CHECK(id < 0);
        } else {
            CHECK(id >= 0 && graph_is_alive(ng->g, id));
            CHECK(strcmp(named_vertex_data(ng, id), name) == 0);
        }
    }
    CHECK(named_vertex_id(ng, "ann") < 0);
    CHECK(named_vertex_id(ng, "dee") >= 0);

    /* Mapping-level behaviour: assign / remove / reuse. */
    Mapping *m = mapping_create(KEYOPS_INT);
    CHECK(m != NULL);
    int k0 = 10, k1 = 20, k2 = 30;
    CHECK(mapping_get_or_create(m, &k0) == 0 && mapping_get_or_create(m, &k1) == 1);
    CHECK(mapping_assign(m, &k2, 1) == -1);          /* ID taken */
    CHECK(mapping_assign(m, &k0, 2) == -1);          /* key present */
    CHECK(mapping_remove(m, 0) && !mapping_remove(m, 0));
    CHECK(mapping_find(m, &k0) < 0 && mapping_count(m) == 1);
    CHECK(mapping_assign(m, &k2, 0) == 0 && mapping_find(m, &k2) == 0);
    CHECK(mapping_find(m, &k1) == 1);
    mapping_destroy(m);
    named_graph_destroy(ng);
}

/* ------------------------------------------------------------------ */
/*  6. Graph-building wrappers                                          */
/* ------------------------------------------------------------------ */

static bool same_graph(Graph *a, Graph *b)
{
    if (graph_id_bound(a) != graph_id_bound(b) || graph_edge_count(a) != graph_edge_count(b)) {
        return false;
    }
    for (int u = 0; u < (int)graph_id_bound(a); ++u) {
        if (graph_is_alive(a, u) != graph_is_alive(b, u)) {
            return false;
        }
        for (int v = 0; v < (int)graph_id_bound(a); ++v) {
            double wa = 0, wb = 0;
            const bool ea = graph_has_edge(a, u, v, &wa), eb = graph_has_edge(b, u, v, &wb);
            if (ea != eb || (ea && wa != wb)) {
                return false;
            }
        }
    }
    return true;
}

static void test_builders(void)
{
    for (int rep = 0; rep < 2; ++rep) {
        const GraphRep r = rep == 0 ? GRAPH_LIST : GRAPH_MATRIX;
        Graph *g = graph_gen_erdos_renyi(r, 25, 0.2, 77, true, 1.0, 9.0);
        CHECK(g != NULL && graph_representation(g) == r);
        graph_set_adaptive(g, false);

        Graph *copy = graph_copy(g);
        CHECK(copy != NULL && same_graph(g, copy) && graph_representation(copy) == r);

        Graph *t = NULL, *tt = NULL;
        CHECK(graph_transpose(g, &t) == GRAPH_OK && t != NULL);
        CHECK(graph_edge_count(t) == graph_edge_count(g));
        CHECK(graph_transpose(t, &tt) == GRAPH_OK && same_graph(g, tt));

        Graph *other = graph_gen_erdos_renyi(r == GRAPH_LIST ? GRAPH_MATRIX : GRAPH_LIST,
                                             25, 0.2, 78, true, 1.0, 9.0);
        Graph *u = NULL, *i = NULL;
        CHECK(graph_union(g, other, &u) == GRAPH_OK);       /* crosses representations */
        CHECK(graph_intersection(g, other, &i) == GRAPH_OK);
        CHECK(graph_edge_count(u) + graph_edge_count(i) == graph_edge_count(g) + graph_edge_count(other));
        bool iso = false;
        CHECK(graph_is_isomorphic(g, copy, &iso, NULL) == GRAPH_OK && iso);

        Graph *comp = NULL;
        CHECK(graph_complement(g, &comp) == GRAPH_OK);
        CHECK(graph_edge_count(comp) + graph_edge_count(g) == 25u * 24u);

        int keep[] = { 0, 2, 4, 6 };
        Graph *sub = NULL;
        CHECK(graph_induced_subgraph(g, keep, 4, &sub) == GRAPH_OK && graph_vertex_count(sub) == 4);

        /* Failure leaves *out NULL and records the code. */
        Graph *bad = (Graph *)1;
        CHECK(graph_transpose(NULL, &bad) == GRAPH_ERR_INVALID_ARGUMENT && bad == NULL);

        graph_destroy(g); graph_destroy(copy); graph_destroy(t); graph_destroy(tt);
        graph_destroy(other); graph_destroy(u); graph_destroy(i); graph_destroy(comp);
        graph_destroy(sub);
    }

    /* Generators produce what they promise in either representation. */
    for (int rep = 0; rep < 2; ++rep) {
        const GraphRep r = rep == 0 ? GRAPH_LIST : GRAPH_MATRIX;
        Graph *p = graph_gen_path(r, 10, false);
        Graph *c = graph_gen_cycle(r, 10, false);
        Graph *k = graph_gen_complete(r, 8, false);
        Graph *tr = graph_gen_random_tree(r, 30, 5);
        CHECK(p && c && k && tr);
        CHECK(graph_edge_count(p) == 9 && graph_edge_count(c) == 10 && graph_edge_count(k) == 28);
        CHECK(graph_is_tree(tr) && graph_is_tree(p) && !graph_is_tree(c));
        graph_destroy(p); graph_destroy(c); graph_destroy(k); graph_destroy(tr);
    }
    graph_clear_error();
    CHECK(graph_gen_erdos_renyi(GRAPH_LIST, 10, 2.0, 1, false, 1.0, 2.0) == NULL);   /* p > 1 */
    CHECK(graph_last_error() != GRAPH_OK);

    /* Round trips through files and the binary format, into either representation. */
    Graph *g = graph_gen_erdos_renyi(GRAPH_LIST, 20, 0.2, 3, false, 1.0, 9.0);
    FILE *f = tmpfile();
    CHECK(g != NULL && f != NULL);
    CHECK(graph_save_binary(g, f) == GRAPH_OK);
    rewind(f);
    Graph *back = NULL;
    CHECK(graph_load_binary(f, GRAPH_MATRIX, &back) == GRAPH_OK);
    CHECK(graph_representation(back) == GRAPH_MATRIX && same_graph(g, back));
    fclose(f);
    f = tmpfile();
    CHECK(graph_write_edge_list(g, f) == GRAPH_OK);
    rewind(f);
    Graph *back2 = NULL;
    CHECK(graph_read_edge_list(f, false, GRAPH_LIST, &back2) == GRAPH_OK);
    CHECK(graph_edge_count(back2) == graph_edge_count(g));
    fclose(f);
    f = tmpfile();
    fputs("not a graph\n", f);
    rewind(f);
    Graph *none = (Graph *)1;
    CHECK(graph_load_binary(f, GRAPH_LIST, &none) == GRAPH_ERR_IO && none == NULL);
    fclose(f);
    graph_destroy(g); graph_destroy(back); graph_destroy(back2);

    double cost[9] = { 4, 1, 3, 2, 0, 5, 3, 2, 2 };
    int assign[3];
    double total;
    CHECK(graph_hungarian(cost, 3, 3, assign, &total) == GRAPH_OK && near(total, 5.0));
}

/* ------------------------------------------------------------------ */
/*  7. Crossover report                                                 */
/* ------------------------------------------------------------------ */

static void test_crossover(void)
{
    const size_t sizes[]      = { 60, 150 };
    const double densities[]  = { 0.02, 0.2, 0.6 };
    FILE *f = tmpfile();
    CHECK(f != NULL);
    CHECK(graph_crossover_report(f, sizes, 2, densities, 3, 1, 1, false));
    long len = ftell(f);
    CHECK(len > 500);
    rewind(f);
    char *text = calloc((size_t)len + 1, 1);
    CHECK(text != NULL && fread(text, 1, (size_t)len, f) == (size_t)len);
    CHECK(strstr(text, "Crossover summary") != NULL);
    CHECK(strstr(text, "PageRank") != NULL && strstr(text, "Triangle count") != NULL);
    CHECK(strstr(text, "RESULTS DIFFER") == NULL);           /* both representations agreed */
    free(text);
    fclose(f);

    f = tmpfile();
    CHECK(graph_crossover_report(f, sizes, 1, densities, 2, 1, 1, true));
    rewind(f);
    char line[256];
    CHECK(fgets(line, sizeof line, f) && strncmp(line, "algorithm,n,density", 19) == 0);
    fclose(f);

    graph_clear_error();
    CHECK(!graph_crossover_report(NULL, sizes, 1, densities, 1, 1, 1, false));
    CHECK(graph_last_error() == GRAPH_ERR_INVALID_ARGUMENT);
}

int main(void)
{
    test_error_codes();
    test_cross_check();
    test_cost_model();
    test_frozen();
    test_named_removal();
    test_builders();
    test_crossover();
    puts("All facade tests passed.");
    return 0;
}
