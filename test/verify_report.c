/*
 * verify_report.c — a verification report: for each thing the library does, print what the answer
 * SHOULD be next to what the library GAVE, and whether they match.
 *
 *   make verify        # prints the whole table, exits non-zero if any row fails
 *
 * Where the expected values come from (never from the library under test):
 *   - closed forms for graph families of any size (a path of n vertices has n-1 edges and n-1 bridges;
 *     a k x k grid has diameter 2(k-1); K_n has C(n,3) triangles; ...);
 *   - textbook examples whose answers are published (the max-flow network of CLRS, the Kruskal example
 *     of the literature, the Petersen graph's famous lack of a Hamiltonian cycle, ...), re-checked by a
 *     separate brute-force script;
 *   - for the list and matrix representations: each must give the same answer.
 *
 * Sections:  1 graph families by size   2 textbook examples   3 advanced algorithms   4 social queries
 *            5 files too big for memory (streaming and sharding)   6 memory limits   7 list == matrix
 */

#define _POSIX_C_SOURCE 200809L

#include <math.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "graph_advanced.h"
#include "graph_algorithms.h"
#include "graph_shard.h"
#include "graph_social.h"
#include "graph_stream.h"

/* ------------------------------------------------------------------ */
/*  The table                                                           */
/* ------------------------------------------------------------------ */

static int g_rows, g_failed, g_section_rows, g_section_failed;
static bool g_color;

#define GREEN (g_color ? "\033[32m" : "")
#define RED   (g_color ? "\033[31m" : "")
#define BOLD  (g_color ? "\033[1m" : "")
#define DIM   (g_color ? "\033[2m" : "")
#define RESET (g_color ? "\033[0m" : "")

static void section(const char *title)
{
    if (g_section_rows > 0) {
        printf("  %s-> %d checks, %d failed%s\n", DIM, g_section_rows, g_section_failed, RESET);
    }
    g_section_rows = g_section_failed = 0;
    printf("\n%s=== %s%s\n", BOLD, title, RESET);
    printf("%s  %-5s %-58s %-34s %s%s\n", DIM, "", "check", "expected", "actual", RESET);
}

static void subsection(const char *title)
{
    printf("%s--- %s%s\n", DIM, title, RESET);
}

static void row(bool ok, const char *name, const char *expected, const char *actual)
{
    g_rows++; g_section_rows++;
    if (!ok) { g_failed++; g_section_failed++; }
    printf("  %s%-5s%s %-58.58s %-34.34s %s\n", ok ? GREEN : RED, ok ? "ok" : "FAIL", RESET, name, expected, actual);
}

static void check_size(const char *name, size_t expected, size_t actual)
{
    char e[40], a[40];
    snprintf(e, sizeof e, "%zu", expected);
    snprintf(a, sizeof a, "%zu", actual);
    row(expected == actual, name, e, a);
}

static void check_int(const char *name, long expected, long actual)
{
    char e[40], a[40];
    snprintf(e, sizeof e, "%ld", expected);
    snprintf(a, sizeof a, "%ld", actual);
    row(expected == actual, name, e, a);
}

static void check_num(const char *name, double expected, double actual)
{
    char e[40], a[40];
    snprintf(e, sizeof e, "%.6g", expected);
    snprintf(a, sizeof a, "%.6g", actual);
    row(fabs(expected - actual) <= 1e-9 * (1.0 + fabs(expected)), name, e, a);
}

static void check_bool(const char *name, bool expected, bool actual)
{
    row(expected == actual, name, expected ? "yes" : "no", actual ? "yes" : "no");
}

static void check_status(const char *name, GraphStatus expected, GraphStatus actual)
{
    row(expected == actual, name, graph_status_string(expected), graph_status_string(actual));
}

/* A property check: `expected` describes the property, `ok` is whether it held, `actual` how it looks. */
static void check_prop(const char *name, const char *expected, bool ok, const char *actual_fmt, ...)
{
    char a[80];
    va_list ap;
    va_start(ap, actual_fmt);
    vsnprintf(a, sizeof a, actual_fmt, ap);
    va_end(ap);
    row(ok, name, expected, a);
}

static void *xmalloc(size_t n)
{
    void *p = malloc(n > 0 ? n : 1);
    if (p == NULL) { fprintf(stderr, "out of memory\n"); exit(2); }
    return p;
}

static Graph *build(GraphRep rep, size_t n, bool directed, const double (*e)[3], size_t m)
{
    Graph *g = graph_create_sized(n, directed);
    if (g == NULL) { fprintf(stderr, "cannot create graph\n"); exit(2); }
    graph_set_adaptive(g, false);
    for (size_t i = 0; i < m; i++) graph_add_edge(g, (int)e[i][0], (int)e[i][1], e[i][2]);
    if (!graph_convert(g, rep)) { fprintf(stderr, "cannot convert graph\n"); exit(2); }
    return g;
}

static const char *rep_name(GraphRep r) { return r == GRAPH_MATRIX ? "matrix" : "list"; }

/* ------------------------------------------------------------------ */
/*  1. Graph families, by size                                          */
/* ------------------------------------------------------------------ */

static void family_path(GraphRep rep, size_t n, bool deep)
{
    char title[96], name[96];
    snprintf(title, sizeof title, "path of %zu vertices (%s)", n, rep_name(rep));
    subsection(title);
    Graph *g = graph_gen_path(rep, n, false);
    check_size("vertices", n, graph_vertex_count(g));
    check_size("edges (n-1)", n - 1, graph_edge_count(g));
    check_bool("is a tree", true, graph_is_tree(g));
    int *comp = xmalloc(n * sizeof *comp);
    check_size("connected components", 1, graph_connected_components(g, comp));
    free(comp);
    GraphEdge *br = xmalloc(n * sizeof *br);
    size_t nbr = 0, ncut = 0;
    bool *cut = xmalloc(n);
    graph_bridges(g, br, &nbr);
    check_size("bridges (every edge)", n - 1, nbr);
    graph_articulation_points(g, cut, &ncut);
    check_size("cut vertices (all but the two ends)", n - 2, ncut);
    double *dist = xmalloc(n * sizeof *dist);
    int *par = xmalloc(n * sizeof *par);
    graph_dijkstra(g, 0, dist, par);
    check_num("distance from one end to the other (n-1)", (double)(n - 1), dist[n - 1]);
    GraphEdge *mst = xmalloc(n * sizeof *mst);
    size_t mc = 0;
    check_num("minimum spanning tree weight (n-1)", (double)(n - 1), graph_mst_kruskal(g, mst, &mc));
    bool bip = false;
    graph_is_bipartite(g, &bip, NULL);
    check_bool("bipartite", true, bip);
    if (deep) {
        int diam = 0, rad = 0;
        graph_diameter_radius(g, &diam, &rad);
        snprintf(name, sizeof name, "diameter / radius (n-1 / %zu)", n / 2);
        check_prop(name, "diameter n-1, radius n/2", (size_t)diam == n - 1 && (size_t)rad == n / 2, "%d / %d", diam, rad);
    }
    free(br); free(cut); free(dist); free(par); free(mst);
    graph_destroy(g);
}

static void family_cycle(GraphRep rep, size_t n, bool deep)
{
    char title[96];
    snprintf(title, sizeof title, "cycle of %zu vertices (%s)", n, rep_name(rep));
    subsection(title);
    Graph *g = graph_gen_cycle(rep, n, false);
    check_size("edges (n)", n, graph_edge_count(g));
    check_bool("contains a cycle", true, graph_has_cycle(g));
    GraphEdge *br = xmalloc(n * sizeof *br);
    bool *cut = xmalloc(n);
    size_t nbr = 9, ncut = 9;
    graph_bridges(g, br, &nbr);
    graph_articulation_points(g, cut, &ncut);
    check_size("bridges", 0, nbr);
    check_size("cut vertices", 0, ncut);
    bool bip = false;
    graph_is_bipartite(g, &bip, NULL);
    check_bool(n % 2 == 0 ? "bipartite (even length)" : "bipartite (odd length: no)", n % 2 == 0, bip);
    double *dist = xmalloc(n * sizeof *dist);
    int *par = xmalloc(n * sizeof *par);
    graph_dijkstra(g, 0, dist, par);
    check_num("distance to the opposite vertex (n/2)", (double)(n / 2), dist[n / 2]);
    int *color = xmalloc(n * sizeof *color);
    size_t nc = 0;
    graph_color_dsatur(g, color, &nc);
    check_size(n % 2 == 0 ? "colours needed (even cycle)" : "colours needed (odd cycle)", n % 2 == 0 ? 2 : 3, nc);
    int *trail = xmalloc((n + 2) * sizeof *trail);
    size_t len = 0;
    check_status("Euler circuit exists", GRAPH_OK, graph_euler_trail(g, true, trail, &len));
    check_size("  ... and visits n+1 stops", n + 1, len);
    if (deep) {
        int diam = 0, rad = 0;
        graph_diameter_radius(g, &diam, &rad);
        check_int("diameter (n/2)", (long)(n / 2), diam);
    }
    free(br); free(cut); free(dist); free(par); free(color); free(trail);
    graph_destroy(g);
}

static void family_star(GraphRep rep, size_t n, bool deep)
{
    char title[96];
    snprintf(title, sizeof title, "star of %zu vertices (%s)", n, rep_name(rep));
    subsection(title);
    Graph *g = graph_gen_star(rep, n);
    check_size("edges (n-1)", n - 1, graph_edge_count(g));
    check_size("hub degree (n-1)", n - 1, graph_degree(g, 0));
    GraphNetworkStats st;
    graph_network_stats(g, &st);
    check_size("largest degree of anyone", n - 1, st.max_degree);
    bool *in = xmalloc(n);
    size_t sz = 0;
    graph_independent_set_greedy(g, in, &sz);
    check_size("independent set (all leaves, n-1)", n - 1, sz);
    if (deep) {   /* triangle counting costs the sum of squared degrees: a million-leaf hub is 10^12 */
        long long tri = -1;
        graph_triangle_count(g, &tri);
        check_int("triangles", 0, (long)tri);
    }
    size_t ncut = 0;
    bool *cut = xmalloc(n);
    graph_articulation_points(g, cut, &ncut);
    check_size("cut vertices (only the hub)", 1, ncut);
    if (deep) {
        int diam = 0, rad = 0;
        graph_diameter_radius(g, &diam, &rad);
        check_prop("diameter 2, radius 1", "2 / 1", diam == 2 && rad == 1, "%d / %d", diam, rad);
    }
    free(in); free(cut);
    graph_destroy(g);
}

static void family_grid(GraphRep rep, size_t s, bool deep)
{
    char title[96];
    snprintf(title, sizeof title, "grid %zu x %zu = %zu vertices (%s)", s, s, s * s, rep_name(rep));
    subsection(title);
    const size_t n = s * s;
    Graph *g = graph_gen_grid(rep, s, s);
    check_size("vertices (s*s)", n, graph_vertex_count(g));
    check_size("edges (2s(s-1))", 2 * s * (s - 1), graph_edge_count(g));
    bool bip = false;
    graph_is_bipartite(g, &bip, NULL);
    check_bool("bipartite (checkerboard)", true, bip);
    long long tri = -1;
    graph_triangle_count(g, &tri);
    check_int("triangles", 0, (long)tri);
    int *core = xmalloc(n * sizeof *core);
    graph_core_numbers(g, core);
    bool all2 = true;
    for (size_t v = 0; v < n; v++) if (core[v] != 2) all2 = false;
    check_bool("every vertex has core number 2", true, all2);
    int *match = xmalloc(n * sizeof *match);
    size_t msize = 0;
    graph_bipartite_matching(g, match, &msize);
    check_size("largest matching (floor(n/2))", n / 2, msize);
    int *color = xmalloc(n * sizeof *color);
    size_t nc = 0;
    graph_color_dsatur(g, color, &nc);
    check_size("colours (DSatur is exact on bipartite graphs)", 2, nc);
    double dist_expected = 2.0 * (double)(s - 1);
    double *dist = xmalloc(n * sizeof *dist);
    int *par = xmalloc(n * sizeof *par);
    graph_dijkstra(g, 0, dist, par);
    check_num("corner to corner distance (2(s-1))", dist_expected, dist[n - 1]);
    if (deep) {
        int diam = 0, rad = 0;
        graph_diameter_radius(g, &diam, &rad);
        check_int("diameter (2(s-1))", (long)dist_expected, diam);
    }
    double flow = 0;
    graph_max_flow(g, 0, (int)n - 1, GRAPH_FLOW_DINIC, &flow, NULL);
    check_num("max flow corner to corner (each corner has 2 edges)", 2.0, flow);
    free(core); free(match); free(color); free(dist); free(par);
    graph_destroy(g);
}

static void family_complete(GraphRep rep, size_t n)
{
    char title[96];
    snprintf(title, sizeof title, "complete graph K%zu (%s)", n, rep_name(rep));
    subsection(title);
    Graph *g = graph_gen_complete(rep, n, false);
    check_size("edges (n(n-1)/2)", n * (n - 1) / 2, graph_edge_count(g));
    long long tri = 0;
    graph_triangle_count(g, &tri);
    check_int("triangles (C(n,3))", (long)(n * (n - 1) * (n - 2) / 6), (long)tri);
    int *clique = xmalloc(n * sizeof *clique);
    size_t cs = 0;
    graph_max_clique(g, clique, &cs);
    check_size("largest clique (everyone)", n, cs);
    int *color = xmalloc(n * sizeof *color);
    size_t nc = 0;
    graph_color_dsatur(g, color, &nc);
    check_size("colours needed (n)", n, nc);
    int diam = 0, rad = 0;
    graph_diameter_radius(g, &diam, &rad);
    check_prop("diameter / radius", "1 / 1", diam == 1 && rad == 1, "%d / %d", diam, rad);
    double avg = 0;
    graph_clustering_coefficient(g, NULL, &avg);
    check_num("clustering coefficient", 1.0, avg);
    free(clique); free(color);
    graph_destroy(g);
}

static void family_tree(GraphRep rep, size_t n)
{
    char title[96];
    snprintf(title, sizeof title, "random tree of %zu vertices (%s)", n, rep_name(rep));
    subsection(title);
    Graph *g = graph_gen_random_tree(rep, n, 4);
    check_size("edges (n-1)", n - 1, graph_edge_count(g));
    check_bool("is a tree", true, graph_is_tree(g));
    check_bool("has no cycle", false, graph_has_cycle(g));
    GraphEdge *br = xmalloc(n * sizeof *br);
    size_t nbr = 0;
    graph_bridges(g, br, &nbr);
    check_size("bridges (every edge)", n - 1, nbr);
    free(br);
    graph_destroy(g);
}

static void family_ring_of_cliques(GraphRep rep, size_t k)
{
    char title[96];
    snprintf(title, sizeof title, "ring of %zu groups of 5 friends = %zu vertices (%s)", k, 5 * k, rep_name(rep));
    subsection(title);
    const size_t n = 5 * k;
    Graph *g = graph_create_sized(n, false);
    graph_set_adaptive(g, false);
    for (size_t c = 0; c < k; c++) {
        for (size_t a = 0; a < 5; a++) for (size_t b = a + 1; b < 5; b++) graph_add_edge(g, (int)(c * 5 + a), (int)(c * 5 + b), 1.0);
        graph_add_edge(g, (int)(c * 5), (int)(((c + 1) % k) * 5 + 4), 1.0);
    }
    graph_convert(g, rep);
    check_size("edges (10 per group + 1 ring link per group)", 11 * k, graph_edge_count(g));
    long long tri = 0;
    graph_triangle_count(g, &tri);
    check_int("triangles (10 per group)", (long)(10 * k), (long)tri);
    int *comp = xmalloc(n * sizeof *comp);
    check_size("connected components", 1, graph_connected_components(g, comp));
    GraphGroups cl = {0};
    graph_maximal_cliques(g, &cl);
    check_size("maximal cliques (k groups + k link edges)", 2 * k, cl.count);
    graph_groups_free(&cl);
    int *clique = xmalloc(n * sizeof *clique);
    size_t cs = 0;
    graph_max_clique(g, clique, &cs);
    check_size("largest clique", 5, cs);
    GraphEdge *br = xmalloc(n * sizeof *br);
    size_t nbr = 9;
    graph_bridges(g, br, &nbr);
    check_size("bridges (the ring has no weak link)", 0, nbr);
    free(comp); free(clique); free(br);
    graph_destroy(g);
}

static void family_dag(GraphRep rep, size_t n)
{
    char title[96];
    snprintf(title, sizeof title, "directed chain 0 -> 1 -> ... of %zu vertices (%s)", n, rep_name(rep));
    subsection(title);
    Graph *g = graph_gen_path(rep, n, true);
    check_bool("is a DAG", true, graph_is_dag(g));
    check_bool("strongly connected", false, graph_is_strongly_connected(g));
    int *order = xmalloc(n * sizeof *order);
    bool ok = graph_topological_sort(g, order);
    bool identity = ok;
    for (size_t i = 0; i < n && identity; i++) if (order[i] != (int)i) identity = false;
    check_bool("topological order is 0,1,2,...", true, identity);
    double *dist = xmalloc(n * sizeof *dist);
    int *par = xmalloc(n * sizeof *par);
    int end = -1;
    graph_dag_longest_path(g, dist, par, &end);
    check_num("longest path (n-1)", (double)(n - 1), end >= 0 ? dist[end] : -1);
    free(order); free(dist); free(par);
    graph_destroy(g);
}

static void part1_families(void)
{
    section("1. Graph families, by size (answers from closed formulas)");
    const GraphRep reps[2] = {GRAPH_LIST, GRAPH_MATRIX};
    for (int r = 0; r < 2; r++) {
        family_path(reps[r], 10, true);
        family_cycle(reps[r], 10, true);
        family_cycle(reps[r], 11, true);
        family_star(reps[r], 10, true);
        family_grid(reps[r], 5, true);
        family_complete(reps[r], 12);
        family_tree(reps[r], 50);
        family_ring_of_cliques(reps[r], 6);
        family_dag(reps[r], 10);
    }
    /* Larger sizes, list only (the matrix is for small and dense graphs; see section 6). */
    family_path(GRAPH_LIST, 1000, true);
    family_cycle(GRAPH_LIST, 1000, true);
    family_star(GRAPH_LIST, 1000, true);
    family_grid(GRAPH_LIST, 30, true);
    family_complete(GRAPH_LIST, 100);
    family_tree(GRAPH_LIST, 1000);
    family_ring_of_cliques(GRAPH_LIST, 200);
    family_path(GRAPH_LIST, 100000, false);
    family_cycle(GRAPH_LIST, 100000, false);
    family_star(GRAPH_LIST, 100000, false);
    family_grid(GRAPH_LIST, 316, false);
    family_tree(GRAPH_LIST, 100000);
    family_ring_of_cliques(GRAPH_LIST, 20000);
    family_dag(GRAPH_LIST, 100000);
    family_path(GRAPH_LIST, 1000000, false);
    family_cycle(GRAPH_LIST, 1000000, false);
    family_star(GRAPH_LIST, 1000000, false);
    family_grid(GRAPH_LIST, 1000, false);
    family_tree(GRAPH_LIST, 1000000);
    family_dag(GRAPH_LIST, 1000000);
}

/* ------------------------------------------------------------------ */
/*  2. Textbook examples                                                */
/* ------------------------------------------------------------------ */

static void part2_textbook(void)
{
    section("2. Textbook examples (published answers, re-checked by a separate script)");
    const GraphRep reps[2] = {GRAPH_LIST, GRAPH_MATRIX};
    for (int r = 0; r < 2; r++) {
        char title[96];
        snprintf(title, sizeof title, "(%s)", rep_name(reps[r]));
        subsection(title);

        /* Dijkstra: s=0 t=1 y=2 x=3 z=4; shortest distances from s are t=8 y=5 x=9 z=7 */
        const double dj[][3] = {{0, 1, 10}, {0, 2, 5}, {1, 3, 1}, {1, 2, 2}, {2, 1, 3}, {2, 3, 9}, {2, 4, 2}, {3, 4, 4}, {4, 0, 7}, {4, 3, 6}};
        Graph *g = build(reps[r], 5, true, dj, 10);
        double d[5];
        int p[5];
        graph_dijkstra(g, 0, d, p);
        check_prop("Dijkstra from s: distances to t, y, x, z", "8, 5, 9, 7", d[1] == 8 && d[2] == 5 && d[3] == 9 && d[4] == 7,
                   "%g, %g, %g, %g", d[1], d[2], d[3], d[4]);
        GraphPath path = {0};
        graph_a_star(g, 0, 3, NULL, NULL, &path);
        check_num("A* s to x: cost", 9.0, path.cost);
        graph_path_free(&path);
        double fw[25];
        int nx[25];
        graph_floyd_warshall(g, fw, nx);
        check_prop("Floyd-Warshall: s row equals Dijkstra", "0, 8, 5, 9, 7", fw[0] == 0 && fw[1] == 8 && fw[2] == 5 && fw[3] == 9 && fw[4] == 7,
                   "%g, %g, %g, %g, %g", fw[0], fw[1], fw[2], fw[3], fw[4]);
        graph_destroy(g);

        /* Bellman-Ford with a negative edge */
        const double bf[][3] = {{0, 1, 4}, {0, 2, 5}, {2, 1, -3}, {1, 3, 2}};
        g = build(reps[r], 4, true, bf, 4);
        double bd[4];
        int bp[4];
        check_status("Bellman-Ford with a negative edge", GRAPH_OK, graph_bellman_ford(g, 0, bd, bp));
        check_prop("  distances (0->2->1 beats 0->1)", "0, 2, 5, 4", bd[0] == 0 && bd[1] == 2 && bd[2] == 5 && bd[3] == 4, "%g, %g, %g, %g", bd[0], bd[1], bd[2], bd[3]);
        double jd[16];
        int jn[16];
        graph_johnson(g, jd, jn);
        check_prop("  Johnson agrees from vertex 0", "0, 2, 5, 4", jd[0] == 0 && jd[1] == 2 && jd[2] == 5 && jd[3] == 4, "%g, %g, %g, %g", jd[0], jd[1], jd[2], jd[3]);
        graph_destroy(g);
        const double nc[][3] = {{0, 1, 1}, {1, 2, -3}, {2, 0, 1}};
        g = build(reps[r], 3, true, nc, 3);
        check_status("Bellman-Ford on a negative cycle", GRAPH_ERR_NEGATIVE_CYCLE, graph_bellman_ford(g, 0, bd, bp));
        graph_destroy(g);

        /* Maximum flow: the network of CLRS section 26.1 has maximum flow 23 */
        const double fl[][3] = {{0, 1, 16}, {0, 2, 13}, {1, 3, 12}, {2, 1, 4}, {2, 4, 14}, {3, 2, 9}, {3, 5, 20}, {4, 3, 7}, {4, 5, 4}};
        g = build(reps[r], 6, true, fl, 9);
        double val = 0;
        graph_max_flow(g, 0, 5, GRAPH_FLOW_DINIC, &val, NULL);
        check_num("maximum flow, Dinic", 23.0, val);
        graph_max_flow(g, 0, 5, GRAPH_FLOW_EDMONDS_KARP, &val, NULL);
        check_num("maximum flow, Edmonds-Karp", 23.0, val);
        GraphEdge cut[16];
        size_t ncut = 0;
        graph_min_cut(g, 0, 5, &val, cut, &ncut);
        double capsum = 0;
        for (size_t i = 0; i < ncut; i++) { double w = 0; graph_has_edge(g, cut[i].u, cut[i].v, &w); capsum += w; }
        check_prop("minimum cut = maximum flow, and its edges weigh that much", "23 / 23", val == 23 && capsum == 23, "%g / %g", val, capsum);
        graph_destroy(g);

        /* Kruskal's example from the literature: weight 39 */
        const double mk[][3] = {{0, 1, 7}, {0, 3, 5}, {1, 2, 8}, {1, 3, 9}, {1, 4, 7}, {2, 4, 5}, {3, 4, 15}, {3, 5, 6}, {4, 5, 8}, {4, 6, 9}, {5, 6, 11}};
        g = build(reps[r], 7, false, mk, 11);
        GraphEdge me[8];
        size_t mc = 0;
        check_num("minimum spanning tree weight, Kruskal", 39.0, graph_mst_kruskal(g, me, &mc));
        check_size("  edges in the tree (n-1)", 6, mc);
        check_num("minimum spanning tree weight, Prim", 39.0, graph_mst_prim(g, me, &mc));
        graph_destroy(g);

        /* Strongly connected components: {0,1,2} and {3,4,5} */
        const double sc[][3] = {{0, 1, 1}, {1, 2, 1}, {2, 0, 1}, {1, 3, 1}, {3, 4, 1}, {4, 5, 1}, {5, 3, 1}};
        g = build(reps[r], 6, true, sc, 7);
        int comp[6];
        size_t ncomp = 0;
        graph_strongly_connected_components(g, comp, &ncomp);
        check_prop("strongly connected components", "2: {0,1,2} {3,4,5}", ncomp == 2 && comp[0] == comp[1] && comp[1] == comp[2] && comp[3] == comp[4] && comp[4] == comp[5] && comp[0] != comp[3],
                   "%zu groups", ncomp);
        graph_destroy(g);

        /* Two triangles sharing vertex 2 (a bowtie): one cut vertex, two blocks, no bridges */
        const double bt[][3] = {{0, 1, 1}, {1, 2, 1}, {0, 2, 1}, {2, 3, 1}, {3, 4, 1}, {2, 4, 1}};
        g = build(reps[r], 5, false, bt, 6);
        bool cutv[5];
        size_t nc2 = 0, nbr = 0;
        graph_articulation_points(g, cutv, &nc2);
        GraphEdge bridges[8];
        graph_bridges(g, bridges, &nbr);
        GraphGroups blocks = {0};
        graph_biconnected_components(g, &blocks);
        check_prop("bowtie: cut vertices / bridges / blocks", "1 (vertex 2) / 0 / 2", nc2 == 1 && cutv[2] && nbr == 0 && blocks.count == 2,
                   "%zu / %zu / %zu", nc2, nbr, blocks.count);
        graph_groups_free(&blocks);
        graph_destroy(g);

        /* Petersen: 10 vertices, 3 friends each; famous for having a Hamiltonian path but no Hamiltonian cycle */
        Graph *pg = graph_create_sized(10, false);
        graph_set_adaptive(pg, false);
        for (int i = 0; i < 5; i++) { graph_add_edge(pg, i, (i + 1) % 5, 1); graph_add_edge(pg, i, i + 5, 1); graph_add_edge(pg, 5 + i, 5 + (i + 2) % 5, 1); }
        graph_convert(pg, reps[r]);
        check_size("Petersen graph: edges", 15, graph_edge_count(pg));
        int ham[10];
        check_status("  Hamiltonian path exists", GRAPH_OK, graph_hamiltonian(pg, false, ham));
        check_status("  Hamiltonian cycle does not", GRAPH_ERR_NOT_FOUND, graph_hamiltonian(pg, true, ham));
        bool bip = true;
        graph_is_bipartite(pg, &bip, NULL);
        check_bool("  bipartite (it has 5-cycles)", false, bip);
        long long tri = -1;
        graph_triangle_count(pg, &tri);
        check_int("  triangles", 0, (long)tri);
        int diam = 0, rad = 0;
        graph_diameter_radius(pg, &diam, &rad);
        check_prop("  diameter / radius", "2 / 2", diam == 2 && rad == 2, "%d / %d", diam, rad);
        graph_destroy(pg);

        /* Hungarian assignment: minimum total 5 */
        const double cost[9] = {4, 1, 3, 2, 0, 5, 3, 2, 2};
        int assign[3];
        double total = 0;
        graph_hungarian(cost, 3, 3, assign, &total);
        check_num("assignment problem, minimum total cost", 5.0, total);

        /* Isomorphism: a 5-cycle against a relabelled 5-cycle, and against a path */
        const double c5[][3] = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 4, 1}, {4, 0, 1}};
        const double c5r[][3] = {{3, 0, 1}, {0, 4, 1}, {4, 1, 1}, {1, 2, 1}, {2, 3, 1}};
        const double p5[][3] = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 4, 1}};
        Graph *a = build(reps[r], 5, false, c5, 5), *b = build(reps[r], 5, false, c5r, 5), *c = build(reps[r], 5, false, p5, 4);
        bool iso = false;
        graph_is_isomorphic(a, b, &iso, NULL);
        check_bool("isomorphism: C5 vs relabelled C5", true, iso);
        graph_is_isomorphic(a, c, &iso, NULL);
        check_bool("isomorphism: C5 vs a path", false, iso);
        graph_destroy(a); graph_destroy(b); graph_destroy(c);

        /* Transitive closure of a 4-vertex chain: 6 arcs; reduction brings back the 3 of the chain */
        Graph *chain = graph_gen_path(reps[r], 4, true), *closure = NULL, *reduced = NULL;
        graph_transitive_closure(chain, &closure);
        graph_transitive_reduction(closure, &reduced);
        check_prop("transitive closure / reduction of a 4-chain", "6 / 3 arcs", graph_edge_count(closure) == 6 && graph_edge_count(reduced) == 3,
                   "%zu / %zu arcs", graph_edge_count(closure), graph_edge_count(reduced));
        graph_destroy(chain); graph_destroy(closure); graph_destroy(reduced);
    }
}

/* ------------------------------------------------------------------ */
/*  3. Advanced algorithms                                              */
/* ------------------------------------------------------------------ */

static void part3_advanced(void)
{
    section("3. Advanced algorithms (graph_advanced.h)");
    const GraphRep reps[2] = {GRAPH_LIST, GRAPH_MATRIX};
    for (int r = 0; r < 2; r++) {
        char title[96];
        snprintf(title, sizeof title, "(%s)", rep_name(reps[r]));
        subsection(title);

        const double dag[][3] = {{0, 1, 2}, {1, 2, 3}, {0, 2, 1}, {2, 3, 4}};
        Graph *g = build(reps[r], 4, true, dag, 4);
        double dist[4];
        int par[4];
        graph_adv_longest_path(g, 0, dist, par);
        check_prop("longest path from 0 in a DAG", "0, 2, 5, 9", dist[0] == 0 && dist[1] == 2 && dist[2] == 5 && dist[3] == 9, "%g, %g, %g, %g", dist[0], dist[1], dist[2], dist[3]);
        graph_destroy(g);

        const double tr[][3] = {{0, 1, 1}, {0, 2, 1}, {1, 3, 1}, {1, 4, 1}};
        g = build(reps[r], 5, false, tr, 4);
        int out[5];
        size_t cnt = 0;
        graph_adv_tree_traversal(g, 0, GRAPH_ADV_TREE_PREORDER, out, &cnt);
        check_prop("tree traversal, pre-order", "0 1 3 4 2", out[0] == 0 && out[1] == 1 && out[2] == 3 && out[3] == 4 && out[4] == 2, "%d %d %d %d %d", out[0], out[1], out[2], out[3], out[4]);
        graph_adv_tree_traversal(g, 0, GRAPH_ADV_TREE_POSTORDER, out, &cnt);
        check_prop("tree traversal, post-order", "3 4 1 2 0", out[0] == 3 && out[1] == 4 && out[2] == 1 && out[3] == 2 && out[4] == 0, "%d %d %d %d %d", out[0], out[1], out[2], out[3], out[4]);
        graph_adv_tree_traversal(g, 0, GRAPH_ADV_TREE_LEVELORDER, out, &cnt);
        check_prop("tree traversal, level-order", "0 1 2 3 4", out[0] == 0 && out[1] == 1 && out[2] == 2 && out[3] == 3 && out[4] == 4, "%d %d %d %d %d", out[0], out[1], out[2], out[3], out[4]);
        graph_destroy(g);

        const double wp[][3] = {{0, 1, 5}, {1, 2, 3}, {0, 2, 2}, {2, 3, 4}};
        g = build(reps[r], 4, false, wp, 4);
        double width[4];
        graph_adv_widest_path(g, 0, width, NULL);
        check_prop("widest path from 0 (the best bottleneck)", "5, 3, 3", width[1] == 5 && width[2] == 3 && width[3] == 3, "%g, %g, %g", width[1], width[2], width[3]);
        graph_destroy(g);

        const double zc[][3] = {{0, 1, 1}, {1, 2, 2}, {2, 0, -3}};
        g = build(reps[r], 3, true, zc, 3);
        GraphAdvCycleKind kind;
        int cyc[3];
        size_t clen = 0;
        graph_adv_zero_weight_cycle(g, 1e-9, &kind, cyc, &clen);
        check_prop("zero-weight cycle (1 + 2 - 3)", "found, 3 vertices", kind == GRAPH_ADV_CYCLE_ZERO && clen == 3, "%s, %zu vertices", kind == GRAPH_ADV_CYCLE_ZERO ? "found" : "not found", clen);
        graph_destroy(g);

        const double k4[][3] = {{0, 1, 1}, {0, 2, 1}, {0, 3, 1}, {1, 2, 1}, {1, 3, 1}, {2, 3, 1}};
        const double tri[][3] = {{0, 1, 1}, {1, 2, 1}, {0, 2, 1}};
        const double path4[][3] = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}};
        Graph *big = build(reps[r], 4, false, k4, 6), *pat = build(reps[r], 3, false, tri, 3), *pth = build(reps[r], 4, false, path4, 3);
        bool found = false;
        int map[3];
        graph_adv_subgraph_isomorphism(pat, big, true, false, map, &found);
        check_bool("a triangle occurs inside K4", true, found);
        graph_adv_subgraph_isomorphism(pat, pth, false, false, NULL, &found);
        check_bool("a triangle occurs inside a path", false, found);
        int order[4] = {0, 1, 2, 3};
        long long crossings = 0;
        graph_adv_crossing_count(big, order, 4, &crossings);
        check_int("K4 drawn on a circle: edge crossings", 1, (long)crossings);
        graph_destroy(big); graph_destroy(pat); graph_destroy(pth);

        Graph *p7 = graph_gen_path(reps[r], 7, false);
        int color[7];
        size_t ncolors = 0, rounds = 0;
        graph_adv_color_refinement(p7, color, &ncolors, &rounds);
        check_size("colour refinement of a 7-path: classes (4)", 4, ncolors);
        double harm[7];
        graph_adv_harmonic_centrality(p7, harm);
        const double h1 = 1.0 + 1.0 / 2 + 1.0 / 3 + 1.0 / 4 + 1.0 / 5 + 1.0 / 6;          /* an end: 1 + 1/2 + ... + 1/6 */
        check_num("harmonic centrality of an end of the 7-path", h1 / 6.0, harm[0]);
        graph_destroy(p7);

        Graph *gr = graph_gen_grid(reps[r], 5, 5);
        int round[25], parent[25];
        size_t messages = 0, reached = 0;
        graph_adv_flood(gr, 0, 100, round, parent, &messages, &reached);
        check_prop("flood across a 5 x 5 grid from a corner", "25 reached, 8 rounds", reached == 25 && round[24] == 8, "%zu reached, %d rounds", reached, round[24]);
        graph_destroy(gr);

        const double c5[][3] = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 4, 1}, {4, 0, 1}};
        const double c5r[][3] = {{3, 0, 2}, {0, 4, 2}, {4, 1, 2}, {1, 2, 2}, {2, 3, 2}};
        Graph *x = build(reps[r], 5, false, c5, 5), *y = build(reps[r], 5, false, c5r, 5);
        bool iso = false;
        graph_adv_graph_isomorphism(x, y, false, &iso, NULL);
        check_bool("isomorphism ignoring weights: C5 vs relabelled C5 (weights 2)", true, iso);
        graph_adv_graph_isomorphism(x, y, true, &iso, NULL);
        check_bool("isomorphism matching weights (1 vs 2)", false, iso);
        graph_destroy(x); graph_destroy(y);

        /* dynamic connectivity follows edits and representation changes */
        const double two[][3] = {{0, 1, 1}, {2, 3, 1}};
        g = build(reps[r], 4, false, two, 2);
        GraphAdvDynConn *dc = graph_adv_dynconn_create(g);
        size_t comps = 0;
        bool con = false;
        graph_adv_dynconn_components(dc, &comps);
        check_size("dynamic connectivity: groups at the start", 2, comps);
        graph_add_edge(g, 1, 2, 1.0);
        graph_adv_dynconn_connected(dc, 0, 3, &con);
        check_bool("  after adding 1-2: are 0 and 3 connected?", true, con);
        graph_remove_edge(g, 1, 2);
        graph_adv_dynconn_connected(dc, 0, 3, &con);
        check_bool("  after removing it again", false, con);
        graph_adv_dynconn_destroy(dc);
        graph_destroy(g);
    }
}

/* ------------------------------------------------------------------ */
/*  4. Social queries                                                   */
/* ------------------------------------------------------------------ */

static void part4_social(void)
{
    section("4. Social queries (graph_social.h): a small network worked out by hand");
    printf("  friendships: 0-1, 0-2, 1-3, 2-3, 1-4, 6-7;  user 5 has no friends\n");
    const double e[][3] = {{0, 1, 1}, {0, 2, 1}, {1, 3, 1}, {2, 3, 1}, {1, 4, 1}, {6, 7, 1}};
    const GraphRep reps[2] = {GRAPH_LIST, GRAPH_MATRIX};
    for (int r = 0; r < 2; r++) {
        char title[96];
        snprintf(title, sizeof title, "(%s)", rep_name(reps[r]));
        subsection(title);
        Graph *g = build(reps[r], 8, false, e, 6);
        int users[8], hops[8];
        size_t count = 0;
        bool trunc = false;
        graph_social_within_hops(g, 0, 2, 8, users, hops, &count, &trunc);
        int ones = 0, twos = 0;
        for (size_t i = 0; i < count; i++) { ones += hops[i] == 1; twos += hops[i] == 2; }
        check_prop("people within 2 hops of user 0", "4: two at 1 hop, two at 2", count == 4 && ones == 2 && twos == 2, "%zu: %d at 1 hop, %d at 2", count, ones, twos);
        graph_social_within_hops(g, 0, 1, 8, users, hops, &count, &trunc);
        check_size("people within 1 hop of user 0 (just friends)", 2, count);
        graph_social_within_hops(g, 5, 3, 8, users, hops, &count, &trunc);
        check_size("people within 3 hops of user 5 (knows nobody)", 0, count);
        graph_social_within_hops(g, 0, 5, 2, users, hops, &count, &trunc);
        check_prop("result limit of 2 keeps the nearest and says it was cut", "2, truncated", count == 2 && trunc, "%zu, %s", count, trunc ? "truncated" : "complete");
        GraphSuggestion sug[4];
        graph_social_suggest_friends(g, 0, 4, sug, &count);
        check_prop("people you may know (user 0)", "3 (2 mutual), 4 (1 mutual)", count == 2 && sug[0].user == 3 && sug[0].mutual == 2 && sug[1].user == 4 && sug[1].mutual == 1,
                   "%d (%zu), %d (%zu)", count > 0 ? sug[0].user : -1, count > 0 ? sug[0].mutual : 0, count > 1 ? sug[1].user : -1, count > 1 ? sug[1].mutual : 0);
        GraphPath p = {0};
        graph_social_path(g, 0, 4, SIZE_MAX, &p);
        check_prop("shortest chain of friends 0 -> 4", "0 - 1 - 4 (2 hops)", p.length == 3 && p.vertices[1] == 1 && p.cost == 2, "%zu people, %g hops", p.length, p.cost);
        graph_path_free(&p);
        check_status("chain of friends 0 -> 6 (different groups)", GRAPH_ERR_NOT_FOUND, graph_social_path(g, 0, 6, SIZE_MAX, &p));
        check_status("chain of friends 0 -> 4 within 1 hop only", GRAPH_ERR_NOT_FOUND, graph_social_path(g, 0, 4, 1, &p));
        graph_destroy(g);
    }
}

/* ------------------------------------------------------------------ */
/*  5. Files too big for memory: streaming and sharding                 */
/* ------------------------------------------------------------------ */

static FILE *ring_file(size_t k)          /* a ring of k groups of 5 friends, as an edge-list file */
{
    FILE *f = tmpfile();
    for (size_t c = 0; c < k; c++) {
        for (size_t a = 0; a < 5; a++) for (size_t b = a + 1; b < 5; b++) fprintf(f, "%zu %zu\n", c * 5 + a, c * 5 + b);
        fprintf(f, "%zu %zu\n", c * 5, ((c + 1) % k) * 5 + 4);
    }
    fflush(f);
    rewind(f);
    return f;
}

static size_t files_in(const char *dir)
{
    char cmd[300];
    snprintf(cmd, sizeof cmd, "ls -A %s 2>/dev/null | wc -l", dir);
    FILE *p = popen(cmd, "r");
    size_t n = 0;
    if (p != NULL) { if (fscanf(p, "%zu", &n) != 1) n = 0; pclose(p); }
    return n;
}

static void part5_files(void)
{
    section("5. Edge-list files: streaming and the sharding engine");
    printf("  file A: 200 groups of 5 friends in a ring (1000 users, 2200 friendships)\n");
    printf("          every user has 4 friends, except the 400 at the ring links who have 5\n");
    FILE *f = ring_file(200);

    subsection("streaming: one pass, nothing but a few integers per user held in memory");
    GraphStreamInfo info;
    graph_stream_scan(f, &info);
    check_size("users found", 1000, info.vertices);
    check_size("friendships found", 2200, info.edges);
    rewind(f);
    size_t *deg = xmalloc(1000 * sizeof *deg);
    graph_stream_degrees(f, 1000, false, deg, NULL);
    size_t d4 = 0, d5 = 0;
    for (size_t v = 0; v < 1000; v++) { d4 += deg[v] == 4; d5 += deg[v] == 5; }
    check_prop("users with 4 / 5 friends", "600 / 400", d4 == 600 && d5 == 400, "%zu / %zu", d4, d5);
    rewind(f);
    int *comp = xmalloc(1000 * sizeof *comp);
    size_t groups = 0;
    graph_stream_components(f, 1000, comp, &groups);
    check_size("connected groups", 1, groups);

    subsection("sharding engine forced on, 512 KB budget: friendships live on disk");
    char dir[] = "/tmp/verify_shards_XXXXXX";
    if (mkdtemp(dir) == NULL) { fprintf(stderr, "mkdtemp failed\n"); exit(2); }
    ShardConfig cfg = shard_default_config();
    cfg.memory_budget = 512 * 1024;
    cfg.dir = dir;
    cfg.force = true;
    Graph *mem = NULL;
    ShardedGraph *sg = NULL;
    check_status("open with sharding forced", GRAPH_OK, graph_open_edge_list(f, &cfg, &mem, &sg));
    check_bool("  came back as shards (not an in-memory graph)", true, sg != NULL && mem == NULL);
    check_bool("  shard files exist in the directory", true, files_in(dir) > 0);
    check_size("users", 1000, shard_vertex_count(sg));
    check_size("friendships", 2200, shard_edge_count(sg));
    size_t sd4 = 0, sd5 = 0;
    for (int v = 0; v < 1000; v++) { sd4 += shard_degree(sg, v) == 4; sd5 += shard_degree(sg, v) == 5; }
    check_prop("users with 4 / 5 friends", "600 / 400", sd4 == 600 && sd5 == 400, "%zu / %zu", sd4, sd5);
    size_t sgroups = 0;
    shard_components(sg, comp, &sgroups);
    check_size("connected groups", 1, sgroups);
    double *rank = xmalloc(1000 * sizeof *rank);
    size_t rounds = 0;
    GraphStatus st = shard_pagerank(sg, 0.85, 1e-12, 500, rank, &rounds);
    double sum = 0;
    for (int v = 0; v < 1000; v++) sum += rank[v];
    check_status("PageRank converges", GRAPH_OK, st);
    check_num("  ranks add up to", 1.0, sum);
    /* the same graph held in memory: the answers must be identical */
    rewind(f);
    Graph *ref = NULL;
    graph_read_edge_list(f, false, GRAPH_LIST, &ref);
    double *rank2 = xmalloc(1000 * sizeof *rank2);
    size_t r2 = 0;
    graph_pagerank(ref, 0.85, 1e-12, 500, rank2, &r2);
    double worst = 0;
    for (int v = 0; v < 1000; v++) worst = fmax(worst, fabs(rank[v] - rank2[v]));
    check_prop("  same ranks as the in-memory graph (largest difference)", "< 1e-9", worst < 1e-9, "%.1e", worst);
    int *level = xmalloc(1000 * sizeof *level);
    size_t reached = 0;
    shard_bfs(sg, 0, 2, level, &reached);
    double *dist = xmalloc(1000 * sizeof *dist);
    int *par = xmalloc(1000 * sizeof *par);
    graph_dijkstra(ref, 0, dist, par);
    size_t want = 0;
    for (int v = 0; v < 1000; v++) want += dist[v] <= 2.0;
    check_size("people within 2 hops of user 0 (same as in memory)", want, reached);
    shard_destroy(sg);
    check_size("shard files left behind after shard_destroy", 0, files_in(dir));
    rmdir(dir);
    graph_destroy(ref);
    free(deg); free(comp); free(rank); free(rank2); free(level); free(dist); free(par);
    fclose(f);

    printf("  file B: '0 1', '1 2', '3 4' and users 5, 6, 7 who know nobody ('# vertices 8')\n");
    f = tmpfile();
    fputs("# vertices 8\n0 1\n1 2\n3 4\n", f);
    rewind(f);
    graph_stream_scan(f, &info);
    check_size("users found (silent ones are counted)", 8, info.vertices);
    rewind(f);
    int comp8[8];
    graph_stream_components(f, 8, comp8, &groups);
    check_size("connected groups ({0,1,2} {3,4} {5} {6} {7})", 5, groups);
    fclose(f);

    subsection("mode switch: the same call picks in-memory or shards from the budget");
    f = ring_file(200);
    cfg = shard_default_config();
    cfg.memory_budget = (size_t)256 << 20;
    Graph *g = NULL;
    ShardedGraph *s2 = NULL;
    graph_open_edge_list(f, &cfg, &g, &s2);
    check_bool("256 MiB budget: loaded as an ordinary graph", true, g != NULL && s2 == NULL);
    graph_destroy(g);
    cfg.memory_budget = 100 * 1024;
    g = NULL;
    graph_open_edge_list(f, &cfg, &g, &s2);
    check_bool("100 KB budget: switched to shards by itself", true, g == NULL && s2 != NULL);
    shard_destroy(s2);
    cfg.memory_budget = 4 * 1024;
    check_status("4 KB budget: even the per-user state does not fit", GRAPH_ERR_TOO_LARGE, graph_open_edge_list(f, &cfg, &g, &s2));
    fclose(f);
}

/* ------------------------------------------------------------------ */
/*  6. Memory limits                                                    */
/* ------------------------------------------------------------------ */

static void part6_memory(void)
{
    section("6. Memory budget and the matrix size limit");
    check_size("default memory budget (MiB)", 4096, graph_memory_budget() >> 20);
    check_size("largest weighted matrix (sqrt(8 * 0.5 * 4 GiB))", 131072, graph_matrix_max_vertices(true));
    check_size("largest unweighted matrix (sqrt(8 * 0.9 * 4 GiB))", 175851, graph_matrix_max_vertices(false));

    subsection("a graph with too many users for the matrix stays a list");
    Graph *g = graph_create_sized(200000, false);
    check_bool("convert 200 000 users to a matrix succeeds", false, graph_convert(g, GRAPH_MATRIX));
    check_status("  the reason", GRAPH_ERR_TOO_LARGE, graph_last_error());
    check_bool("  still a list, still usable", true, graph_representation(g) == GRAPH_LIST && graph_add_edge(g, 1, 2, 1.0));
    graph_destroy(g);
    GraphConfig cfg = graph_default_config(false);
    cfg.start = GRAPH_MATRIX;
    check_bool("start a graph of 1 000 000 users as a matrix", false, graph_create_with(&cfg, 1000000) != NULL);

    subsection("unweighted mode: no weights stored, every edge is 1");
    cfg = graph_default_config(false);
    cfg.unweighted = true;
    g = graph_create_with(&cfg, 100);
    graph_add_edge(g, 1, 2, 1.0);
    check_bool("add_edge with weight 2.5", false, graph_add_edge(g, 3, 4, 2.5));
    check_status("  the reason", GRAPH_ERR_INVALID_WEIGHT, graph_last_error());
    double w = 0;
    graph_has_edge(g, 1, 2, &w);
    check_num("weight reported for the edge that exists", 1.0, w);
    graph_destroy(g);

    subsection("a lowered budget refuses cleanly instead of crashing");
    const size_t saved = graph_memory_budget();
    graph_set_memory_budget((size_t)1 << 20);
    Graph *huge = graph_create_sized(5000000, false);
    check_bool("create 5 000 000 users under a 1 MB budget", false, huge != NULL);
    check_status("  the reason", GRAPH_ERR_TOO_LARGE, graph_last_error());
    Graph *small = graph_create_sized(20000, false);
    size_t added = 0;
    bool refused = false;
    for (int i = 1; i < 20000 && !refused; i++)
        for (int k = 1; k <= 8 && !refused; k++) {
            const int v = (i * k * 31 + 7) % 20000;
            if (v == i || graph_has_edge(small, i, v, NULL)) continue;
            if (graph_add_edge(small, i, v, 1.0)) added++; else refused = graph_last_error() == GRAPH_ERR_TOO_LARGE;
        }
    check_prop("adding friendships until the budget is reached", "refused (too large)", refused, "refused after %zu edges", added);
    check_size("  the graph kept exactly the edges that were accepted", added, graph_edge_count(small));
    graph_destroy(small);
    graph_set_memory_budget(saved);

    subsection("shrink-to-fit returns slack and keeps the contents");
    Graph *sg = graph_gen_barabasi_albert(GRAPH_LIST, 5000, 4, 5);
    const size_t edges = graph_edge_count(sg);
    double w1 = 0, w2 = 0;
    graph_has_edge(sg, 1, 0, &w1);
    check_bool("shrink_to_fit succeeds", true, graph_shrink_to_fit(sg));
    graph_has_edge(sg, 1, 0, &w2);
    check_prop("  edge count and a sample edge unchanged", "same", graph_edge_count(sg) == edges && w1 == w2, "%zu edges", graph_edge_count(sg));
    graph_destroy(sg);
}

/* ------------------------------------------------------------------ */
/*  7. List and matrix must agree                                       */
/* ------------------------------------------------------------------ */

static void part7_agree(void)
{
    section("7. The two representations agree (a random weighted graph, 150 users)");
    printf("  expected = the list's answer, actual = the matrix's answer\n");
    unsigned long long x = 88172645463325252ULL;
    Graph *a = graph_create_sized(150, false);
    graph_set_adaptive(a, false);
    for (int i = 0; i < 700; i++) {
        x ^= x << 13; x ^= x >> 7; x ^= x << 17;
        const int u = (int)(x % 150), v = (int)((x >> 20) % 150);
        if (u != v) graph_add_edge(a, u, v, (double)(1 + (x >> 40) % 9));
    }
    Graph *b = graph_copy(a);
    graph_convert(b, GRAPH_MATRIX);
    check_bool("setup: one graph is a list, the copy a matrix", true, graph_representation(a) == GRAPH_LIST && graph_representation(b) == GRAPH_MATRIX);
    check_size("edges", graph_edge_count(a), graph_edge_count(b));

    int ca[150], cb[150];
    check_size("connected components", graph_connected_components(a, ca), graph_connected_components(b, cb));
    long long ta = 0, tb = 0;
    graph_triangle_count(a, &ta); graph_triangle_count(b, &tb);
    check_int("triangles", (long)ta, (long)tb);
    int da = 0, db = 0, ra = 0, rb = 0;
    graph_diameter_radius(a, &da, &ra); graph_diameter_radius(b, &db, &rb);
    check_prop("diameter / radius", "same", da == db && ra == rb, "%d / %d", db, rb);
    double dista[150], distb[150], sa = 0, sb = 0;
    int pa[150], pb[150];
    graph_dijkstra(a, 0, dista, pa); graph_dijkstra(b, 0, distb, pb);
    for (int v = 0; v < 150; v++) { if (dista[v] < 1e300) sa += dista[v]; if (distb[v] < 1e300) sb += distb[v]; }
    check_num("sum of shortest distances from user 0", sa, sb);
    GraphEdge ea[150], eb[150];
    size_t na = 0, nb = 0;
    check_num("minimum spanning tree weight", graph_mst_prim(a, ea, &na), graph_mst_prim(b, eb, &nb));
    double bca[150], bcb[150], ba = 0, bb = 0;
    graph_betweenness_centrality(a, bca); graph_betweenness_centrality(b, bcb);
    for (int v = 0; v < 150; v++) { ba += bca[v]; bb += bcb[v]; }
    check_num("sum of betweenness centrality", ba, bb);
    double pra[150], prb[150], worst = 0;
    size_t it = 0;
    graph_pagerank(a, 0.85, 1e-12, 500, pra, &it); graph_pagerank(b, 0.85, 1e-12, 500, prb, &it);
    for (int v = 0; v < 150; v++) worst = fmax(worst, fabs(pra[v] - prb[v]));
    check_prop("PageRank, largest difference", "< 1e-9", worst < 1e-9, "%.1e", worst);
    int ka[150], kb[150];
    graph_core_numbers(a, ka); graph_core_numbers(b, kb);
    check_prop("core numbers of every user", "identical", memcmp(ka, kb, sizeof ka) == 0, "%s", memcmp(ka, kb, sizeof ka) == 0 ? "identical" : "differ");
    size_t sza = 0, szb = 0;
    int clique[150];
    graph_max_clique(a, clique, &sza); graph_max_clique(b, clique, &szb);
    check_size("largest clique", sza, szb);
    GraphGroups mca = {0}, mcb = {0};
    graph_maximal_cliques(a, &mca); graph_maximal_cliques(b, &mcb);
    check_size("maximal cliques", mca.count, mcb.count);
    graph_groups_free(&mca); graph_groups_free(&mcb);
    int col_a[150], col_b[150];
    size_t nca = 0, ncb = 0;
    graph_color_dsatur(a, col_a, &nca); graph_color_dsatur(b, col_b, &ncb);
    check_size("colours used by DSatur", nca, ncb);
    double ja = 0, jb = 0;
    graph_jaccard(a, 1, 2, &ja); graph_jaccard(b, 1, 2, &jb);
    check_num("Jaccard similarity of users 1 and 2", ja, jb);
    graph_destroy(a); graph_destroy(b);
}

int main(int argc, char **argv)
{
    g_color = isatty(STDOUT_FILENO) != 0;
    const bool quick = argc > 1 && strcmp(argv[1], "quick") == 0;   /* skip the million-vertex graphs */
    printf("%sVerification report%s: what the answer should be, what the library gave, and whether they match\n",
           BOLD, RESET);
    if (!quick) {
        part1_families();
    } else {
        section("1. Graph families, by size (answers from closed formulas)");
        family_path(GRAPH_LIST, 10, true);
        family_grid(GRAPH_LIST, 5, true);
        family_complete(GRAPH_LIST, 12);
        family_ring_of_cliques(GRAPH_LIST, 6);
    }
    part2_textbook();
    part3_advanced();
    part4_social();
    part5_files();
    part6_memory();
    part7_agree();
    printf("  %s-> %d checks, %d failed%s\n", DIM, g_section_rows, g_section_failed, RESET);
    printf("\n%s%s%d checks, %d failed%s\n", BOLD, g_failed == 0 ? GREEN : RED, g_rows, g_failed, RESET);
    return g_failed == 0 ? 0 : 1;
}
