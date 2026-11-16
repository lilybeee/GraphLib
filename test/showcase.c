/*
 * showcase.c — run every algorithm of the library on small, fixed graphs and
 * print, for each function, the EXPECTED outcome next to the ACTUAL outcome.
 *
 * Expected values were worked out independently of the C code (by hand and by
 * a separate brute-force script), never by calling the library.  Every call is
 * made through the unified Graph API, once with the graph stored as a list and
 * once as a matrix, so the output also shows that both representations agree.
 *
 *   make showcase        # prints the table, exits non-zero if any row fails
 *
 * Where the answer is not unique (a colouring, a community split, a Hamiltonian
 * order) the expected column states the property that must hold and the check
 * verifies that property instead of one particular answer.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "graph_algorithms.h"
#include "graph_viz.h"

/* ------------------------------------------------------------------ */
/*  Result table                                                        */
/* ------------------------------------------------------------------ */

#define MAX_ROWS 400
#define TXT      220

typedef struct {
    char name[96];
    char expected[TXT];
    char actual[2][TXT];     /* [0] list, [1] matrix */
    bool ok[2];
} Row;

static Row rows[MAX_ROWS];
static int nrows, cursor, pass;      /* pass 0 = list, 1 = matrix */
static GraphRep cur_rep;

static void record(const char *name, const char *expected, const char *actual, bool ok)
{
    Row *r;
    if (pass == 0) {
        if (nrows >= MAX_ROWS) {
            fprintf(stderr, "too many rows\n");
            exit(2);
        }
        r = &rows[nrows++];
        snprintf(r->name, sizeof r->name, "%s", name);
        snprintf(r->expected, sizeof r->expected, "%s", expected);
    } else {
        r = &rows[cursor++];
        if (strcmp(r->name, name) != 0) {
            fprintf(stderr, "row mismatch between passes: %s vs %s\n", r->name, name);
            exit(2);
        }
    }
    snprintf(r->actual[pass], TXT, "%s", actual);
    r->ok[pass] = ok;
}

static void fmt_d(char *b, size_t n, double v)
{
    if (isinf(v)) { snprintf(b, n, v > 0 ? "inf" : "-inf"); }
    else if (isnan(v)) { snprintf(b, n, "nan"); }
    else if (v == floor(v) && fabs(v) < 1e9) { snprintf(b, n, "%.0f", v); }
    else { snprintf(b, n, "%.4f", v); }
}

static void fmt_ints(char *b, size_t n, const int *a, size_t len)
{
    size_t o = (size_t)snprintf(b, n, "[");
    for (size_t i = 0; i < len && o + 8 < n; ++i) { o += (size_t)snprintf(b + o, n - o, i ? " %d" : "%d", a[i]); }
    snprintf(b + o, n - o, "]");
}

static void fmt_dbls(char *b, size_t n, const double *a, size_t len)
{
    size_t o = (size_t)snprintf(b, n, "[");
    for (size_t i = 0; i < len && o + 12 < n; ++i) {
        char t[32];
        fmt_d(t, sizeof t, a[i]);
        o += (size_t)snprintf(b + o, n - o, i ? " %s" : "%s", t);
    }
    snprintf(b + o, n - o, "]");
}

static bool near(double a, double b, double tol)
{
    if (isinf(a) || isinf(b)) { return a == b; }
    return fabs(a - b) <= tol * fmax(1.0, fmax(fabs(a), fabs(b)));
}

static void chk_int(const char *name, long exp, long act)
{
    char e[TXT], a[TXT];
    snprintf(e, sizeof e, "%ld", exp);
    snprintf(a, sizeof a, "%ld", act);
    record(name, e, a, exp == act);
}

static void chk_bool(const char *name, bool exp, bool act)
{
    record(name, exp ? "true" : "false", act ? "true" : "false", exp == act);
}

static void chk_dbl(const char *name, double exp, double act, double tol)
{
    char e[TXT], a[TXT];
    fmt_d(e, sizeof e, exp);
    fmt_d(a, sizeof a, act);
    record(name, e, a, near(exp, act, tol));
}

static void chk_status(const char *name, GraphStatus exp, GraphStatus act)
{
    char e[TXT], a[TXT];
    snprintf(e, sizeof e, "%s", graph_status_string(exp));
    snprintf(a, sizeof a, "%s", graph_status_string(act));
    record(name, e, a, exp == act);
}

static void chk_ints(const char *name, const int *exp, const int *act, size_t n)
{
    char e[TXT], a[TXT];
    fmt_ints(e, sizeof e, exp, n);
    fmt_ints(a, sizeof a, act, n);
    record(name, e, a, memcmp(exp, act, n * sizeof(int)) == 0);
}

static void chk_dbls(const char *name, const double *exp, const double *act, size_t n, double tol)
{
    char e[TXT], a[TXT];
    fmt_dbls(e, sizeof e, exp, n);
    fmt_dbls(a, sizeof a, act, n);
    bool ok = true;
    for (size_t i = 0; i < n; ++i) { ok = ok && near(exp[i], act[i], tol); }
    record(name, e, a, ok);
}

static void chk_prop(const char *name, const char *expected, bool ok, const char *fmt, ...)
    __attribute__((format(printf, 4, 5)));

#include <stdarg.h>
static void chk_prop(const char *name, const char *expected, bool ok, const char *fmt, ...)
{
    char a[TXT];
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(a, sizeof a, fmt, ap);
    va_end(ap);
    record(name, expected, a, ok);
}

/* ------------------------------------------------------------------ */
/*  The graphs                                                          */
/* ------------------------------------------------------------------ */

typedef struct { int u, v; double w; } E;

static Graph *build(bool directed, size_t n, const E *edges, size_t m)
{
    GraphConfig cfg = graph_default_config(directed);
    cfg.start    = cur_rep;
    cfg.adaptive = false;                 /* keep the representation fixed for this pass */
    Graph *g = graph_create_with(&cfg, n);
    if (g == NULL) { fprintf(stderr, "out of memory\n"); exit(2); }
    for (size_t i = 0; i < m; ++i) {
        if (!graph_add_edge(g, edges[i].u, edges[i].v, edges[i].w)) {
            fprintf(stderr, "bad test edge %d-%d\n", edges[i].u, edges[i].v);
            exit(2);
        }
    }
    return g;
}

/* U: undirected, 8 vertices; two triangles {0,1,2} {3,4,5}, a pendant 6, an isolated 7. */
static const E U_EDGES[] = { {0,1,2}, {0,2,3}, {1,2,1}, {1,3,4}, {2,4,5}, {3,4,2}, {3,5,6}, {4,5,1}, {5,6,3} };
/* D: directed; cycle 0->1->2->0 then a tail 2->3->4->5 with a shortcut 3->5. */
static const E D_EDGES[] = { {0,1,1}, {1,2,2}, {2,0,3}, {2,3,1}, {3,4,2}, {4,5,1}, {3,5,5} };
/* G: a DAG with a redundant edge 0->3. */
static const E G_EDGES[] = { {0,1,3}, {0,2,2}, {1,3,4}, {2,3,1}, {3,4,2}, {0,3,1} };
/* F: a flow network (weights are capacities). */
static const E F_EDGES[] = { {0,1,10}, {0,2,5}, {1,2,15}, {1,3,5}, {2,3,10} };
/* N: negative weights, no negative cycle;  C: negative cycle. */
static const E N_EDGES[] = { {0,1,4}, {0,2,5}, {1,2,-3}, {2,3,2} };
static const E C_EDGES[] = { {0,1,1}, {1,2,-3}, {2,0,1} };
/* B: bipartite {0,1,2} | {3,4,5}. */
static const E B_EDGES[] = { {0,3,1}, {0,4,1}, {1,3,1}, {2,4,1}, {2,5,1} };
/* W: bowtie (two triangles sharing vertex 0): an Euler circuit exists. */
static const E W_EDGES[] = { {0,1,1}, {1,2,1}, {2,0,1}, {0,3,1}, {3,4,1}, {4,0,1} };
/* Z: weights 0 and 1 only. */
static const E Z_EDGES[] = { {0,1,0}, {1,2,1}, {0,2,1}, {2,3,0} };

#define COUNT(a) (sizeof(a) / sizeof((a)[0]))
#define MK(dir, n, arr) build((dir), (n), (arr), COUNT(arr))

/* Does path[0..len) walk real edges of g? Returns total weight, or -1. */
static double walk_cost(Graph *g, const int *p, size_t len)
{
    double total = 0;
    for (size_t i = 0; i + 1 < len; ++i) {
        double w;
        if (!graph_has_edge(g, p[i], p[i + 1], &w)) { return -1; }
        total += w;
    }
    return total;
}

static bool proper_coloring(Graph *g, const int *color, size_t n, size_t used)
{
    int top = -1;
    for (size_t u = 0; u < n; ++u) {
        if (!graph_is_alive(g, (int)u)) { continue; }
        top = color[u] > top ? color[u] : top;
        for (size_t v = u + 1; v < n; ++v) {
            if (graph_is_alive(g, (int)v) && graph_has_edge(g, (int)u, (int)v, NULL) && color[u] == color[v]) {
                return false;
            }
        }
    }
    return (size_t)(top + 1) == used;
}

static double zero_h(int v, int t, void *ctx) { (void)v; (void)t; (void)ctx; return 0.0; }
/* Admissible: never overestimates (every edge of U costs >= 1, so hops <= cost). */
static double hop_h(int v, int t, void *ctx) { (void)ctx; return v == t ? 0.0 : 1.0; }
static bool visit_collect(int v, int parent, void *ctx)
{
    (void)parent;
    int *list = ctx;
    list[++list[0]] = v;
    return true;
}

/* ------------------------------------------------------------------ */
/*  The suite (run once per representation)                             */
/* ------------------------------------------------------------------ */

static bool valid_cycle(Graph *g, const int *c, size_t len)
{
    if (len < 3) { return false; }
    for (size_t i = 0; i < len; ++i) {
        if (!graph_has_edge(g, c[i], c[(i + 1) % len], NULL)) { return false; }
    }
    return true;
}

static size_t count_str(const char *hay, const char *needle)
{
    size_t c = 0;
    for (const char *p = hay; (p = strstr(p, needle)) != NULL; p += strlen(needle)) { ++c; }
    return c;
}

/** Runs write(g, tmpfile) and returns the text it produced (caller frees), or NULL. */
static char *capture(GraphStatus (*write)(Graph *, FILE *), Graph *g)
{
    FILE *f = tmpfile();
    if (f == NULL || write(g, f) != GRAPH_OK) { if (f) { fclose(f); } return NULL; }
    fflush(f);
    const long len = ftell(f);
    rewind(f);
    char *text = calloc((size_t)len + 1, 1);
    if (text != NULL && fread(text, 1, (size_t)len, f) != (size_t)len) { free(text); text = NULL; }
    fclose(f);
    return text;
}

static GraphStatus w_graphml(Graph *g, FILE *f) { return graph_write_graphml(g, f); }
static GraphStatus w_gml(Graph *g, FILE *f)     { return graph_write_gml(g, f); }
static GraphStatus w_json(Graph *g, FILE *f)    { return graph_write_json(g, f); }
static GraphStatus w_dot(Graph *g, FILE *f)     { return graph_write_dot(g, f, "U", NULL, true); }

static void suite(void)
{
    char buf[TXT], buf2[TXT];
    Graph *U = MK(false, 8, U_EDGES);
    Graph *D = MK(true, 6, D_EDGES);
    Graph *G = MK(true, 5, G_EDGES);
    Graph *F = MK(true, 4, F_EDGES);
    Graph *NG = MK(true, 4, N_EDGES);
    Graph *CG = MK(true, 3, C_EDGES);
    Graph *B = MK(false, 6, B_EDGES);
    Graph *W = MK(false, 5, W_EDGES);
    Graph *Z = MK(false, 4, Z_EDGES);
    Graph *C6 = graph_gen_cycle(cur_rep, 6, false);
    Graph *DC6 = graph_gen_cycle(cur_rep, 6, true);
    Graph *P5 = graph_gen_path(cur_rep, 5, false);
    Graph *S5 = graph_gen_star(cur_rep, 5);
    const int seven[7] = { 0, 1, 2, 3, 4, 5, 6 };
    Graph *Uc = NULL;                                   /* U without the isolated vertex 7 */
    if (graph_induced_subgraph(U, seven, 7, &Uc) != GRAPH_OK) { fprintf(stderr, "subgraph failed\n"); exit(2); }
    double dv[64], dv2[64];
    int iv[64], iv2[64];
    size_t cnt = 0, cnt2 = 0;
    GraphStatus st;

    /* ================= baseline algorithms (cached by the facade) ================= */
    {
        int order[8], parent[8], comp[8], topo[5];
        const int e_order[] = { 0, 1, 2, 3, 4, 5, 6 };
        chk_int("graph_bfs(U, 0): vertices reached", 7, (long)graph_bfs(U, 0, order, parent));
        chk_ints("graph_bfs(U, 0): visit order", e_order, order, 7);
        chk_int("graph_dfs(U, 0): vertices reached", 7, (long)graph_dfs(U, 0, order, parent));
        const int e_comp[8] = { 0, 0, 0, 0, 0, 0, 0, 1 };
        chk_int("graph_connected_components(U): count", 2, (long)graph_connected_components(U, comp));
        chk_ints("graph_connected_components(U): labels", e_comp, comp, 8);
        chk_bool("graph_has_cycle(U)", true, graph_has_cycle(U));
        chk_bool("graph_has_cycle(P5)", false, graph_has_cycle(P5));
        chk_bool("graph_has_cycle(D)", true, graph_has_cycle(D));
        const bool ok = graph_topological_sort(G, topo);
        int pos[5] = { 0 };
        for (int i = 0; i < 5; ++i) { pos[topo[i]] = i; }
        bool valid = ok;
        for (size_t i = 0; i < COUNT(G_EDGES); ++i) { valid = valid && pos[G_EDGES[i].u] < pos[G_EDGES[i].v]; }
        fmt_ints(buf, sizeof buf, topo, 5);
        chk_prop("graph_topological_sort(G)", "every edge u->v has u before v", valid, "%s", buf);
        graph_clear_error();
        chk_bool("graph_topological_sort(D) on a cyclic graph", false, graph_topological_sort(D, topo));
        chk_status("  ... graph_last_error()", GRAPH_ERR_CYCLE, graph_last_error());
        const double e_dist[8] = { 0, 2, 3, 6, 8, 9, 12, INFINITY };
        chk_bool("graph_dijkstra(U, 0) succeeds", true, graph_dijkstra(U, 0, dv, iv));
        chk_dbls("graph_dijkstra(U, 0): distances", e_dist, dv, 8, 1e-9);
        const size_t plen = graph_path_from_parent(iv, 8, 0, 6, iv2, 8);
        chk_prop("graph_path_from_parent(dijkstra tree, 0 -> 6)", "a path 0..6 of cost 12", plen > 1 && iv2[0] == 0 &&
                 iv2[plen - 1] == 6 && walk_cost(U, iv2, plen) == 12.0, "len %zu, cost %.0f", plen, walk_cost(U, iv2, plen));
        GraphEdge mst[8];
        graph_clear_error();
        chk_dbl("graph_mst_kruskal(U) [disconnected]", INFINITY, graph_mst_kruskal(U, mst, &cnt), 0);
        chk_status("  ... graph_last_error()", GRAPH_ERR_NOT_FOUND, graph_last_error());
        chk_dbl("graph_mst_kruskal(Uc) [connected part]: weight", 13, graph_mst_kruskal(Uc, mst, &cnt), 1e-9);
        chk_int("graph_mst_kruskal(Uc): edges", 6, (long)cnt);
    }

    /* ================= analysis ================= */
    {
        const double e_deg[8] = { 2.0/7, 3.0/7, 3.0/7, 3.0/7, 3.0/7, 3.0/7, 1.0/7, 0 };
        chk_status("graph_degree_centrality(U) status", GRAPH_OK, graph_degree_centrality(U, dv));
        chk_dbls("graph_degree_centrality(U)", e_deg, dv, 8, 1e-9);
        const double e_clo[8] = { 0.395604, 0.514286, 0.514286, 0.571429, 0.571429, 0.514286, 0.342857, 0 };
        (void)graph_closeness_centrality(U, dv);
        chk_dbls("graph_closeness_centrality(U) [Wasserman-Faust]", e_clo, dv, 8, 1e-5);
        const double e_btw[8] = { 0, 2.5, 2.5, 3.5, 3.5, 5, 0, 0 };
        (void)graph_betweenness_centrality(U, dv);
        chk_dbls("graph_betweenness_centrality(U)", e_btw, dv, 8, 1e-9);
        long long tri = -1;
        chk_status("graph_triangle_count(U) status", GRAPH_OK, graph_triangle_count(U, &tri));
        chk_int("graph_triangle_count(U)", 2, (long)tri);
        chk_status("graph_triangle_count(D) on a directed graph", GRAPH_ERR_UNSUPPORTED, graph_triangle_count(D, &tri));
        double avg;
        const double e_cl[8] = { 1, 1.0/3, 1.0/3, 1.0/3, 1.0/3, 1.0/3, 0, 0 };
        (void)graph_clustering_coefficient(U, dv, &avg);
        chk_dbls("graph_clustering_coefficient(U): local", e_cl, dv, 8, 1e-9);
        chk_dbl("graph_clustering_coefficient(U): average", 1.0 / 3.0, avg, 1e-9);
        bool bip = true;
        chk_status("graph_is_bipartite(U) status", GRAPH_OK, graph_is_bipartite(U, &bip, iv));
        chk_bool("graph_is_bipartite(U) [has triangles]", false, bip);
        (void)graph_is_bipartite(B, &bip, iv);
        bool sides = bip;
        for (size_t i = 0; i < COUNT(B_EDGES); ++i) { sides = sides && iv[B_EDGES[i].u] != iv[B_EDGES[i].v]; }
        chk_prop("graph_is_bipartite(B)", "true, endpoints of every edge get different colours", bip && sides,
                 "%s, colours %d%d%d|%d%d%d", bip ? "true" : "false", iv[0], iv[1], iv[2], iv[3], iv[4], iv[5]);
        st = graph_eccentricity(U, iv);
        const int e_ecc[8] = { 4, 3, 3, 2, 2, 3, 4, 0 };
        chk_prop("graph_eccentricity(U)", "[4 3 3 2 2 3 4 0] (hops to the farthest reachable vertex)",
                 st == GRAPH_OK && memcmp(iv, e_ecc, sizeof e_ecc) == 0, "%s", (fmt_ints(buf, sizeof buf, iv, 8), buf));
        int diam = -1, rad = -1;
        chk_status("graph_diameter_radius(U) [not connected]", GRAPH_ERR_UNSUPPORTED, graph_diameter_radius(U, &diam, &rad));
        chk_status("graph_diameter_radius(C6) status", GRAPH_OK, graph_diameter_radius(C6, &diam, &rad));
        chk_int("graph_diameter_radius(C6): diameter", 3, diam);
        chk_int("graph_diameter_radius(C6): radius", 3, rad);
        (void)graph_diameter_radius(P5, &diam, &rad);
        chk_int("graph_diameter_radius(P5): diameter", 4, diam);
        chk_int("graph_diameter_radius(P5): radius", 2, rad);
        Graph *T = NULL;
        chk_status("graph_transpose(D) status", GRAPH_OK, graph_transpose(D, &T));
        bool rev = T != NULL && graph_edge_count(T) == 7;
        for (size_t i = 0; i < COUNT(D_EDGES) && rev; ++i) { rev = graph_has_edge(T, D_EDGES[i].v, D_EDGES[i].u, NULL); }
        chk_prop("graph_transpose(D)", "7 edges, each u->v reversed to v->u", rev, "%zu edges, %s", graph_edge_count(T), rev ? "all reversed" : "mismatch");
        graph_destroy(T);
        chk_status("graph_common_neighbors(U, 1, 2) status", GRAPH_OK, graph_common_neighbors(U, 1, 2, iv, &cnt));
        chk_int("graph_common_neighbors(U, 1, 2): count", 1, (long)cnt);
        chk_int("graph_common_neighbors(U, 1, 2): vertex", 0, iv[0]);
        double jac;
        (void)graph_jaccard(U, 1, 2, &jac);
        chk_dbl("graph_jaccard(U, 1, 2) = |{0}| / |{0,1,2,3,4}|", 0.2, jac, 1e-9);
        (void)graph_jaccard(U, 7, 7, &jac);
        chk_dbl("graph_jaccard(U, 7, 7) [both neighbourhoods empty]", 0.0, jac, 1e-9);
        GraphNetworkStats ns;
        chk_status("graph_network_stats(U) status", GRAPH_OK, graph_network_stats(U, &ns));
        chk_dbl("graph_network_stats(U): density 9 / 28", 9.0 / 28.0, ns.density, 1e-9);
        chk_dbl("graph_network_stats(U): reciprocity (undirected)", 1.0, ns.reciprocity, 1e-9);
        chk_int("graph_network_stats(U): max degree", 3, (long)ns.max_degree);
        size_t hist[5];
        const int e_hist[5] = { 1, 1, 1, 5, 0 };
        (void)graph_degree_histogram(U, hist, 5);
        for (int i = 0; i < 5; ++i) { iv[i] = (int)hist[i]; }
        chk_ints("graph_degree_histogram(U, len 5)", e_hist, iv, 5);
        int seen[8] = { 0 };
        const size_t visited = graph_bfs_visit(U, 0, visit_collect, seen);
        const int e_seen[8] = { 7, 0, 1, 2, 3, 4, 5, 6 };
        chk_int("graph_bfs_visit(U, 0): vertices visited", 7, (long)visited);
        chk_ints("graph_bfs_visit(U, 0): callback order (slot 0 = count)", e_seen, seen, 8);
        memset(seen, 0, sizeof seen);
        chk_int("graph_dfs_visit(U, 0): vertices visited", 7, (long)graph_dfs_visit(U, 0, visit_collect, seen));
    }

    /* ================= connectivity ================= */
    {
        chk_bool("graph_is_connected(U) [vertex 7 is isolated]", false, graph_is_connected(U));
        chk_bool("graph_is_connected(C6)", true, graph_is_connected(C6));
        chk_bool("graph_is_strongly_connected(D)", false, graph_is_strongly_connected(D));
        chk_bool("graph_is_strongly_connected(directed 6-cycle)", true, graph_is_strongly_connected(DC6));
        chk_bool("graph_is_dag(G)", true, graph_is_dag(G));
        chk_bool("graph_is_dag(D)", false, graph_is_dag(D));
        chk_bool("graph_is_tree(P5)", true, graph_is_tree(P5));
        chk_bool("graph_is_tree(U)", false, graph_is_tree(U));
        chk_bool("graph_is_forest(P5)", true, graph_is_forest(P5));
        chk_bool("graph_is_forest(U) [has cycles]", false, graph_is_forest(U));
        const int e_scc[6] = { 3, 3, 3, 2, 1, 0 };
        chk_status("graph_strongly_connected_components(D) status", GRAPH_OK, graph_strongly_connected_components(D, iv, &cnt));
        chk_int("graph_strongly_connected_components(D): count", 4, (long)cnt);
        chk_ints("graph_strongly_connected_components(D): labels (reverse topological)", e_scc, iv, 6);
        bool cut[8];
        chk_status("graph_articulation_points(U) status", GRAPH_OK, graph_articulation_points(U, cut, &cnt));
        chk_int("graph_articulation_points(U): count", 1, (long)cnt);
        for (int i = 0; i < 8; ++i) { iv[i] = cut[i]; }
        const int e_cut[8] = { 0, 0, 0, 0, 0, 1, 0, 0 };
        chk_ints("graph_articulation_points(U): is_cut flags (vertex 5)", e_cut, iv, 8);
        GraphEdge br[9];
        chk_status("graph_bridges(U) status", GRAPH_OK, graph_bridges(U, br, &cnt));
        chk_prop("graph_bridges(U)", "1 bridge: 5-6", cnt == 1 && ((br[0].u == 5 && br[0].v == 6) || (br[0].u == 6 && br[0].v == 5)),
                 "%zu bridge(s), first %d-%d", cnt, cnt ? br[0].u : -1, cnt ? br[0].v : -1);
        GraphGroups gr;
        chk_status("graph_biconnected_components(U) status", GRAPH_OK, graph_biconnected_components(U, &gr));
        size_t big = 0, small = 0;
        for (size_t i = 0; i < gr.count; ++i) {
            const size_t sz = gr.offsets[i + 1] - gr.offsets[i];
            big = sz > big ? sz : big;
            small = (sz < small || small == 0) ? sz : small;
        }
        chk_prop("graph_biconnected_components(U)", "2 blocks: {0,1,2,3,4,5} and {5,6}", gr.count == 2 && big == 6 && small == 2,
                 "%zu blocks, sizes %zu and %zu", gr.count, big, small);
        graph_groups_free(&gr);
        size_t len = graph_find_cycle(U, iv);
        chk_prop("graph_find_cycle(U)", "a real cycle (>= 3 vertices, consecutive vertices adjacent)", valid_cycle(U, iv, len), "length %zu", len);
        chk_int("graph_find_cycle(P5)", 0, (long)graph_find_cycle(P5, iv));
        chk_int("graph_find_cycle(D): length of 0->1->2->0", 3, (long)graph_find_cycle(D, iv));
        Graph *R = NULL;
        chk_status("graph_transitive_closure(D) status", GRAPH_OK, graph_transitive_closure(D, &R));
        chk_prop("graph_transitive_closure(D)", "18 arcs (u -> v iff v reachable from u, u != v)", graph_edge_count(R) == 18,
                 "%zu arcs", graph_edge_count(R));
        const bool r05 = graph_has_edge(R, 0, 5, NULL), r50 = graph_has_edge(R, 5, 0, NULL);
        chk_prop("  ... 0 reaches 5, 5 does not reach 0", "0->5 yes, 5->0 no", r05 && !r50, "0->5 %s, 5->0 %s", r05 ? "yes" : "no", r50 ? "yes" : "no");
        graph_destroy(R);
        chk_status("graph_transitive_closure(G) status", GRAPH_OK, graph_transitive_closure(G, &R));
        chk_int("graph_transitive_closure(G): arcs", 9, (long)graph_edge_count(R));
        graph_destroy(R);
        chk_status("graph_transitive_reduction(G) status", GRAPH_OK, graph_transitive_reduction(G, &R));
        chk_prop("graph_transitive_reduction(G)", "5 edges; the redundant 0->3 is gone", graph_edge_count(R) == 5 && !graph_has_edge(R, 0, 3, NULL),
                 "%zu edges, 0->3 %s", graph_edge_count(R), graph_has_edge(R, 0, 3, NULL) ? "kept" : "removed");
        graph_destroy(R);
    }

    /* ================= shortest paths ================= */
    {
        const double e_bf[4] = { 0, 4, 1, 3 };
        chk_status("graph_bellman_ford(N, 0) status", GRAPH_OK, graph_bellman_ford(NG, 0, dv, iv));
        chk_dbls("graph_bellman_ford(N, 0): distances (negative edge 1->2)", e_bf, dv, 4, 1e-9);
        chk_status("graph_bellman_ford(C, 0) [negative cycle]", GRAPH_ERR_NEGATIVE_CYCLE, graph_bellman_ford(CG, 0, dv, iv));
        static double fw[64];
        static int nx[64];
        chk_status("graph_floyd_warshall(U) status", GRAPH_OK, graph_floyd_warshall(U, fw, nx));
        chk_dbl("graph_floyd_warshall(U): d(0,6)", 12, fw[0 * 8 + 6], 1e-9);
        chk_dbl("graph_floyd_warshall(U): d(1,5)", 7, fw[1 * 8 + 5], 1e-9);
        chk_dbl("graph_floyd_warshall(U): d(0,7) [unreachable]", INFINITY, fw[0 * 8 + 7], 0);
        const size_t plen = graph_path_from_next(nx, 8, 0, 6, iv, 8);
        chk_prop("graph_path_from_next(floyd tree, 0 -> 6)", "a path 0..6 of cost 12", plen > 1 && iv[0] == 0 && iv[plen - 1] == 6 &&
                 walk_cost(U, iv, plen) == 12.0, "len %zu, cost %.0f", plen, walk_cost(U, iv, plen));
        chk_status("graph_floyd_warshall(N) status", GRAPH_OK, graph_floyd_warshall(NG, fw, nx));
        chk_dbl("graph_floyd_warshall(N): d(1,3) = -3 + 2", -1, fw[1 * 4 + 3], 1e-9);
        chk_dbl("graph_floyd_warshall(N): d(0,3) = 4 - 3 + 2", 3, fw[0 * 4 + 3], 1e-9);
        chk_status("graph_johnson(N) status", GRAPH_OK, graph_johnson(NG, fw, nx));
        chk_dbl("graph_johnson(N): d(1,3)", -1, fw[1 * 4 + 3], 1e-9);
        chk_dbl("graph_johnson(N): d(0,3)", 3, fw[0 * 4 + 3], 1e-9);
        chk_status("graph_johnson(C) [negative cycle]", GRAPH_ERR_NEGATIVE_CYCLE, graph_johnson(CG, fw, nx));
        GraphPath p;
        chk_status("graph_a_star(U, 0 -> 6, h = 0) status", GRAPH_OK, graph_a_star(U, 0, 6, zero_h, NULL, &p));
        chk_dbl("graph_a_star(U, 0 -> 6, h = 0): cost", 12, p.cost, 1e-9);
        graph_path_free(&p);
        chk_status("graph_a_star(U, 0 -> 6, h = 1 hop) status", GRAPH_OK, graph_a_star(U, 0, 6, hop_h, NULL, &p));
        chk_prop("graph_a_star(U, 0 -> 6, admissible h)", "optimal: cost 12, path starts at 0 ends at 6",
                 p.cost == 12.0 && p.vertices[0] == 0 && p.vertices[p.length - 1] == 6 && walk_cost(U, p.vertices, p.length) == 12.0,
                 "cost %.0f, %zu vertices", p.cost, p.length);
        graph_path_free(&p);
        chk_status("graph_a_star(U, 0 -> 7) [7 is isolated]", GRAPH_ERR_NOT_FOUND, graph_a_star(U, 0, 7, NULL, NULL, &p));
        const int src[2] = { 0, 6 };
        const double e_ms[8] = { 0, 2, 3, 6, 4, 3, 0, INFINITY };
        chk_status("graph_dijkstra_multi(U, {0, 6}) status", GRAPH_OK, graph_dijkstra_multi(U, src, 2, dv, iv, iv2));
        chk_dbls("graph_dijkstra_multi(U, {0, 6}): distance to nearest source", e_ms, dv, 8, 1e-9);
        chk_prop("graph_dijkstra_multi(U, {0, 6}): origin", "vertex 4 -> source 6, vertex 3 -> source 0", iv2[4] == 6 && iv2[3] == 0,
                 "vertex 4 -> %d, vertex 3 -> %d", iv2[4], iv2[3]);
        const double e_z[4] = { 0, 0, 1, 1 };
        chk_status("graph_zero_one_bfs(Z, 0) status", GRAPH_OK, graph_zero_one_bfs(Z, 0, dv, iv));
        chk_dbls("graph_zero_one_bfs(Z, 0): distances", e_z, dv, 4, 1e-9);
        chk_status("graph_zero_one_bfs(U) [weights are not 0/1]", GRAPH_ERR_INVALID_ARGUMENT, graph_zero_one_bfs(U, 0, dv, iv));
        GraphPath paths[4];
        chk_status("graph_k_shortest_paths(U, 0 -> 5, k = 4) status", GRAPH_OK, graph_k_shortest_paths(U, 0, 5, 4, paths, &cnt));
        chk_int("graph_k_shortest_paths(U, 0 -> 5): paths found", 4, (long)cnt);
        for (size_t i = 0; i < cnt && i < 4; ++i) { dv[i] = paths[i].cost; }
        const double e_k[4] = { 9, 9, 9, 11 };
        chk_dbls("graph_k_shortest_paths(U, 0 -> 5): costs", e_k, dv, 4, 1e-9);
        for (size_t i = 0; i < cnt && i < 4; ++i) { graph_path_free(&paths[i]); }
    }

    /* ================= spanning trees, flow, matching ================= */
    {
        GraphEdge edges[8];
        chk_dbl("graph_mst_prim(Uc): weight", 13, graph_mst_prim(Uc, edges, &cnt), 1e-9);
        chk_int("graph_mst_prim(Uc): edges", 6, (long)cnt);
        chk_dbl("graph_mst_prim(U) [disconnected]", INFINITY, graph_mst_prim(U, edges, &cnt), 0);
        double flow = -1;
        bool side[4];
        chk_status("graph_max_flow(F, 0 -> 3, Dinic) status", GRAPH_OK, graph_max_flow(F, 0, 3, GRAPH_FLOW_DINIC, &flow, side));
        chk_dbl("graph_max_flow(F, 0 -> 3, Dinic): value", 15, flow, 1e-9);
        chk_prop("  ... source side of the minimum cut", "contains 0, not 3", side[0] && !side[3], "0:%d 3:%d", side[0], side[3]);
        chk_status("graph_max_flow(F, 0 -> 3, Edmonds-Karp) status", GRAPH_OK, graph_max_flow(F, 0, 3, GRAPH_FLOW_EDMONDS_KARP, &flow, NULL));
        chk_dbl("graph_max_flow(F, 0 -> 3, Edmonds-Karp): value", 15, flow, 1e-9);
        GraphEdge cut[5];
        chk_status("graph_min_cut(F, 0 -> 3) status", GRAPH_OK, graph_min_cut(F, 0, 3, &flow, cut, &cnt));
        double cut_sum = 0;
        for (size_t i = 0; i < cnt; ++i) { cut_sum += cut[i].weight; }
        chk_dbl("graph_min_cut(F, 0 -> 3): capacity", 15, flow, 1e-9);
        chk_prop("  ... the cut edges add up to the capacity", "sum of cut edge weights = 15", near(cut_sum, 15, 1e-9), "%zu edges, sum %.0f", cnt, cut_sum);
        int match[6];
        chk_status("graph_bipartite_matching(B) status", GRAPH_OK, graph_bipartite_matching(B, match, &cnt));
        chk_int("graph_bipartite_matching(B): size", 3, (long)cnt);
        bool pairs = true;
        for (int v = 0; v < 6; ++v) { pairs = pairs && match[v] >= 0 && match[match[v]] == v && graph_has_edge(B, v, match[v], NULL); }
        chk_prop("  ... the pairs are real edges, each vertex used once", "perfect matching", pairs, "%s", pairs ? "perfect matching" : "invalid");
        chk_status("graph_bipartite_matching(U) [not bipartite]", GRAPH_ERR_UNSUPPORTED, graph_bipartite_matching(U, match, &cnt));
        const double cost[9] = { 4, 1, 3, 2, 0, 5, 3, 2, 2 };
        const int e_asg[3] = { 1, 0, 2 };
        chk_status("graph_hungarian(3x3) status", GRAPH_OK, graph_hungarian(cost, 3, 3, iv, &dv[0]));
        chk_dbl("graph_hungarian(3x3): minimum cost 1 + 2 + 2", 5, dv[0], 1e-9);
        chk_ints("graph_hungarian(3x3): row -> column", e_asg, iv, 3);
    }

    /* ================= combinatorial ================= */
    {
        chk_status("graph_euler_trail(W bowtie, closed) status", GRAPH_OK, graph_euler_trail(W, true, iv, &cnt));
        bool eul = cnt == 7 && iv[0] == iv[6];
        for (size_t i = 0; i + 1 < cnt && eul; ++i) { eul = graph_has_edge(W, iv[i], iv[i + 1], NULL); }
        chk_prop("graph_euler_trail(W, closed)", "7 vertices, closed, every step an edge", eul, "length %zu", cnt);
        chk_status("graph_euler_trail(U) [6 odd-degree vertices]", GRAPH_ERR_NOT_FOUND, graph_euler_trail(U, false, iv, &cnt));
        chk_status("graph_hamiltonian(C6, cycle) status", GRAPH_OK, graph_hamiltonian(C6, true, iv));
        bool ham = true;
        int used[6] = { 0 };
        for (int i = 0; i < 6; ++i) { used[iv[i]]++; ham = ham && graph_has_edge(C6, iv[i], iv[(i + 1) % 6], NULL); }
        for (int i = 0; i < 6; ++i) { ham = ham && used[i] == 1; }
        chk_prop("graph_hamiltonian(C6, cycle)", "visits all 6 vertices once, returns to start", ham, "%s", (fmt_ints(buf, sizeof buf, iv, 6), buf));
        chk_status("graph_hamiltonian(P5, path) status", GRAPH_OK, graph_hamiltonian(P5, false, iv));
        chk_status("graph_hamiltonian(S5 star, path)", GRAPH_ERR_NOT_FOUND, graph_hamiltonian(S5, false, iv));
        const size_t n8 = 8;
        (void)graph_color_greedy(U, iv, &cnt);
        chk_prop("graph_color_greedy(U)", "proper colouring using 3 or 4 colours (chromatic number is 3)", proper_coloring(U, iv, n8, cnt) && cnt >= 3 && cnt <= 4,
                 "%zu colours", cnt);
        (void)graph_color_welsh_powell(U, iv, &cnt);
        chk_prop("graph_color_welsh_powell(U)", "proper colouring using 3 or 4 colours", proper_coloring(U, iv, n8, cnt) && cnt >= 3 && cnt <= 4,
                 "%zu colours", cnt);
        (void)graph_color_dsatur(U, iv, &cnt);
        chk_prop("graph_color_dsatur(U)", "proper colouring using exactly 3 colours (DSatur is optimal here)", proper_coloring(U, iv, n8, cnt) && cnt == 3,
                 "%zu colours", cnt);
        GraphGroups cl;
        chk_status("graph_maximal_cliques(U) status", GRAPH_OK, graph_maximal_cliques(U, &cl));
        chk_int("graph_maximal_cliques(U): count ({0,1,2} {1,3} {2,4} {3,4,5} {5,6} {7})", 6, (long)cl.count);
        graph_groups_free(&cl);
        chk_status("graph_max_clique(U) status", GRAPH_OK, graph_max_clique(U, iv, &cnt));
        chk_int("graph_max_clique(U): size", 3, (long)cnt);
        bool in_set[8];
        chk_status("graph_vertex_cover_approx(U) status", GRAPH_OK, graph_vertex_cover_approx(U, in_set, &cnt));
        bool cover = true;
        for (size_t i = 0; i < COUNT(U_EDGES); ++i) { cover = cover && (in_set[U_EDGES[i].u] || in_set[U_EDGES[i].v]); }
        chk_prop("graph_vertex_cover_approx(U)", "covers every edge, size between optimum 4 and 2 x 4 = 8", cover && cnt >= 4 && cnt <= 8, "size %zu, %s", cnt,
                 cover ? "covers all edges" : "misses an edge");
        chk_status("graph_independent_set_greedy(U) status", GRAPH_OK, graph_independent_set_greedy(U, in_set, &cnt));
        bool indep = true, maximal = true;
        for (size_t i = 0; i < COUNT(U_EDGES); ++i) { indep = indep && !(in_set[U_EDGES[i].u] && in_set[U_EDGES[i].v]); }
        for (int v = 0; v < 8; ++v) {
            bool blocked = in_set[v];
            for (size_t i = 0; i < COUNT(U_EDGES); ++i) {
                blocked = blocked || (U_EDGES[i].u == v && in_set[U_EDGES[i].v]) || (U_EDGES[i].v == v && in_set[U_EDGES[i].u]);
            }
            maximal = maximal && blocked;
        }
        chk_prop("graph_independent_set_greedy(U)", "no two chosen vertices adjacent, and no vertex can be added", indep && maximal,
                 "size %zu, %s%s", cnt, indep ? "independent" : "NOT independent", maximal ? ", maximal" : ", not maximal");
        int par[5], end;
        const double e_lp[5] = { 0, 3, 2, 7, 9 };
        chk_status("graph_dag_longest_path(G) status", GRAPH_OK, graph_dag_longest_path(G, dv, par, &end));
        chk_dbls("graph_dag_longest_path(G): longest path ending at each vertex", e_lp, dv, 5, 1e-9);
        chk_int("graph_dag_longest_path(G): end vertex", 4, end);
        chk_status("graph_dag_longest_path(U) [undirected]", GRAPH_ERR_UNSUPPORTED, graph_dag_longest_path(U, dv, par, &end));
        const int e_cp[4] = { 0, 1, 3, 4 };
        chk_status("graph_critical_path(G) status", GRAPH_OK, graph_critical_path(G, iv, &cnt, &dv[0]));
        chk_ints("graph_critical_path(G): path", e_cp, iv, 4);
        chk_dbl("graph_critical_path(G): cost 3 + 4 + 2", 9, dv[0], 1e-9);
    }

    /* ================= network science ================= */
    {
        size_t its;
        const double e_prU[8] = { 0.097643, 0.127385, 0.161392, 0.204583, 0.141453, 0.179749, 0.066815, 0.020979 };
        chk_status("graph_pagerank(U) status", GRAPH_OK, graph_pagerank(U, 0.85, 1e-12, 500, dv, &its));
        chk_dbls("graph_pagerank(U, damping 0.85; edge weights are transition weights)", e_prU, dv, 8, 1e-5);
        double sum = 0;
        for (int i = 0; i < 8; ++i) { sum += dv[i]; }
        chk_dbl("  ... ranks sum to 1", 1.0, sum, 1e-9);
        const double e_prD[6] = { 0.199678, 0.219147, 0.235696, 0.099507, 0.073587, 0.172385 };
        (void)graph_pagerank(D, 0.85, 1e-12, 500, dv, &its);
        chk_dbls("graph_pagerank(D, damping 0.85; edge weights are transition weights)", e_prD, dv, 6, 1e-5);
        const double e_hub[6] = { 0, 0, 0, 0.850651, 0.525731, 0 }, e_auth[6] = { 0, 0, 0, 0, 0.525731, 0.850651 };
        chk_status("graph_hits(D) status", GRAPH_OK, graph_hits(D, 1e-12, 1000, dv, dv2));
        chk_dbls("graph_hits(D): hub scores", e_hub, dv, 6, 1e-5);
        chk_dbls("graph_hits(D): authority scores", e_auth, dv2, 6, 1e-5);
        const double e_ev[8] = { 0.304254, 0.422261, 0.422261, 0.445560, 0.445560, 0.368926, 0.132912, 0 };
        chk_status("graph_eigenvector_centrality(U) status", GRAPH_OK, graph_eigenvector_centrality(U, 1e-12, 5000, dv));
        chk_dbls("graph_eigenvector_centrality(U): L2-normalised", e_ev, dv, 8, 1e-5);
        const int e_core[8] = { 2, 2, 2, 2, 2, 2, 1, 0 };
        chk_status("graph_core_numbers(U) status", GRAPH_OK, graph_core_numbers(U, iv));
        chk_ints("graph_core_numbers(U)", e_core, iv, 8);
        chk_status("graph_label_propagation(U) status", GRAPH_OK, graph_label_propagation(U, iv, &cnt));
        chk_prop("graph_label_propagation(U)", "vertex 7 alone, 2 or 3 communities in total", cnt >= 2 && cnt <= 3 && iv[7] != iv[0] && iv[7] != iv[6],
                 "%zu communities", cnt);
        double q;
        chk_status("graph_louvain(U) status", GRAPH_OK, graph_louvain(U, iv, &q, &cnt));
        chk_prop("graph_louvain(U)", "modularity >= 0.20 (the optimum is 0.2606)", q >= 0.20 && q <= 0.2607, "modularity %.4f, %zu communities", q, cnt);
        double q2;
        chk_status("graph_modularity(U, {0,1,2}{3,4,5,6}{7}) status", GRAPH_OK, graph_modularity(U, (const int[]){ 0, 0, 0, 1, 1, 1, 1, 2 }, &q2));
        chk_dbl("graph_modularity(U, {0,1,2}{3,4,5,6}{7})", 0.141975, q2, 1e-5);
        double r;
        chk_status("graph_degree_assortativity(U) status", GRAPH_OK, graph_degree_assortativity(U, &r));
        chk_dbl("graph_degree_assortativity(U)", -0.173913, r, 1e-5);
        chk_status("graph_degree_assortativity(C6) [regular graph: undefined]", GRAPH_ERR_UNSUPPORTED, graph_degree_assortativity(C6, &r));
    }

    /* ================= graph operations ================= */
    {
        Graph *X = NULL;
        Graph *cp = graph_copy(U);
        chk_prop("graph_copy(U)", "9 edges, 8 vertices, same edges and weights", cp && graph_edge_count(cp) == 9 && graph_vertex_count(cp) == 8 &&
                 walk_cost(cp, (const int[]){ 0, 1, 2 }, 3) == 3.0, "%zu edges, %zu vertices", graph_edge_count(cp), graph_vertex_count(cp));
        graph_destroy(cp);
        const int tri3[3] = { 0, 1, 2 };
        chk_status("graph_induced_subgraph(U, {0,1,2}) status", GRAPH_OK, graph_induced_subgraph(U, tri3, 3, &X));
        chk_prop("graph_induced_subgraph(U, {0,1,2})", "3 vertices, 3 edges (the triangle)", graph_vertex_count(X) == 3 && graph_edge_count(X) == 3,
                 "%zu vertices, %zu edges", graph_vertex_count(X), graph_edge_count(X));
        graph_destroy(X);
        static const E O_EDGES[] = { {0,1,9}, {0,7,5}, {6,7,1} };
        Graph *O = MK(false, 8, O_EDGES);
        Graph *un = NULL, *in = NULL;
        chk_status("graph_union(U, O) status", GRAPH_OK, graph_union(U, O, &un));
        double w01 = -1;
        (void)graph_has_edge(un, 0, 1, &w01);
        chk_prop("graph_union(U, O)", "11 edges (9 + 2 new); shared edge 0-1 keeps U's weight 2", graph_edge_count(un) == 11 && w01 == 2.0,
                 "%zu edges, w(0-1) = %.0f", graph_edge_count(un), w01);
        chk_status("graph_intersection(U, O) status", GRAPH_OK, graph_intersection(U, O, &in));
        chk_int("graph_intersection(U, O): edges (only 0-1 is in both)", 1, (long)graph_edge_count(in));
        graph_destroy(un); graph_destroy(in); graph_destroy(O);
        chk_status("graph_complement(U) status", GRAPH_OK, graph_complement(U, &X));
        chk_int("graph_complement(U): 28 possible - 9 present", 19, (long)graph_edge_count(X));
        graph_destroy(X);
        /* Relabel U by v -> (3v + 1) mod 8: an isomorphic copy. */
        E perm[COUNT(U_EDGES)];
        for (size_t i = 0; i < COUNT(U_EDGES); ++i) {
            perm[i] = (E){ (3 * U_EDGES[i].u + 1) % 8, (3 * U_EDGES[i].v + 1) % 8, U_EDGES[i].w };
        }
        Graph *PU = build(false, 8, perm, COUNT(perm));
        bool iso = false;
        chk_status("graph_is_isomorphic(U, relabelled U) status", GRAPH_OK, graph_is_isomorphic(U, PU, &iso, iv));
        chk_bool("graph_is_isomorphic(U, relabelled U)", true, iso);
        bool mapping_ok = true;
        for (size_t i = 0; i < COUNT(U_EDGES); ++i) { mapping_ok = mapping_ok && graph_has_edge(PU, iv[U_EDGES[i].u], iv[U_EDGES[i].v], NULL); }
        chk_prop("  ... the returned mapping preserves every edge", "all 9 edges map to edges", mapping_ok, "%s", mapping_ok ? "all map" : "broken");
        (void)graph_remove_edge(PU, perm[0].u, perm[0].v);
        (void)graph_is_isomorphic(U, PU, &iso, NULL);
        chk_bool("graph_is_isomorphic(U, copy minus one edge)", false, iso);
        graph_destroy(PU);
    }

    /* ================= generators ================= */
    {
        Graph *g;
        g = graph_gen_path(cur_rep, 5, false);       chk_int("graph_gen_path(5): edges", 4, (long)graph_edge_count(g));  graph_destroy(g);
        g = graph_gen_cycle(cur_rep, 6, false);      chk_int("graph_gen_cycle(6): edges", 6, (long)graph_edge_count(g)); graph_destroy(g);
        g = graph_gen_complete(cur_rep, 5, false);   chk_int("graph_gen_complete(5): edges", 10, (long)graph_edge_count(g)); graph_destroy(g);
        g = graph_gen_star(cur_rep, 5);              chk_int("graph_gen_star(5): edges", 4, (long)graph_edge_count(g)); graph_destroy(g);
        g = graph_gen_grid(cur_rep, 3, 4);           chk_int("graph_gen_grid(3, 4): edges 3*3 + 2*4", 17, (long)graph_edge_count(g)); graph_destroy(g);
        g = graph_gen_erdos_renyi(cur_rep, 6, 1.0, 1, false, 1, 1); chk_int("graph_gen_erdos_renyi(6, p = 1): edges", 15, (long)graph_edge_count(g)); graph_destroy(g);
        g = graph_gen_erdos_renyi(cur_rep, 6, 0.0, 1, false, 1, 1); chk_int("graph_gen_erdos_renyi(6, p = 0): edges", 0, (long)graph_edge_count(g)); graph_destroy(g);
        g = graph_gen_barabasi_albert(cur_rep, 10, 2, 7);           chk_int("graph_gen_barabasi_albert(10, m = 2): edges 3 + 7*2", 17, (long)graph_edge_count(g)); graph_destroy(g);
        g = graph_gen_random_tree(cur_rep, 20, 3);   chk_bool("graph_gen_random_tree(20): is a tree", true, graph_is_tree(g)); graph_destroy(g);
        g = graph_gen_random_dag(cur_rep, 15, 0.3, 3, 1, 5); chk_bool("graph_gen_random_dag(15, 0.3): is a DAG", true, graph_is_dag(g)); graph_destroy(g);
    }

    /* ================= input / output ================= */
    {
        char *t = capture(w_graphml, U);
        chk_prop("graph_write_graphml(U)", "8 <node> and 9 <edge> elements", t && count_str(t, "<node ") == 8 && count_str(t, "<edge ") == 9, "%zu nodes, %zu edges",
                 t ? count_str(t, "<node ") : 0, t ? count_str(t, "<edge ") : 0);
        free(t);
        t = capture(w_gml, U);
        chk_prop("graph_write_gml(U)", "8 `node [` and 9 `edge [` blocks", t && count_str(t, "node [") == 8 && count_str(t, "edge [") == 9, "%zu nodes, %zu edges",
                 t ? count_str(t, "node [") : 0, t ? count_str(t, "edge [") : 0);
        free(t);
        t = capture(w_json, U);
        chk_prop("graph_write_json(U)", "non-empty JSON object listing 9 edges", t && t[0] == '{' && count_str(t, "\"weight\"") == 9, "%zu weight fields",
                 t ? count_str(t, "\"weight\"") : 0);
        free(t);
        t = capture(w_dot, U);
        chk_prop("graph_write_dot(U)", "`graph \"U\" {` with 9 `--` edges", t && strstr(t, "graph \"U\"") && count_str(t, "--") == 9, "%zu edges", t ? count_str(t, "--") : 0);
        free(t);
        FILE *f = tmpfile();
        chk_status("graph_save_binary(U) status", GRAPH_OK, graph_save_binary(U, f));
        rewind(f);
        Graph *back = NULL;
        chk_status("graph_load_binary(...) status", GRAPH_OK, graph_load_binary(f, cur_rep == GRAPH_LIST ? GRAPH_MATRIX : GRAPH_LIST, &back));
        chk_prop("  ... binary round trip (loaded into the other representation)", "8 vertices, 9 edges, same weights",
                 back && graph_vertex_count(back) == 8 && graph_edge_count(back) == 9 && walk_cost(back, (const int[]){ 0, 1, 3, 4, 5, 6 }, 6) == 2 + 4 + 2 + 1 + 3,
                 "%zu vertices, %zu edges", back ? graph_vertex_count(back) : 0, back ? graph_edge_count(back) : 0);
        graph_destroy(back);
        fclose(f);
        f = tmpfile();
        (void)graph_write_edge_list(U, f);
        rewind(f);
        chk_status("graph_read_edge_list(...) status", GRAPH_OK, graph_read_edge_list(f, false, cur_rep, &back));
        chk_int("  ... edge list round trip: edges", 9, (long)graph_edge_count(back));
        graph_destroy(back);
        fclose(f);
        f = tmpfile();
        fputs("alice bob 2.5\nbob carol\n# a comment\ncarol alice 4\n", f);
        rewind(f);
        char **names = NULL;
        size_t nn = 0;
        chk_status("graph_read_named_edge_list(alice/bob/carol) status", GRAPH_OK, graph_read_named_edge_list(f, false, cur_rep, &back, &names, &nn));
        chk_prop("  ... names and weights", "3 vertices (alice bob carol), 3 edges, alice-bob weighs 2.5",
                 nn == 3 && strcmp(names[0], "alice") == 0 && strcmp(names[2], "carol") == 0 && graph_edge_count(back) == 3,
                 "%zu names (%s %s %s), %zu edges", nn, nn > 0 ? names[0] : "", nn > 1 ? names[1] : "", nn > 2 ? names[2] : "", graph_edge_count(back));
        double wab = 0;
        (void)graph_has_edge(back, 0, 1, &wab);
        chk_dbl("  ... weight of alice-bob", 2.5, wab, 1e-12);
        graph_names_free(names, nn);
        graph_destroy(back);
        fclose(f);
    }

    /* ================= visualization ================= */
    {
        double x[8], y[8];
        chk_status("graph_layout(U) status", GRAPH_OK, graph_layout(U, x, y, NULL));
        bool inside = true;
        for (int v = 0; v < 8; ++v) { inside = inside && x[v] >= 0 && x[v] <= 1 && y[v] >= 0 && y[v] <= 1; }
        chk_prop("graph_layout(U): coordinates", "all 8 vertices inside the unit square", inside, "%s", inside ? "all inside" : "outside");
        VizOptions o = graph_viz_default_options();
        o.title = "showcase graph U (weights shown)";
        o.show_weights = true;
        o.width = o.height = 700;
        const char *svg = cur_rep == GRAPH_LIST ? "build/showcase_U_list.svg" : "build/showcase_U_matrix.svg";
        chk_status("graph_visualize(U, .svg) status", GRAPH_OK, graph_visualize(U, svg, &o));
        FILE *f = fopen(svg, "rb");
        char head[8192] = { 0 };
        const size_t got = f ? fread(head, 1, sizeof head - 1, f) : 0;
        if (f) { fclose(f); }
        chk_prop("  ... the SVG file", "starts with <svg, draws 8 circles", got > 100 && strncmp(head, "<svg", 4) == 0 && count_str(head, "<circle") == 8,
                 "%zu bytes read, %zu circles", got, count_str(head, "<circle"));
        chk_status("graph_visualize(U, .png) status", GRAPH_OK, graph_visualize(U, "build/showcase_U.png", &o));
        graph_clear_error();
        chk_status("graph_visualize(U, unwritable path)", GRAPH_ERR_IO, graph_visualize(U, "/nonexistent_dir/x.png", NULL));
    }

    /* ================= error handling on the unified handle ================= */
    {
        graph_clear_error();
        (void)graph_add_edge(U, 0, 1, 1.0);
        chk_status("graph_add_edge(U, 0, 1) again -> graph_last_error()", GRAPH_ERR_DUPLICATE_EDGE, graph_last_error());
        (void)graph_add_edge(U, 0, 99, 1.0);
        chk_status("graph_add_edge(U, 0, 99) -> graph_last_error()", GRAPH_ERR_INVALID_VERTEX, graph_last_error());
        (void)graph_add_edge(U, 0, 2, NAN);
        chk_status("graph_add_edge(U, weight NaN) -> graph_last_error()", GRAPH_ERR_INVALID_WEIGHT, graph_last_error());
        (void)graph_remove_edge(U, 0, 7);
        chk_status("graph_remove_edge(U, 0, 7) -> graph_last_error()", GRAPH_ERR_EDGE_NOT_FOUND, graph_last_error());
        graph_freeze(U);
        (void)graph_add_edge(U, 0, 7, 1.0);
        chk_status("graph_add_edge on a frozen graph -> graph_last_error()", GRAPH_ERR_FROZEN, graph_last_error());
        graph_thaw(U);
        chk_status("graph_pagerank(NULL)", GRAPH_ERR_INVALID_ARGUMENT, graph_pagerank(NULL, 0.85, 1e-9, 10, dv, NULL));
    }
    (void)buf2; (void)iv2; (void)dv2; (void)cnt2; (void)st; 

    graph_destroy(U); graph_destroy(D); graph_destroy(G); graph_destroy(F); graph_destroy(NG);
    graph_destroy(CG); graph_destroy(B); graph_destroy(W); graph_destroy(Z); graph_destroy(C6);
    graph_destroy(DC6); graph_destroy(P5); graph_destroy(S5); graph_destroy(Uc);
}

/* ------------------------------------------------------------------ */
/*  Output                                                              */
/* ------------------------------------------------------------------ */

int main(void)
{
    for (pass = 0; pass < 2; ++pass) {
        cursor = 0;
        cur_rep = pass == 0 ? GRAPH_LIST : GRAPH_MATRIX;
        suite();
    }
    int failed = 0;
    for (int i = 0; i < nrows; ++i) {
        const Row *r = &rows[i];
        const bool ok = r->ok[0] && r->ok[1];
        failed += !ok;
        printf("%s %s\n", ok ? "[PASS]" : "[FAIL]", r->name);
        printf("         expected: %s\n", r->expected);
        if (strcmp(r->actual[0], r->actual[1]) == 0) {
            printf("         actual  : %s    (list = matrix)\n", r->actual[0]);
        } else {
            printf("         actual  : list   %s\n", r->actual[0]);
            printf("                   matrix %s\n", r->actual[1]);
        }
    }
    printf("\n%d checks, %d passed, %d failed.\n", nrows, nrows - failed, failed);
    printf("Each check ran on a list-backed and a matrix-backed graph; picture of graph U: "
           "build/showcase_U_list.svg, build/showcase_U.png\n");
    return failed == 0 ? 0 : 1;
}
