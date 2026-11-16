/*
 * bench_matrix.c — timing and memory sweep for the matrix algorithms.
 *
 * Prints CSV on stdout:  algorithm,n,density,edges,time_ns,matrix_bytes
 * Build with optimisation for meaningful numbers:
 *     make bench CFLAGS="-std=c11 -O2"
 * Each timing is the best of three runs.
 */

#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "matrix_algorithms.h"

static double now_ns(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec * 1e9 + (double)ts.tv_nsec;
}

typedef struct {
    const MatrixGraph *g;
    const MatrixGraph *dag;
    size_t             n;
    int               *ia;      /* scratch of n ints   */
    int               *ib;
    double            *da;      /* scratch of n doubles */
    bool              *ba;
    MatrixEdge        *edges;
} Ctx;

static void run_bfs(Ctx *c)        { (void)matrix_bfs(c->g, 0, c->ia, c->ib); }
static void run_dfs(Ctx *c)        { (void)matrix_dfs(c->g, 0, c->ia, c->ib); }
static void run_components(Ctx *c) { (void)matrix_connected_components(c->g, c->ia); }
static void run_dijkstra(Ctx *c)   { (void)matrix_dijkstra(c->g, 0, c->da, c->ib); }
static void run_kruskal(Ctx *c)    { size_t k; (void)matrix_mst_kruskal(c->g, c->edges, &k); }
static void run_prim(Ctx *c)       { size_t k; (void)matrix_mst_prim(c->g, c->edges, &k); }
static void run_bellman(Ctx *c)    { (void)matrix_bellman_ford(c->g, 0, c->da, c->ib); }
static void run_scc(Ctx *c)        { (void)matrix_strongly_connected_components(c->dag, c->ia, NULL); }
static void run_bridges(Ctx *c)    { size_t k; (void)matrix_bridges(c->g, c->edges, &k); }
static void run_cut_vertices(Ctx *c) { size_t k; (void)matrix_articulation_points(c->g, c->ba, &k); }
static void run_pagerank(Ctx *c)   { (void)matrix_pagerank(c->g, 0.85, 1e-9, 100, c->da, NULL); }
static void run_core(Ctx *c)       { (void)matrix_core_numbers(c->g, c->ia); }
static void run_coloring(Ctx *c)   { (void)matrix_color_dsatur(c->g, c->ia, NULL); }
static void run_louvain(Ctx *c)    { (void)matrix_louvain(c->g, c->ia, NULL, NULL); }
static void run_triangles(Ctx *c)  { long long t; (void)matrix_triangle_count(c->g, &t); }
static void run_flow(Ctx *c)       { double v; (void)matrix_max_flow(c->dag, 0, (int)c->n - 1, MATRIX_FLOW_DINIC, &v, NULL); }

typedef struct {
    const char *name;
    void      (*fn)(Ctx *);
    size_t      max_n;        /* skip when the graph is larger than this */
} Case;

int main(void)
{
    static const Case cases[] = {
        { "bfs", run_bfs, 100000 },           { "dfs", run_dfs, 100000 },
        { "components", run_components, 100000 }, { "dijkstra", run_dijkstra, 100000 },
        { "kruskal", run_kruskal, 100000 },   { "prim", run_prim, 4000 },
        { "bellman_ford", run_bellman, 1000 }, { "scc", run_scc, 100000 },
        { "bridges", run_bridges, 100000 },   { "cut_vertices", run_cut_vertices, 100000 },
        { "pagerank", run_pagerank, 100000 }, { "core_numbers", run_core, 100000 },
        { "dsatur", run_coloring, 4000 },     { "louvain", run_louvain, 4000 },
        { "triangles", run_triangles, 2000 }, { "dinic", run_flow, 2000 },
    };
    static const size_t sizes[]      = { 200, 1000, 2000 };
    static const double densities[]  = { 0.001, 0.01, 0.1 };

    puts("algorithm,n,density,edges,time_ns,matrix_bytes");
    for (size_t si = 0; si < sizeof(sizes) / sizeof(sizes[0]); ++si) {
        for (size_t di = 0; di < sizeof(densities) / sizeof(densities[0]); ++di) {
            const size_t n = sizes[si];
            MatrixGraph *g   = matrix_gen_erdos_renyi(n, densities[di], 12345, false, 1.0, 100.0);
            MatrixGraph *dag = matrix_gen_random_dag(n, densities[di], 12345, 1.0, 100.0);
            if (g == NULL || dag == NULL) {
                fputs("allocation failed\n", stderr);
                return 1;
            }
            Ctx c = { g, dag, n, malloc(n * sizeof(int)), malloc(n * sizeof(int)),
                      malloc(n * sizeof(double)), malloc(n * sizeof(bool)),
                      malloc(n * sizeof(MatrixEdge)) };
            if (!c.ia || !c.ib || !c.da || !c.ba || !c.edges) {
                fputs("allocation failed\n", stderr);
                return 1;
            }
            for (size_t k = 0; k < sizeof(cases) / sizeof(cases[0]); ++k) {
                if (n > cases[k].max_n) {
                    continue;
                }
                double best = 1e30;
                for (int rep = 0; rep < 3; ++rep) {
                    const double t0 = now_ns();
                    cases[k].fn(&c);
                    const double dt = now_ns() - t0;
                    if (dt < best) {
                        best = dt;
                    }
                }
                printf("%s,%zu,%g,%zu,%.0f,%zu\n", cases[k].name, n, densities[di],
                       matrix_graph_edge_count(g), best, matrix_graph_memory_bytes(g));
            }
            free(c.ia); free(c.ib); free(c.da); free(c.ba); free(c.edges);
            matrix_graph_destroy(g);
            matrix_graph_destroy(dag);
        }
    }
    return 0;
}
