/*
 * gallery.c — draws one picture per kind of graph and per size, with the library's own
 * visualiser (graph_viz.h), so you can see what the structures the library works on look like.
 * `make gallery` writes them to docs/images/gallery/ and prints a line for each.
 *
 *   small  (SVG, labelled)   path, cycle, star, complete graph, binary tree, grid, Petersen graph,
 *                            complete bipartite graph, layered DAG, ring of cliques
 *   medium (PNG)             grid, random tree, random graph (Erdos-Renyi), small world, ring of cliques
 *   growing                  the same preferential-attachment (hub) graph at 100 / 1 000 / 10 000 / 50 000 users
 *   large  (PNG)             a 100 x 100 grid; a 100 000-user social network (only the connected people)
 *   matrices                 the adjacency matrix of five different kinds of graph, side by side in structure
 */

#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "graph_algorithms.h"
#include "graph_viz.h"

static const char *g_dir = "docs/images/gallery";
static unsigned long long seed = 88172645463325252ULL;
static unsigned rnd(unsigned bound) { seed ^= seed << 13; seed ^= seed >> 7; seed ^= seed << 17; return (unsigned)(seed % bound); }

static double now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec + (double)t.tv_nsec * 1e-9;
}

static void must(const char *what, GraphStatus st)
{
    if (st != GRAPH_OK) {
        fprintf(stderr, "%s failed: %s\n", what, graph_status_string(st));
        exit(1);
    }
}

static Graph *blank(size_t n, bool directed)
{
    Graph *g = graph_create_sized(n, directed);
    if (g == NULL) { fprintf(stderr, "out of memory\n"); exit(1); }
    graph_set_adaptive(g, false);
    return g;
}

/* Draw g to <dir>/<name>.<ext>, print what was drawn and how long it took, then free it. */
static void draw(Graph *g, const char *name, const char *ext, int size, VizColor color, const char *title)
{
    char path[512];
    snprintf(path, sizeof path, "%s/%s.%s", g_dir, name, ext);
    if (strlen(title) > (size_t)size / 17) fprintf(stderr, "WARNING: title of %s is too long for %d px\n", name, size);
    VizOptions o = graph_viz_default_options();
    o.width = o.height = size;
    o.format = strcmp(ext, "svg") == 0 ? VIZ_FORMAT_SVG : VIZ_FORMAT_PNG;
    o.color  = color;
    o.title  = title;
    o.group_name = "Community";
    const double t = now();
    must(name, graph_visualize(g, path, &o));
    printf("  %-28s %7zu vertices %8zu edges  %5.1f s   %s\n", name, graph_vertex_count(g),
           graph_edge_count(g), now() - t, title);
    graph_destroy(g);
}

/* ---- constructions the generators do not provide ---- */

static Graph *binary_tree(size_t levels)               /* complete binary tree: node i has children 2i+1, 2i+2 */
{
    const size_t n = ((size_t)1 << levels) - 1;
    Graph *g = blank(n, false);
    for (size_t i = 1; i < n; i++) graph_add_edge(g, (int)i, (int)((i - 1) / 2), 1.0);
    return g;
}

static Graph *petersen(void)
{
    Graph *g = blank(10, false);
    for (int i = 0; i < 5; i++) {
        graph_add_edge(g, i, (i + 1) % 5, 1.0);                /* outer pentagon */
        graph_add_edge(g, i, i + 5, 1.0);                      /* spokes */
        graph_add_edge(g, 5 + i, 5 + (i + 2) % 5, 1.0);        /* inner pentagram */
    }
    return g;
}

static Graph *complete_bipartite(size_t a, size_t b)
{
    Graph *g = blank(a + b, false);
    for (size_t i = 0; i < a; i++) for (size_t j = 0; j < b; j++) graph_add_edge(g, (int)i, (int)(a + j), 1.0);
    return g;
}

static Graph *layered_dag(size_t layers, size_t width)
{
    Graph *g = blank(layers * width, true);
    for (size_t l = 0; l + 1 < layers; l++)
        for (size_t i = 0; i < width; i++)
            for (size_t j = 0; j < width; j++)
                if ((i + j + l) % 2 == 0 || i == j) graph_add_edge(g, (int)(l * width + i), (int)((l + 1) * width + j), 1.0);
    return g;
}

static Graph *ring_of_cliques(size_t cliques, size_t size)
{
    Graph *g = blank(cliques * size, false);
    for (size_t c = 0; c < cliques; c++) {
        for (size_t a = 0; a < size; a++)
            for (size_t b = a + 1; b < size; b++) graph_add_edge(g, (int)(c * size + a), (int)(c * size + b), 1.0);
        graph_add_edge(g, (int)(c * size), (int)(((c + 1) % cliques) * size + size - 1), 1.0);
    }
    return g;
}

static Graph *small_world(size_t n, size_t k, unsigned rewire_percent)
{
    Graph *g = blank(n, false);
    for (size_t i = 0; i < n; i++)
        for (size_t d = 1; d <= k; d++) {
            size_t j = (i + d) % n;
            if (rnd(100) < rewire_percent) j = rnd((unsigned)n);          /* a shortcut across the ring */
            if (j != i) graph_add_edge(g, (int)i, (int)j, 1.0);
        }
    return g;
}

/* A social network: preferential attachment, every fifth user silent (never connected). */
static Graph *social(size_t users)
{
    Graph *g = blank(users, false);
    int *ends = malloc(users * 24 * sizeof *ends);
    size_t en = 0;
    for (size_t v = 0; v < users; v++) {
        if (v % 5 == 4) continue;
        const size_t links = 1 + rnd(8);
        for (size_t k = 0; k < links && en > 0; k++) {
            const int f = ends[rnd((unsigned)en)];
            if (f != (int)v && graph_add_edge(g, (int)v, f, 1.0)) { ends[en++] = (int)v; ends[en++] = f; }
        }
        if (en == 0) { ends[en++] = (int)v; ends[en++] = (int)v; }
    }
    free(ends);
    return g;
}

/* Only the people who have at least one friend (silent users would be a ring of isolated dots). */
static Graph *connected_people(Graph *g)
{
    const size_t n = graph_id_bound(g);
    int *keep = malloc(n * sizeof *keep);
    size_t count = 0;
    for (size_t v = 0; v < n; v++) if (graph_degree(g, (int)v) > 0) keep[count++] = (int)v;
    Graph *sub = NULL;
    must("induced_subgraph", graph_induced_subgraph(g, keep, count, &sub));
    free(keep);
    return sub;
}

/* Renumber so members of a community are adjacent; then each community is a block on the diagonal. */
static Graph *sorted_by_community(Graph *g)
{
    const size_t n = graph_id_bound(g);
    int *comm = malloc(n * sizeof *comm), *id = malloc(n * sizeof *id);
    size_t groups = 0;
    must("louvain", graph_louvain(g, comm, NULL, &groups));
    size_t next = 0;
    for (size_t c = 0; c < groups; c++) for (size_t v = 0; v < n; v++) if ((size_t)comm[v] == c) id[v] = (int)next++;
    Graph *out = blank(n, graph_is_directed(g));
    for (size_t v = 0; v < n; v++) {
        int nb[512];
        const size_t d = graph_neighbors(g, (int)v, nb, NULL, 512);
        for (size_t i = 0; i < d; i++) if ((int)v < nb[i]) graph_add_edge(out, id[v], id[nb[i]], 1.0);
    }
    free(comm); free(id);
    return out;
}

static void matrix_view(Graph *g, const char *name, const char *title)
{
    char path[512];
    snprintf(path, sizeof path, "%s/%s.png", g_dir, name);
    if (strlen(title) > 640 / 17) fprintf(stderr, "WARNING: title of %s is too long\n", name);
    VizOptions o = graph_viz_default_options();
    o.width = o.height = 640;
    o.format = VIZ_FORMAT_PNG;
    o.layout = VIZ_LAYOUT_MATRIX;
    o.title  = title;
    must(name, graph_visualize(g, path, &o));
    printf("  %-28s %7zu vertices %8zu edges             %s\n", name, graph_vertex_count(g), graph_edge_count(g), title);
    graph_destroy(g);
}

int main(int argc, char **argv)
{
    if (argc > 1) g_dir = argv[1];
    const double t0 = now();

    printf("small graphs (SVG, every vertex labelled)\n");
    draw(graph_gen_path(GRAPH_LIST, 12, false), "path_12", "svg", 640, VIZ_COLOR_NONE, "Path: 12 vertices, all bridges");
    draw(graph_gen_cycle(GRAPH_LIST, 12, false), "cycle_12", "svg", 640, VIZ_COLOR_NONE, "Cycle: 12 vertices, no bridges");
    draw(graph_gen_star(GRAPH_LIST, 12), "star_12", "svg", 640, VIZ_COLOR_DEGREE, "Star: 1 hub, 11 leaves");
    draw(graph_gen_complete(GRAPH_LIST, 8, false), "complete_8", "svg", 640, VIZ_COLOR_NONE, "Complete graph K8: 28 edges");
    draw(binary_tree(5), "binary_tree_31", "svg", 800, VIZ_COLOR_DEGREE, "Complete binary tree: 31 vertices");
    draw(graph_gen_grid(GRAPH_LIST, 6, 6), "grid_6x6", "svg", 640, VIZ_COLOR_NONE, "Grid 6 x 6: 36 vertices");
    draw(petersen(), "petersen_10", "svg", 640, VIZ_COLOR_NONE, "Petersen graph: 3 friends each");
    draw(complete_bipartite(3, 4), "bipartite_3_4", "svg", 640, VIZ_COLOR_NONE, "Complete bipartite K3,4");
    draw(layered_dag(4, 4), "dag_4x4", "svg", 800, VIZ_COLOR_NONE, "DAG in 4 layers (arrows point forward)");
    draw(ring_of_cliques(4, 5), "ring_of_cliques_4x5", "svg", 800, VIZ_COLOR_COMMUNITY, "Four groups of five, joined in a ring");

    printf("medium graphs (PNG)\n");
    draw(graph_gen_grid(GRAPH_LIST, 30, 30), "grid_30x30", "png", 800, VIZ_COLOR_NONE, "Grid 30 x 30: 900 vertices");
    draw(graph_gen_random_tree(GRAPH_LIST, 500, 5), "random_tree_500", "png", 900, VIZ_COLOR_DEGREE, "Random tree: 500 vertices, no cycles");
    /* About 8 friends per vertex keeps a random graph in one piece; with fewer, a few loose vertices
     * make the layout push the main cluster into a corner. */
    draw(graph_gen_erdos_renyi(GRAPH_LIST, 400, 0.02, 7, false, 1, 1), "random_graph_400", "png", 900, VIZ_COLOR_COMMUNITY,
         "Random graph: a hairball, no clusters");
    draw(small_world(500, 3, 8), "small_world_500", "png", 900, VIZ_COLOR_COMMUNITY, "Small world: ring + 8% random shortcuts");
    {
        Graph *ring = ring_of_cliques(20, 6);
        int comm[120];
        size_t found = 0;
        must("louvain", graph_louvain(ring, comm, NULL, &found));
        char title[96];
        snprintf(title, sizeof title, "20 groups of 6; Louvain finds %zu", found);   /* measured, not assumed */
        draw(ring, "ring_of_cliques_20x6", "png", 900, VIZ_COLOR_COMMUNITY, title);
    }

    printf("one kind of graph, growing: preferential attachment (a few very popular hubs)\n");
    const size_t sizes[] = {100, 1000, 10000, 50000};
    for (size_t i = 0; i < 4; i++) {
        char name[64], title[128];
        snprintf(name, sizeof name, "hubs_%zu", sizes[i]);
        snprintf(title, sizeof title, "Hubs: %zu users, 2 friends on joining", sizes[i]);
        draw(graph_gen_barabasi_albert(GRAPH_LIST, sizes[i], 2, 3), name, "png", 1000, VIZ_COLOR_DEGREE, title);
    }

    printf("large graphs (PNG)\n");
    draw(graph_gen_grid(GRAPH_LIST, 100, 100), "grid_100x100", "png", 1000, VIZ_COLOR_NONE, "Grid 100 x 100: 10 000 vertices");
    {
        Graph *all = social(100000);
        Graph *people = connected_people(all);
        printf("  (social network: %zu users, %zu of them with at least one friend)\n", graph_vertex_count(all), graph_vertex_count(people));
        graph_destroy(all);
        draw(people, "social_100000", "png", 1400, VIZ_COLOR_DEGREE,
             "Social network: 100 000 users, 80 000 drawn (red = popular)");
    }

    printf("adjacency matrices: the same picture tells the type of graph\n");
    matrix_view(graph_gen_erdos_renyi(GRAPH_LIST, 400, 0.02, 7, false, 1, 1), "matrix_random", "Random graph: no structure");
    matrix_view(graph_gen_barabasi_albert(GRAPH_LIST, 400, 3, 3), "matrix_hubs", "Hubs: stripes = popular users");
    matrix_view(ring_of_cliques(20, 10), "matrix_cliques", "Ring of cliques: diagonal blocks");
    matrix_view(graph_gen_grid(GRAPH_LIST, 20, 20), "matrix_grid", "Grid 20 x 20: thin bands");
    {
        Graph *net = connected_people(social(2000));
        matrix_view(sorted_by_community(net), "matrix_social", "Social, sorted by community");
        graph_destroy(net);
    }

    printf("done in %.1f s\n", now() - t0);
    return 0;
}
