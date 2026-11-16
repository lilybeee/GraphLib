/*
 * readme_images.c — draws the pictures used in the README, with the library's own visualiser
 * (graph_viz.h).  `make images` writes them to docs/images/.
 *
 *   friends.svg              the 10-person network from quickstart.c, coloured by community
 *   social_communities.png   a social network of 3000 users (the 2400 with friends), coloured by community
 *   social_degree.png        the same network on a blue (few friends) to red (many) temperature scale
 *   social_matrix.png        the same network as an adjacency-matrix heat-map, users sorted by community
 *   grid.png                 a 40 x 40 lattice (structure survives the layout)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "graph_algorithms.h"
#include "graph_viz.h"

static unsigned long long seed = 88172645463325252ULL;
static unsigned rnd(unsigned bound) { seed ^= seed << 13; seed ^= seed >> 7; seed ^= seed << 17; return (unsigned)(seed % bound); }

static void check(const char *what, GraphStatus st)
{
    if (st != GRAPH_OK) {
        fprintf(stderr, "%s failed: %s\n", what, graph_status_string(st));
        exit(1);
    }
}

/* A social network with `groups` communities: most friendships stay inside a person's own group,
 * a few cross to outsiders, and a handful of people are very popular.  Every fifth user is silent
 * (no connections at all), as on a real platform; those users are numbered but not connected. */
static Graph *social_network(size_t users, size_t groups)
{
    Graph *g = graph_create_sized(users, false);
    graph_set_adaptive(g, false);
    for (size_t u = 0; u < users; u++) {
        if (u % 5 == 4) continue;
        const size_t group = u % groups;
        for (int k = 0; k < 5; k++) {
            size_t v;
            if (rnd(100) < 88) {
                /* a friend in the same group; people with a low index are more popular */
                const size_t a = rnd((unsigned)(users / groups)), b = rnd((unsigned)(users / groups));
                v = group + groups * (a < b ? a : b);
            } else {
                v = rnd((unsigned)users);                    /* someone from anywhere */
            }
            if (v != u && v % 5 != 4) graph_add_edge(g, (int)u, (int)v, 1.0);
        }
    }
    return g;
}

/* The same network without the silent users: they would be a ring of isolated dots that hides
 * the structure.  The caller owns the result. */
static Graph *connected_people(Graph *g)
{
    const size_t n = graph_id_bound(g);
    int *keep = malloc(n * sizeof *keep);
    size_t count = 0;
    for (size_t v = 0; v < n; v++) if (graph_degree(g, (int)v) > 0) keep[count++] = (int)v;
    Graph *sub = NULL;
    check("induced_subgraph", graph_induced_subgraph(g, keep, count, &sub));
    free(keep);
    return sub;
}

int main(int argc, char **argv)
{
    const char *dir = argc > 1 ? argv[1] : "docs/images";
    char path[512];

    /* ---- the 10-person network, with names, coloured by community ---- */
    static const char *const NAME[] = {"Ana", "Ben", "Cara", "Dev", "Eli", "Fay", "Gus", "Hana", "Ivo", "Jo"};
    Graph *small = graph_create_sized(10, false);
    const int fr[][2] = {{0, 1}, {0, 2}, {1, 2}, {2, 3}, {4, 5}, {4, 6}, {5, 6}, {6, 7}, {3, 4}, {7, 8}};
    for (size_t i = 0; i < sizeof fr / sizeof *fr; i++) graph_add_edge(small, fr[i][0], fr[i][1], 1.0);
    VizOptions o = graph_viz_default_options();
    o.width = o.height = 640;
    o.color  = VIZ_COLOR_COMMUNITY;
    o.labels = NAME;
    o.title  = "Ten people, one bridge (Dev-Eli), and Jo who knows nobody";
    snprintf(path, sizeof path, "%s/friends.svg", dir);
    check("friends.svg", graph_visualize(small, path, &o));
    graph_destroy(small);

    /* ---- a 3000-user network, three views ---- */
    Graph *full = social_network(3000, 8);
    Graph *net  = connected_people(full);
    printf("social network: %zu users, %zu with friends, %zu friendships\n", graph_vertex_count(full),
           graph_vertex_count(net), graph_edge_count(net));
    graph_destroy(full);

    o = graph_viz_default_options();
    o.width = o.height = 1100;
    o.format = VIZ_FORMAT_PNG;
    o.color  = VIZ_COLOR_COMMUNITY;
    o.group_name = "Community";
    o.title = "2400 people with friends, coloured by community (Louvain)";
    snprintf(path, sizeof path, "%s/social_communities.png", dir);
    check("social_communities.png", graph_visualize(net, path, &o));

    o.color = VIZ_COLOR_DEGREE;
    o.title = "The same network: blue = few friends, red = very popular";
    snprintf(path, sizeof path, "%s/social_degree.png", dir);
    check("social_degree.png", graph_visualize(net, path, &o));

    /* The matrix view draws vertices in ID order, so first renumber people by community: members
     * of a group then sit next to each other and each group shows up as a block on the diagonal. */
    const size_t n = graph_id_bound(net);
    int *community = malloc(n * sizeof *community), *new_id = malloc(n * sizeof *new_id);
    size_t count;
    check("louvain", graph_louvain(net, community, NULL, &count));
    size_t next = 0;
    for (size_t c = 0; c < count; c++)
        for (size_t v = 0; v < n; v++) if ((size_t)community[v] == c) new_id[v] = (int)next++;
    Graph *sorted = graph_create_sized(n, false);
    graph_set_adaptive(sorted, false);
    for (size_t v = 0; v < n; v++) {
        int nb[256];
        const size_t d = graph_neighbors(net, (int)v, nb, NULL, 256);
        for (size_t i = 0; i < d; i++) if ((int)v < nb[i]) graph_add_edge(sorted, new_id[v], new_id[nb[i]], 1.0);
    }
    free(community);
    free(new_id);
    o.layout = VIZ_LAYOUT_MATRIX;
    o.color  = VIZ_COLOR_AUTO;
    o.title  = "The same network as an adjacency matrix, people grouped by community";
    snprintf(path, sizeof path, "%s/social_matrix.png", dir);
    check("social_matrix.png", graph_visualize(sorted, path, &o));
    graph_destroy(sorted);
    graph_destroy(net);

    /* ---- a lattice ---- */
    Graph *grid = graph_gen_grid(GRAPH_LIST, 40, 40);
    o = graph_viz_default_options();
    o.width = o.height = 700;
    o.format = VIZ_FORMAT_PNG;
    o.color  = VIZ_COLOR_NONE;
    o.title  = "A 40 x 40 lattice: structure survives the layout";
    snprintf(path, sizeof path, "%s/grid.png", dir);
    check("grid.png", graph_visualize(grid, path, &o));
    graph_destroy(grid);
    return 0;
}
