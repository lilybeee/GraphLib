/*
 * quickstart.c — the first program to read.  Builds a small friendship network, then asks the
 * library the questions you would ask of any network: who is connected to whom, how are two
 * people linked, who is most influential, which groups exist, whom should I be introduced to.
 *
 *   make examples          # builds and runs this and the other examples
 *   cc -Iinclude examples/quickstart.c build/libgraph.a -lm -o quickstart && ./quickstart
 *
 * Vertices are plain integers (0, 1, 2, ...).  The library never sees names; the array `NAME`
 * below is how this program remembers who is who (the mapping layer in mapping.h does the same
 * for arbitrary keys).
 */

#include <stdio.h>

#include "graph_algorithms.h"   /* the Graph handle plus every algorithm */
#include "graph_social.h"       /* per-person queries: friends of friends, chains of friends */

static const char *NAME[] = {"Ana", "Ben", "Cara", "Dev", "Eli", "Fay", "Gus", "Hana", "Ivo", "Jo"};

enum { ANA, BEN, CARA, DEV, ELI, FAY, GUS, HANA, IVO, JO, PEOPLE };

int main(void)
{
    /* 1. Create an undirected graph with 10 people.  By default it starts as an adjacency list
     *    and may switch to a matrix by itself if the way you use it makes that cheaper. */
    Graph *g = graph_create_sized(PEOPLE, false);

    /* 2. Friendships.  Two groups of four, one bridge between them, one extra friend, and Jo,
     *    who has no connections at all. */
    const int friendships[][2] = {
        {ANA, BEN}, {ANA, CARA}, {BEN, CARA}, {CARA, DEV},     /* group one            */
        {ELI, FAY}, {ELI, GUS}, {FAY, GUS}, {GUS, HANA},       /* group two            */
        {DEV, ELI},                                            /* the bridge           */
        {HANA, IVO},                                           /* Ivo knows only Hana  */
    };
    for (size_t i = 0; i < sizeof friendships / sizeof *friendships; i++) {
        graph_add_edge(g, friendships[i][0], friendships[i][1], 1.0);
    }
    printf("%zu people, %zu friendships, stored as a %s\n", graph_vertex_count(g),
           graph_edge_count(g), graph_rep_name(graph_representation(g)));

    /* 3. Connected components: groups of people who can reach each other. */
    int component[PEOPLE];
    const size_t groups = graph_connected_components(g, component);
    printf("\nconnected groups: %zu\n", groups);
    for (size_t c = 0; c < groups; c++) {
        printf("  group %zu:", c);
        for (int v = 0; v < PEOPLE; v++) {
            if ((size_t)component[v] == c) printf(" %s", NAME[v]);
        }
        printf("\n");
    }

    /* 4. How are Ana and Ivo connected?  A shortest chain of friends. */
    GraphPath chain;
    if (graph_social_path(g, ANA, IVO, SIZE_MAX, &chain) == GRAPH_OK) {
        printf("\nshortest chain Ana -> Ivo (%.0f hops):", chain.cost);
        for (size_t i = 0; i < chain.length; i++) printf(" %s", NAME[chain.vertices[i]]);
        printf("\n");
        graph_path_free(&chain);               /* results that own memory have a *_free function */
    }
    /* Jo knows nobody, so there is no chain, and the library says so instead of guessing. */
    if (graph_social_path(g, ANA, JO, SIZE_MAX, &chain) != GRAPH_OK) {
        printf("chain Ana -> Jo: %s\n", graph_status_string(graph_last_error()));
    }

    /* 5. Who is most influential?  PageRank. */
    double rank[PEOPLE];
    graph_pagerank(g, 0.85, 1e-9, 100, rank, NULL);
    int best = 0;
    for (int v = 1; v < PEOPLE; v++) if (rank[v] > rank[best]) best = v;
    printf("\nmost influential (PageRank): %s (%.3f)\n", NAME[best], rank[best]);

    /* 6. Which communities exist?  Louvain modularity optimisation. */
    int community[PEOPLE];
    double modularity;
    size_t ncomm;
    graph_louvain(g, community, &modularity, &ncomm);
    printf("\ncommunities found: %zu (modularity %.3f)\n", ncomm, modularity);
    for (size_t c = 0; c < ncomm; c++) {
        printf("  community %zu:", c);
        for (int v = 0; v < PEOPLE; v++) {
            if ((size_t)community[v] == c) printf(" %s", NAME[v]);
        }
        printf("\n");
    }

    /* 7. Whom should Ben be introduced to?  People who share the most friends with him. */
    GraphSuggestion suggestion[3];
    size_t count;
    graph_social_suggest_friends(g, BEN, 3, suggestion, &count);
    printf("\npeople Ben may know:\n");
    for (size_t i = 0; i < count; i++) {
        printf("  %s (%zu mutual friend%s)\n", NAME[suggestion[i].user], suggestion[i].mutual,
               suggestion[i].mutual == 1 ? "" : "s");
    }

    /* 8. Errors are values, never crashes: a bad vertex is reported, and the graph is untouched. */
    if (!graph_add_edge(g, ANA, 99, 1.0)) {
        printf("\nadd_edge(Ana, 99) failed: %s\n", graph_status_string(graph_last_error()));
    }

    graph_destroy(g);
    return 0;
}
