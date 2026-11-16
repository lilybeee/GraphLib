/*
 * algorithms.c — the everyday algorithms, each on a graph small enough to check by hand.
 *
 *   make examples
 *
 * Every function here is declared in graph_algorithms.h and works on a Graph, whichever
 * representation it is stored in.  Arrays you pass in are sized by graph_id_bound() (one entry
 * per vertex ID) unless the documentation of that function says otherwise.
 */

#include <math.h>
#include <stdio.h>

#include "graph_algorithms.h"

int main(void)
{
    /* ---- A weighted road map: vertices A..G are 0..6, weights are kilometres ---- */
    const char *town = "ABCDEFG";
    Graph *map = graph_create_sized(7, false);
    const struct { int a, b; double km; } road[] = {
        {0, 1, 7}, {0, 3, 5}, {1, 2, 8}, {1, 3, 9}, {1, 4, 7}, {2, 4, 5},
        {3, 4, 15}, {3, 5, 6}, {4, 5, 8}, {4, 6, 9}, {5, 6, 11}};
    for (size_t i = 0; i < sizeof road / sizeof *road; i++) graph_add_edge(map, road[i].a, road[i].b, road[i].km);

    /* Shortest paths from A (Dijkstra).  `parent` lets you rebuild the route itself. */
    double dist[7];
    int parent[7];
    graph_dijkstra(map, 0, dist, parent);
    printf("Dijkstra from A:");
    for (int v = 0; v < 7; v++) printf("  %c=%.0f", town[v], dist[v]);
    int route[7];
    const size_t hops = graph_path_from_parent(parent, 7, 0, 6, route, 7);
    printf("\n  route A -> G:");
    for (size_t i = 0; i < hops; i++) printf(" %c", town[route[i]]);

    /* Minimum spanning tree: the cheapest set of roads that keeps every town connected. */
    GraphEdge tree[7];
    size_t n_tree = 0;
    const double total = graph_mst_kruskal(map, tree, &n_tree);
    printf("\n\nminimum spanning tree: %.0f km using %zu roads:", total, n_tree);
    for (size_t i = 0; i < n_tree; i++) printf(" %c-%c", town[tree[i].u], town[tree[i].v]);

    /* ---- Ordering tasks: a directed acyclic graph of "must happen before" ---- */
    const char *task[] = {"design", "build", "test", "docs", "release"};
    Graph *plan = graph_create_sized(5, true);
    const int before[][2] = {{0, 1}, {0, 3}, {1, 2}, {2, 4}, {3, 4}};
    for (size_t i = 0; i < 5; i++) graph_add_edge(plan, before[i][0], before[i][1], 1.0);
    int order[5];
    if (graph_topological_sort(plan, order)) {
        printf("\n\nan order that respects every dependency:");
        for (int i = 0; i < 5; i++) printf(" %s", task[order[i]]);
    }
    graph_add_edge(plan, 4, 0, 1.0);                         /* release before design: a cycle */
    /* Make the call first and read the reason afterwards: the order in which C evaluates the
     * arguments of one printf is unspecified, so calling it inside the printf would be a bug. */
    const bool sorted = graph_topological_sort(plan, order);
    printf("\nafter adding release -> design: topological sort %s (%s)", sorted ? "works" : "fails",
           graph_status_string(graph_last_error()));

    /* ---- Strongly connected components: groups that can all reach each other ---- */
    Graph *links = graph_create_sized(6, true);
    const int arrow[][2] = {{0, 1}, {1, 2}, {2, 0}, {2, 3}, {3, 4}, {4, 3}};
    for (size_t i = 0; i < 6; i++) graph_add_edge(links, arrow[i][0], arrow[i][1], 1.0);
    int scc[6];
    size_t n_scc;
    graph_strongly_connected_components(links, scc, &n_scc);
    printf("\n\nstrongly connected components: %zu   (vertex -> component:", n_scc);
    for (int v = 0; v < 6; v++) printf(" %d->%d", v, scc[v]);
    printf(")");

    /* ---- Maximum flow and the minimum cut that limits it ---- */
    Graph *pipes = graph_create_sized(4, true);              /* 0 = source, 3 = sink; weights are capacities */
    const struct { int from, to; double cap; } pipe[] = {{0, 1, 10}, {0, 2, 5}, {1, 2, 15}, {1, 3, 5}, {2, 3, 10}};
    for (size_t i = 0; i < 5; i++) graph_add_edge(pipes, pipe[i].from, pipe[i].to, pipe[i].cap);
    double flow;
    graph_max_flow(pipes, 0, 3, GRAPH_FLOW_DINIC, &flow, NULL);
    GraphEdge cut[5];
    size_t n_cut;
    double cut_value;
    graph_min_cut(pipes, 0, 3, &cut_value, cut, &n_cut);
    printf("\n\nmax flow 0 -> 3: %.0f; the cut that limits it (%.0f) is made of", flow, cut_value);
    for (size_t i = 0; i < n_cut; i++) printf(" %d->%d", cut[i].u, cut[i].v);

    /* ---- Colouring: schedule exams so that no student has two at once ---- */
    Graph *clash = graph_create_sized(5, false);             /* an edge = a student takes both */
    const int share[][2] = {{0, 1}, {0, 2}, {1, 2}, {2, 3}, {3, 4}};
    for (size_t i = 0; i < 5; i++) graph_add_edge(clash, share[i][0], share[i][1], 1.0);
    int slot[5];
    size_t n_slots;
    graph_color_dsatur(clash, slot, &n_slots);
    printf("\n\nexam time slots needed: %zu   (exam -> slot:", n_slots);
    for (int e = 0; e < 5; e++) printf(" %d->%d", e, slot[e]);
    printf(")\n");

    graph_destroy(map); graph_destroy(plan); graph_destroy(links); graph_destroy(pipes); graph_destroy(clash);
    return 0;
}
