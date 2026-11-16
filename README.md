# Graph Algorithms Library (C)

A graph library in C with **one unified API** over two storage representations, an **adjacency list**
and an **adjacency matrix**. You use the graph; the library decides how it is stored, keeps it inside a
memory budget, and gives you about 90 algorithms, a set of advanced ones, per-user social queries,
streaming analysis of files too big to load, and a visualiser that draws graphs of any size.

![How a call travels through the library](docs/images/architecture.svg)

**Contents**
1. [Quick start](#1-quick-start)
2. [How the library works](#2-how-the-library-works)
3. [Everyday algorithms](#3-everyday-algorithms)
4. [Social-network queries](#4-social-network-queries)
5. [Big graphs and limited memory](#5-big-graphs-and-limited-memory)
6. [Advanced algorithms](#6-advanced-algorithms)
7. [Drawing graphs](#7-drawing-graphs)
8. [Choosing a representation](#8-choosing-a-representation)
9. [Performance and test results](#9-performance-and-test-results)
10. [API reference](#10-api-reference)
11. [Repository layout](#11-repository-layout-and-branches) · 12. [Build and test](#12-build-and-test) ·
13. [Conventions](#13-conventions) · 14. [Contributing](#14-contributing-what-goes-on-which-branch-and-how-to-add-a-function) ·
15. [Status](#15-status)

---

## 1. Quick start

**Build** (any C11 compiler; the only dependency is the math library):

```
make                                   # builds build/libgraph.a and every test
cc -Iinclude my_program.c build/libgraph.a -lm -o my_program
make examples                          # builds and runs every program in this README
```

**Your first program.** `examples/quickstart.c` builds a small friendship network and asks it questions.
Vertices are plain integers; the program keeps its own array of names.

```c
#include "graph_algorithms.h"   /* the Graph handle plus every algorithm */
#include "graph_social.h"       /* per-person queries */

Graph *g = graph_create_sized(10, false);          /* 10 people, undirected; starts as an adjacency list */
graph_add_edge(g, ANA, BEN, 1.0);                  /* a friendship (the weight is 1.0)                    */
/* ... nine more friendships: two groups of four, a bridge Dev-Eli, Ivo knows only Hana, Jo knows nobody ... */

int component[10];
size_t groups = graph_connected_components(g, component);      /* who can reach whom */

GraphPath chain;                                               /* how are Ana and Ivo connected? */
if (graph_social_path(g, ANA, IVO, SIZE_MAX, &chain) == GRAPH_OK) {
    /* chain.vertices[0 .. chain.length), chain.cost = number of hops */
    graph_path_free(&chain);                                   /* results that own memory have a *_free */
}

double rank[10];
graph_pagerank(g, 0.85, 1e-9, 100, rank, NULL);                /* who is most influential? */

int community[10]; double modularity; size_t ncomm;
graph_louvain(g, community, &modularity, &ncomm);              /* which groups exist? */

GraphSuggestion suggestion[3]; size_t count;
graph_social_suggest_friends(g, BEN, 3, suggestion, &count);   /* whom should Ben meet? */

if (!graph_add_edge(g, ANA, 99, 1.0))                          /* errors are values, never crashes */
    printf("failed: %s\n", graph_status_string(graph_last_error()));
graph_destroy(g);
```

**Output** of the full program (`./build/examples/quickstart`):

```
10 people, 10 friendships, stored as a list

connected groups: 2
  group 0: Ana Ben Cara Dev Eli Fay Gus Hana Ivo
  group 1: Jo

shortest chain Ana -> Ivo (6 hops): Ana Cara Dev Eli Gus Hana Ivo
chain Ana -> Jo: no such path or solution

most influential (PageRank): Gus (0.144)

communities found: 4 (modularity 0.415)
  community 0: Ana Ben Cara Dev
  community 1: Eli Fay Gus
  community 2: Hana Ivo
  community 3: Jo

people Ben may know:
  Dev (1 mutual friend)

add_edge(Ana, 99) failed: invalid or removed vertex
```

![The ten people, coloured by the communities the library found](docs/images/friends.svg)

---

## 2. How the library works

1. **Vertices are integers.** The core never sees names or user structs, which is what lets a list and a
   matrix behave identically. IDs live in `[0, graph_id_bound())`; a removed vertex leaves a dead ID that
   `graph_add_vertex` reuses. IDs never change on their own, not even when the graph changes representation;
   only `graph_compact()` renumbers, and it reports old -> new IDs. (`mapping.h` maps your own keys to IDs.)
2. **A graph is stored in one representation at a time**, behind a function table (`GraphBackend`,
   `include/backend.h`). Both backends implement the same operations and report their own cost model.
   *List*: a sorted neighbour array per vertex, about 47 bytes per edge, any size. *Matrix*: bit rows plus
   sparse 64x64 weight tiles, very fast `has_edge` and row scans, but it needs V^2/8 bytes of presence bits.
3. **Adaptive mode (default).** An *observer* counts the operations you perform, predicts what that workload
   would cost on each representation, and migrates the graph when the gain is clear, the saving repays the
   rebuild, and the target fits the memory budget and the matrix's size limit. Migration only ever happens at
   the end of a call, so no pointer you were given can dangle. Turn it off with `graph_set_adaptive(g, false)`.
4. **Memory is bounded.** Every graph may use at most `graph_memory_budget()` bytes (4 GiB by default). A
   request that would exceed it (adding an edge, converting, loading, creating a too-large matrix) fails with
   `GRAPH_ERR_TOO_LARGE` and leaves the graph untouched, instead of the operating system killing the program.
5. **Errors are values.** Functions that return a `GraphStatus` hand it back; the rest (`bool`, `size_t`,
   `double`) record the reason for `graph_last_error()`. Nothing in the library calls `exit()`.
6. **Results that own memory** (`GraphPath`, `GraphGroups`, `GraphAdvCycles`) have an explicit `*_free`.
   Whoever creates an object destroys it; the graph never owns your data.

```
 your data <-> Mapping <-> int IDs <-> Graph (graph.h) --> Observer (policy.h)
 ("Paris")    (mapping.h)    0..n-1        |                    | migrate.h
                                            v                    v
                                    GraphBackend vtable:   list  |  matrix
```

---

## 3. Everyday algorithms

`examples/algorithms.c` runs the classics on graphs small enough to check by hand: shortest paths and a
spanning tree on a road map, task ordering on a dependency graph (and what happens when it has a cycle),
strongly connected components, maximum flow with its minimum cut, and exam scheduling by colouring.

```c
double dist[7]; int parent[7];
graph_dijkstra(map, 0, dist, parent);                          /* shortest distance from A to every town  */
size_t hops = graph_path_from_parent(parent, 7, 0, 6, route, 7);   /* rebuild the route A -> G           */

GraphEdge tree[7]; size_t n_tree;
double km = graph_mst_kruskal(map, tree, &n_tree);             /* cheapest roads keeping all towns linked */

int order[5];
if (graph_topological_sort(plan, order)) { /* order[] respects every "must happen before" */ }

int scc[6]; size_t n_scc;
graph_strongly_connected_components(links, scc, &n_scc);

double flow;  graph_max_flow(pipes, 0, 3, GRAPH_FLOW_DINIC, &flow, NULL);
double cut_value; GraphEdge cut[5]; size_t n_cut;
graph_min_cut(pipes, 0, 3, &cut_value, cut, &n_cut);           /* the edges that limit the flow */

int slot[5]; size_t n_slots;
graph_color_dsatur(clash, slot, &n_slots);                     /* fewest time slots with no clashes */
```

Output (`./build/examples/algorithms`; towns are A..G = 0..6):

```
Dijkstra from A:  A=0  B=7  C=15  D=5  E=14  F=11  G=22
  route A -> G: A D F G

minimum spanning tree: 39 km using 6 roads: A-D C-E D-F A-B B-E E-G

an order that respects every dependency: design build docs test release
after adding release -> design: topological sort fails (graph contains a cycle)

strongly connected components: 3   (vertex -> component: 0->1 1->1 2->1 3->0 4->0 5->2)

max flow 0 -> 3: 15; the cut that limits it (15) is made of 0->1 0->2

exam time slots needed: 3   (exam -> slot: 0->1 1->2 2->0 3->1 4->0)
```

How to read it: from A the shortest road to G is 22 km (A-D-F-G); the cheapest network that keeps every town
connected costs 39 km; "design, build, docs, test, release" is one valid order, and adding *release -> design*
creates a cycle, which the library reports as an error instead of returning nonsense; the 0->3 network can
carry at most 15 units because the two pipes leaving the source (10 + 5) are the bottleneck.

---

## 4. Social-network queries

`graph_social.h` answers the questions an application asks about **one person** or **one pair**. Their cost
depends on the neighbourhood they touch, not on the size of the network, so they stay fast at a million users
(40 query sets on a 1,000,000-user network take about 4.5 s including the checking).

| Function | What it does |
|---|---|
| `graph_social_within_hops(g, user, k, cap, users, hops, &count, &truncated)` | everyone within `k` hops, nearest first |
| `graph_social_suggest_friends(g, user, k, out, &count)` | "people you may know": ranked by number of mutual friends |
| `graph_social_path(g, a, b, max_hops, &path)` | shortest chain of friends, searching from both ends |

Users with no connections are fine: they simply have no neighbourhood (`path` returns
`GRAPH_ERR_NOT_FOUND`, `within_hops` returns zero people). The quick start above uses all three; Jo, who knows
nobody, gets "no such path or solution". `examples/social_network.c` (`make social`) builds a growing society of
up to 20,000 people and analyses it end to end.

---

## 5. Big graphs and limited memory

The list backend stores a million users with ten friends each in about **226 MB** (47 bytes per friendship).
The library gives you five tools for going further, in the order you would reach for them
(`examples/big_graph.c`):

1. **The memory budget.** `graph_memory_budget()` (4 GiB default; `graph_set_memory_budget(bytes)` changes it for
   graphs created afterwards, 0 = unlimited). Over-budget requests fail with `GRAPH_ERR_TOO_LARGE`.
2. **Unweighted graphs.** `GraphConfig.unweighted = true` stores no weights at all: every edge reads back as 1.0
   and any other weight is refused. The list drops from 47 to **23 bytes per edge**; the matrix needs no 32 KB
   weight tiles.
3. **Shrink-to-fit.** `graph_shrink_to_fit(g)` returns the growth slack the adjacency lists keep (17-24% of
   the graph) once a graph is built. Together with (2), the same 4-million-edge graph goes from 159 MB to 61 MB.
4. **Streaming.** `graph_stream_*` analyse an edge-list *file* in one or two passes without building the graph,
   keeping only a few integers per user, so the file may be far larger than memory: user and friendship counts,
   degrees, and connected components.
5. **Sharding (opt-in, for small devices).** Lower the budget and the library keeps the friendships on disk and
   processes them in pieces (see "Tiny devices" below).

```c
printf("memory budget: %zu MiB\n", graph_memory_budget() >> 20);
Graph *g = graph_create_sized(200000, false);                  /* 200 000 users ... */
if (!graph_convert(g, GRAPH_MATRIX))                           /* ... far too many for a matrix */
    printf("refused: %s\n", graph_status_string(graph_last_error()));

GraphConfig cfg = graph_default_config(false);
cfg.unweighted = true;
Graph *lean = graph_create_with(&cfg, 200000);

FILE *file = tmpfile();  graph_write_edge_list(g, file);       /* analyse the FILE, not a graph */
rewind(file);  GraphStreamInfo info;  graph_stream_scan(file, &info);
rewind(file);  graph_stream_components(file, info.vertices, component, &groups);
```

Output (`./build/examples/big_graph`):

```
memory budget: 4096 MiB
200000 users, 799987 friendships, stored as a list
convert to matrix refused: input too large
still a list, still usable: user 0 has 12 friends
largest matrix the budget allows: 131072 vertices weighted, 175851 unweighted

unweighted copy: 799987 friendships
add_edge with weight 2.5: invalid weight (not finite, or not 1.0 in an unweighted graph)

streamed file: 200000 users, 799987 friendships
connected groups found while streaming: 745
users with no friends: 744, most friends of any user: 20
```

**The matrix is limited by a calculated size.** Once every row strip exists, the presence bits take V^2/8
bytes, so `graph_matrix_max_vertices(weighted)` is `sqrt(8 * share * budget)`: **131,072 vertices weighted and
175,851 unweighted** under 4 GiB (above, a matrix is refused however few edges it starts with, and adaptive
mode never offers it). Measured on a social graph with about four friends per user:

| Users | List (unweighted) | Unweighted matrix | Weighted matrix |
|---|---|---|---|
| 50,000 | 3 MB | 298 MB | refused (a 32 KB tile per edge) |
| 100,000 | 7 MB | 1.19 GB | refused |
| 150,000 | 11 MB | 2.7 GB | refused |
| 200,000 | 14 MB | refused (over the limit) | refused |

For sparse graphs such as social networks the list is the right choice at every size; the matrix pays off for
dense or small graphs, where `has_edge` and row scans are a single bit test (see `make crossover`).

### Tiny devices: the sharding engine (`graph_shard.h`, `examples/sharded.c`)

A device with little memory (a Raspberry Pi, say) cannot hold a big graph, and the matrix cannot help (a million
users would need ~125 GB of presence bits). So instead of squeezing the graph into RAM, you **lower the memory
limit and switch the sharding engine on**: the friendships stay in *shard files* on disk, each holding the
adjacency lists of a range of user IDs, and an algorithm becomes a **sweep**: load a shard, process its rows, free
it, move on. Only the per-user numbers (4-24 bytes per user: degrees, PageRank values, BFS levels) stay in memory,
so a shard can name any neighbour without that neighbour's shard being loaded. Shards are read-only, so evicting
one is just freeing it, and the only "lock" is a pin that keeps the shard in use from being evicted.

```c
graph_set_memory_budget(10u << 20);                         /* the limit: 10 MiB */
ShardConfig cfg = shard_default_config();                   /* budget = the one just set; cfg.force = true
                                                               would shard even a graph that fits */
Graph *g = NULL;  ShardedGraph *sg = NULL;
graph_open_edge_list(file, &cfg, &g, &sg);                  /* exactly one of g / sg is returned */
if (sg) {
    shard_components(sg, component, &groups);               /* each of these is a sweep over the shards */
    shard_pagerank(sg, 0.85, 1e-9, 100, rank, &rounds);
    shard_bfs(sg, user, 2, level, &reached);                /* only shards holding the frontier are read */
    shard_destroy(sg);                                      /* deletes the shard files */
}
```

Output (`./build/examples/sharded`; times vary):

```
budget 10 MiB -> sharded (the graph would not fit in memory)
200000 users, 692625 friendships in 13 shard files; the cache holds 7 at a time

connected groups: 40454  (0.04 s)
most influential user: 6 (PageRank 0.00067, 1145 friends), 35 rounds  (0.47 s)
people within 2 hops of them: 26818  (0.00 s)
disk traffic: 243 shard loads, 127 MB read

budget 4 GiB  -> loaded in memory
in-memory PageRank: 35 rounds (2.31 s); largest difference from the sharded result: 1.1e-18
```

What to know:

- **Same answers.** The engine merges repeated edges exactly as the in-memory loader does, so sharded and in-memory
  results agree (here to 1e-18). It ignores edge weights (every edge counts as 1).
- **Speed.** Every algorithm pass reads the shards from disk, so on a device with a slow SD card the sharded mode
  is slower, which is the price of fitting. On the machine above the disk reads are served from the operating
  system's file cache, so it is not slower here; do not read the timings as a statement about a Pi.
  PageRank alternates the direction of its sweep each round, which cut disk reads on a million-user run from
  331 shard loads (1,010 MB) to 137 (416 MB).
- **Limits.** What fits is set by the per-user arrays: about 24 bytes per user, so 64 MB allows roughly 1.4
  million users, and a million users with five million friendships was processed in 200 MB of address space (the
  test caps it). The engine offers sweeps (degrees, components, PageRank, hop-limited BFS); algorithms that
  jump unpredictably between users (Dijkstra, A*, max flow) need the graph in memory.
- **Mode switch.** `graph_open_edge_list` loads an ordinary `Graph` when the graph fits the budget and shards it
  when it does not, so the same program runs on a big machine and a small one.

---

## 6. Advanced algorithms

`graph_advanced.h` (`graph_adv_*`) adds algorithms beyond the everyday set. It is generated from
`list_advanced.h` by `tools/gen_graph_advanced.py`, so it always matches the backends. `examples/advanced.c`:

```c
graph_adv_widest_path(roads, 0, width, NULL);                         /* heaviest load that fits every road */
graph_adv_tree_traversal(tree, 0, GRAPH_ADV_TREE_POSTORDER, order, &n);  /* pre / in / post / level order  */
graph_adv_graph_isomorphism(a, b, false, &same, NULL);                /* same shape, different numbering?   */
graph_adv_subgraph_isomorphism(triangle, roads, true, false, NULL, &found);   /* contains this pattern?    */
graph_adv_johnson(roads, dist, pred);                                 /* all pairs; pred = vertex before v   */
graph_adv_apsp_path(dist, pred, 5, 0, 4, &route);                     /* rebuild one route from `pred`       */
GraphAdvDynConn *dc = graph_adv_dynconn_create(roads);                /* "are u, v connected?" as it changes  */
```

Output (`./build/examples/advanced`):

```
heaviest load from 0 to 4: 3 t
pre-order  : 0 1 3 4 2 5 6
post-order : 3 4 1 5 6 2 0
level-order: 0 1 2 3 4 5 6
5-cycle vs relabelled 5-cycle isomorphic: yes
roads contain a triangle (0-1-2): yes
shortest road 0 -> 4: 11 km via 0 2 1 3 4

connected 0-4 now: yes
connected 0-4 after closing road 3-4: no
```

The full list: longest path, reverse-delete MST, Rocha-Thatte cycle enumeration, tree traversals, zero-weight
cycle detection, Coffman-Graham scheduling, crossing count and heuristic, colour refinement, colour-coding path,
harmonic and eigenvector centrality, flood, force and circle layouts, a GNN forward pass, IDA*, iterative-deepening
DFS, jump point search, widest path, Lex-BFS, isomorphism and subgraph isomorphism, Johnson, and a dynamic
connectivity tracker. PageRank, Floyd-Warshall, k-shortest paths and A* also exist in `matrix_advanced.h` /
`list_advanced.h` but are not exposed here, because the everyday API already has them at the same speed.
`graph_adv_graph_isomorphism` is about 50x faster than `graph_is_isomorphic` at 400 vertices and also matches
edge weights if asked.

---

## 7. Drawing graphs

`graph_visualize(g, "out.png", NULL)` draws any graph (SVG for small, PNG for large, from a built-in encoder with
no libraries). `examples/readme_images.c` (`make images`) produced the pictures below with the library itself.

```c
VizOptions o = graph_viz_default_options();
o.color  = VIZ_COLOR_COMMUNITY;        /* or _DEGREE, _COMPONENT, _NONE */
o.title  = "2400 people with friends, coloured by community (Louvain)";
graph_visualize(net, "social_communities.png", &o);

o.layout = VIZ_LAYOUT_MATRIX;          /* the adjacency matrix as a density heat-map */
graph_visualize(sorted, "social_matrix.png", &o);
```

**By community.** Louvain finds the eight planted groups; each gets a colour, a name tag, and a legend entry
with its size and most-connected member.

![Social network coloured by community](docs/images/social_communities.png)

**By popularity.** A temperature scale: blue for few friends, dark red for very popular, on a log scale, with the
most-connected people named.

![Social network coloured by degree](docs/images/social_degree.png)

**As a matrix.** With people renumbered by community, each group is a block on the diagonal, and the sparse
specks off the diagonal are the friendships that cross communities. This view stays readable when a node-link
drawing would be a solid blob.

![The same network as an adjacency matrix](docs/images/social_matrix.png)

**Structure survives the layout.** A 40 x 40 lattice drawn by the force-directed layout:

![A 40 by 40 lattice](docs/images/grid.png)

### Gallery: what each kind of graph looks like (`examples/gallery.c`, `make gallery`)

One picture per kind and size of graph, drawn by the library's own visualiser (`make gallery` writes them to
`docs/images/gallery/` in about half a minute and prints a line for each). Small graphs are SVG with every vertex
labelled; larger ones are PNG with the most-connected users named.

**Small graphs.** Path, cycle, star, complete graph, binary tree, grid, Petersen graph, complete bipartite graph,
a layered DAG and a ring of cliques:

| | | |
|---|---|---|
| ![path](docs/images/gallery/path_12.svg) | ![cycle](docs/images/gallery/cycle_12.svg) | ![star](docs/images/gallery/star_12.svg) |
| ![complete graph K8](docs/images/gallery/complete_8.svg) | ![binary tree](docs/images/gallery/binary_tree_31.svg) | ![grid 6x6](docs/images/gallery/grid_6x6.svg) |
| ![Petersen graph](docs/images/gallery/petersen_10.svg) | ![complete bipartite K3,4](docs/images/gallery/bipartite_3_4.svg) | ![layered DAG](docs/images/gallery/dag_4x4.svg) |

**Medium graphs.** A random tree, a random graph (a hairball: no structure to find), a small world (a ring of
neighbours with a few long shortcuts) and twenty groups of friends in a ring (Louvain recovers all twenty):

| | |
|---|---|
| ![random tree](docs/images/gallery/random_tree_500.png) | ![random graph](docs/images/gallery/random_graph_400.png) |
| ![small world](docs/images/gallery/small_world_500.png) | ![ring of cliques](docs/images/gallery/ring_of_cliques_20x6.png) |

**One kind of graph, growing.** Preferential attachment (each newcomer befriends the already popular) at 100,
1,000, 10,000 and 50,000 users: a handful of hubs appear and the picture turns from a readable network into a
dense core with a fringe:

| | |
|---|---|
| ![hubs, 100 users](docs/images/gallery/hubs_100.png) | ![hubs, 1000 users](docs/images/gallery/hubs_1000.png) |
| ![hubs, 10000 users](docs/images/gallery/hubs_10000.png) | ![hubs, 50000 users](docs/images/gallery/hubs_50000.png) |

**Large graphs.** A 100 x 100 grid, and a social network of 100,000 users (only the 80,000 who have a friend are
drawn; the silent fifth would be a ring of isolated dots):

| | |
|---|---|
| ![grid 100x100](docs/images/gallery/grid_100x100.png) | ![social network of 100000 users](docs/images/gallery/social_100000.png) |

**Adjacency matrices: the type of graph is visible in the picture.** A random graph is evenly scattered dots, a
hub graph has dense stripes along the popular users, a ring of cliques is solid blocks on the diagonal, a grid is
thin bands, and a social network sorted by community shows its communities as blocks:

| | | |
|---|---|---|
| ![random graph matrix](docs/images/gallery/matrix_random.png) | ![hub graph matrix](docs/images/gallery/matrix_hubs.png) | ![ring of cliques matrix](docs/images/gallery/matrix_cliques.png) |
| ![grid matrix](docs/images/gallery/matrix_grid.png) | ![social network matrix](docs/images/gallery/matrix_social.png) | |

---

## 8. Choosing a representation

| You want | Do |
|---|---|
| Don't care | `graph_create(directed)`. Adaptive, starts as a list. |
| Start on a particular one | `GraphConfig c = graph_default_config(d); c.start = GRAPH_MATRIX; graph_create_with(&c, n);` (refused if `n` is above `graph_matrix_max_vertices`) |
| Never convert automatically | `c.adaptive = false`, or `graph_set_adaptive(g, false)` at any time. |
| Convert by hand | `graph_convert(g, GRAPH_MATRIX)`; fails with `GRAPH_ERR_TOO_LARGE` if it would not fit. In adaptive mode the observer may move the graph again later; turn adaptation off first if the choice should stick. |
| Half the memory | `c.unweighted = true` (every edge is 1.0; other weights are refused). |
| Know what is in use | `graph_representation(g)` (inspection only; never needed to operate). |

Adaptive tuning (check interval, hysteresis, cooldown, memory budget, log stream) is `GraphConfig.observer`;
the cost models live at the bottom of `list_backend.c` and `matrix_backend.c`.

### Comparing the two (`graph_compare`)

```c
GraphCompareResult r;
graph_compare(g, GRAPH_ALGO_DIJKSTRA, /*source*/ 0, /*repetitions*/ 5, &r);
graph_compare_print(&r, stdout);
/*  Algorithm : Dijkstra (source 0)
      list  :      2516127 ns       749752 bytes
      matrix:      6802731 ns     34062968 bytes
      faster: list (2.70x)   results: equivalent  */
```

It copies the graph into both representations, times the algorithm on each
(best of N), and checks the answers are equivalent. Your graph, its
statistics and its representation are untouched. `graph_compare_all` runs every
algorithm. "Equivalent" allows different neighbour orders (BFS/DFS compare the
set of vertices reached, topological sort compares validity).


---

## 9. Performance and test results

**One million users.** The scale suite (`test/test_scale.c`, `make build/scale_test`) runs each algorithm on a
generated network of up to 1,000,000 vertices, checks the answer against an independent implementation or a
known closed form (never against the library itself), caps each case at 4 GiB of address space, and records its
peak memory. On the list backend all 44 cases at 1,000,000 vertices pass; the slowest are shown, the highest
peak memory of any case was 986 MB:

![Time and peak memory per algorithm at one million vertices](docs/images/scale_1M.svg)

Some algorithms are deliberately capped at smaller sizes because they are not linear (all-pairs paths, closeness,
transitive closure, exact search, Lex-BFS, colour refinement on paths, isomorphism). Over the whole run:
234 of 234 list cases (1,000 to 1,000,000 vertices) and 101 of 101 matrix cases (1,000 and 10,000 vertices) pass.


### A million users: the benchmark (`examples/million_users.c`, `make bench-million`)

One program builds a social network of **1,000,000 users and 3,103,057 friendships** (communities, hubs, and a
fifth of the users with no friends), runs every whole-graph analysis on it, checks each answer with independent code,
measures the per-user queries, holds the graph three ways, runs the low-memory sharded mode, draws 17 pictures, and
writes a report (`docs/images/million/REPORT.md`, `results.csv`). **This is the number to compare against.**

| | |
|---|---|
| **Suite score** (geometric mean of the analyses' best-run rates) | **6.3 million friendships per second** |
| Analysis suite time (7 scored analyses, best runs added up) | 14.96 s |
| Independent checks passed | 8 of 8 |
| Peak resident memory after the analyses | 295 MB |
| Hardware | i7-1165G7 laptop, 31 GB, one core used (the library is single-threaded), built with `-O2` |

| Analysis | best (s) | median (s) | worst (s) | M friendships/s |
|---|---:|---:|---:|---:|
| degree histogram | **0.004** | 0.005 | 0.005 | 716.7 |
| connected components | **0.122** | 0.124 | 0.136 | 25.5 |
| k-core decomposition | **0.489** | 0.573 | 0.813 | 6.4 |
| triangles + clustering | **2.985** | 3.045 | 3.767 | 1.0 |
| PageRank (to 1e-8) | **7.726** | 11.310 | 13.516 | 0.4 |
| Louvain communities | **3.295** | 3.565 | 4.160 | 0.9 |
| hop distances from the top hub | **0.341** | 0.352 | 0.364 | 9.1 |

Per-user queries: "people you may know" takes about 0.01 ms (roughly 60,000 per second) and a shortest chain of
friends about 0.25 ms. Holding the network costs 131 MB weighted, 64 MB unweighted and 54 MB shrunk (18 bytes per
friendship); a dense matrix would need 116 GB of presence bits and is refused. On a 64 MiB budget the sharding engine
cut the friendships into 9 shard files, kept 28 MB resident, and ran PageRank in about 1.6 s (110 shard loads; earlier runs took up to 2.7 s on the loaded machine).

| | | |
|---|---|---|
| ![degrees](docs/images/million/01_degree_distribution.png) | ![PageRank](docs/images/million/05_pagerank_top20.png) | ![hops](docs/images/million/09_hops.png) |
| ![memory](docs/images/million/10_memory.png) | ![timings](docs/images/million/11_timings.png) | ![communities](docs/images/million/16_community_graph.png) |
| ![backbone](docs/images/million/13_backbone_communities.png) | ![matrix](docs/images/million/15_backbone_matrix.png) | ![ego network](docs/images/million/12_ego_network.png) |

The whole network (the 800,000 users who have friends) is `docs/images/million/17_whole_network.png`; the other charts
(connected groups, k-core, clustering, PageRank against friends, community sizes) are in the same folder, each
described in `REPORT.md`.

**The result cache, measured.** A second copy of the network has the cache on (`cache_entries = 32`, `cache_bytes` raised
to 256 MiB, because the default 4 MiB cannot hold a million-entry answer). For each of three cached calls the report gives
the first call, the median of five repeats and the speed-up, read with `graph_cache_stats()` (hits, misses, stores,
invalidations). On the last run a repeat was 56 to 110 times faster (components 127 ms to 1.3 ms, BFS 140 ms to 2.5 ms,
Dijkstra 376 ms to 3.4 ms). A repeat is a copy of the stored answer, so it is not free. Each repeat is compared byte for
byte with the first answer, and after one new friendship the next call must recompute and find one group fewer, so a stale
answer would fail the run. `graph_mst_kruskal` is not in the table: this network is disconnected, so it has no spanning tree
and nothing is cached.

**How to run and compare.** `make bench-million` rebuilds the library with `-O2` in `build-o2/` and runs the example
(about three minutes); `make bench-million BASELINE=old/results.csv` adds a "vs baseline" column. Each analysis runs five
times, the score uses the best run, and the report records the machine's load average.

**Things the benchmark taught us, so you can trust the number:**

- **Build flags matter, but less than you would think.** The default `make` has no optimiser. Measured on this workload,
  `-O2` is only **1.3 to 1.7 times faster** than `-O0` (these algorithms are memory-bound), so never compare numbers
  across flags. Preparing the `-O2` build also found and fixed a real defect (`mu_mul` allowed allocation sizes above
  `PTRDIFF_MAX`, which GCC rejects at `-O2`) and an undefined-behaviour bug in a test; the whole suite now passes at `-O0`
  and `-O2`.
- **The result cache can fake a result.** `graph_connected_components` and the other core calls are served from the
  library's own cache on repeat calls (a 110 ms call becomes 1.3 ms; the cache's hit counters confirm it). That is a real
  feature, but it made a first version of this benchmark report an impossible 3,500 million friendships per second. The
  benchmark graph therefore has `cache_entries = 0`. The cache gets its own section in the report instead (below).
- **Noise is large on a shared laptop.** Five complete runs scored 3.7 to 7.1 million friendships per second because other
  programs were running and the CPU changes frequency (PageRank's five repeats within one run ranged from 7.7 to 13.5 s).
  That is why the score uses the best of five runs and the report records the load average; for a comparison, run on an
  idle machine.
- **The degree distribution is heavy-tailed, not a power law.** The generator caps how many links a newcomer makes, so the
  tail bends down faster than a power law; the chart says so.

**The verification report (`make verify`).** `test/verify_report.c` prints, for every area of the library, what the
answer *should* be next to what the library *gave*, and whether they match. Expected values never come from the
library: they are closed formulas for whole families of graphs (a path of n vertices has n-1 edges and n-1
bridges; a k x k grid has diameter 2(k-1); K_n has C(n,3) triangles), published textbook answers (the CLRS
max-flow network has flow 23; the Kruskal example from the literature weighs 39; the Petersen graph has a
Hamiltonian path but no Hamiltonian cycle), each re-checked by a separate script, and, for the two
representations, each other. Seven sections: graph families from 10 to **1,000,000 vertices**, textbook
examples, the advanced algorithms, social queries, edge-list files (streaming and sharding), the memory limits,
and list-versus-matrix agreement. All **424 checks pass in about six seconds**
(`make verify-quick` skips the million-vertex graphs). An excerpt of the real output:

```
Verification report: what the answer should be, what the library gave, and whether they match

=== 1. Graph families, by size (answers from closed formulas)
        check                                                      expected                           actual
--- grid 1000 x 1000 = 1000000 vertices (list)
  ok    vertices (s*s)                                             1000000                            1000000
  ok    edges (2s(s-1))                                            1998000                            1998000
  ok    bipartite (checkerboard)                                   yes                                yes
  ok    triangles                                                  0                                  0
  ok    every vertex has core number 2                             yes                                yes
  ok    largest matching (floor(n/2))                              500000                             500000
  ok    colours (DSatur is exact on bipartite graphs)              2                                  2
  ok    corner to corner distance (2(s-1))                         1998                               1998
  ok    max flow corner to corner (each corner has 2 edges)        2                                  2
  ...
=== 2. Textbook examples (published answers, re-checked by a separate script)
        check                                                      expected                           actual
  ok    Dijkstra from s: distances to t, y, x, z                   8, 5, 9, 7                         8, 5, 9, 7
  ok    maximum flow, Dinic                                        23                                 23
  ok    maximum flow, Edmonds-Karp                                 23                                 23
  ok    minimum spanning tree weight, Kruskal                      39                                 39
  ok      Hamiltonian path exists                                  ok                                 ok
  ok      Hamiltonian cycle does not                               no such path or solution           no such path or solution
  ...
--- sharding engine forced on, 512 KB budget: friendships live on disk
  ok    open with sharding forced                                  ok                                 ok
  ok      came back as shards (not an in-memory graph)             yes                                yes
  ok      shard files exist in the directory                       yes                                yes
  ok    users                                                      1000                               1000
  ok    friendships                                                2200                               2200
  ok    users with 4 / 5 friends                                   600 / 400                          600 / 400
  ok    connected groups                                           1                                  1
  ok    PageRank converges                                         ok                                 ok
  ok      ranks add up to                                          1                                  1
  ok      same ranks as the in-memory graph (largest difference)   < 1e-9                             4.3e-19
  ok    people within 2 hops of user 0 (same as in memory)         11                                 11
  ...
--- a lowered budget refuses cleanly instead of crashing
  ok    create 5 000 000 users under a 1 MB budget                 no                                 no
  ok      the reason                                               input too large                    input too large
  ok    adding friendships until the budget is reached             refused (too large)                refused after 9840 edges
  ok      the graph kept exactly the edges that were accepted      9840                               9840
--- shrink-to-fit returns slack and keeps the contents
  ...
424 checks, 0 failed
```

---

## 10. API reference

`include/graph.h` (core) and `include/graph_algorithms.h` (the algorithm library; its wrappers are
generated from `list_algorithms.h` by `tools/gen_graph_algorithms.py`).

| Area | Functions |
|---|---|
| Lifecycle | `graph_create`, `graph_create_sized`, `graph_create_with`, `graph_destroy`, `graph_copy` |
| Representation | `graph_set_adaptive`, `graph_is_adaptive`, `graph_convert`, `graph_representation`, `graph_migration_count` |
| Sharing | `graph_freeze`, `graph_thaw`, `graph_is_frozen` (read-only handle, safe to share between threads) |
| Queries | `graph_vertex_count`, `graph_id_bound`, `graph_edge_count`, `graph_is_directed`, `graph_is_alive`, `graph_degree`, `graph_has_edge`, `graph_neighbors` |
| Mutation | `graph_add_vertex`, `graph_remove_vertex`, `graph_add_edge`, `graph_remove_edge`, `graph_needs_compaction`, `graph_compact` |
| Baseline algorithms (cached) | `graph_bfs`, `graph_dfs`, `graph_connected_components`, `graph_has_cycle`, `graph_topological_sort`, `graph_dijkstra`, `graph_mst_kruskal` |
| Algorithm library (~90) | see "Algorithms" below |
| Advanced algorithms | `graph_adv_*` in `graph_advanced.h` (isomorphism, Johnson, IDA*, jump point search, dynamic connectivity, ...) |
| Social queries | `graph_social_within_hops`, `graph_social_suggest_friends`, `graph_social_path` (`graph_social.h`) |
| Streaming (beyond RAM) | `graph_stream_scan`, `graph_stream_degrees`, `graph_stream_components` (`graph_stream.h`) |
| Memory | `graph_memory_budget`, `graph_set_memory_budget`, `graph_matrix_max_vertices`, `graph_shrink_to_fit` (`memory_budget.h`) |
| Building graphs | generators `graph_gen_*(rep, ...)`, `graph_transpose`, `graph_union`, `graph_intersection`, `graph_complement`, `graph_induced_subgraph`, readers `graph_read_edge_list`, `graph_load_binary` |
| Benchmark | `graph_compare`, `graph_compare_all`, `graph_compare_print`, `graph_crossover_report` |
| Visualization | `graph_visualize` (SVG / PNG file), `graph_visualize_to`, `graph_layout` (`graph_viz.h`) |
| Errors | `graph_last_error`, `graph_clear_error`, `graph_status_string` |
| Names | `NamedGraph` (`mapping.h`): `named_add_edge`, `named_remove_vertex`, `named_remove_edge`, `named_compact` |

Array sizes: arrays indexed by vertex ID need `graph_id_bound()` entries; arrays listing vertices need
`graph_vertex_count()`. The seven baseline algorithms are cached per graph version and invalidated on
any change; the library algorithms are not cached. A normal handle is **not** thread-safe (even reads
update statistics and may migrate the graph); freeze it to share it.

### Visualization (`graph_viz.h`)

```c
graph_visualize(g, "out.png", NULL);          /* format from the extension; any size */

VizOptions o = graph_viz_default_options();
o.layout = VIZ_LAYOUT_MATRIX;                  /* adjacency-matrix heat-map instead */
graph_visualize(g, "matrix.png", &o);
```

* **Any size.** The layout is multilevel force-directed: the graph is coarsened by edge matching,
  the small graph is laid out, then each finer level is refined, and repulsion only looks at nearby grid
  cells. 10k vertices take ~2 s and 100k vertices ~25 s (`-O2`). Structure such as lattices and
  clusters survives; a dense scale-free core is still a hairball, which is a property of the graph.
* **Formats.** SVG (zoomable, labels, edge weights, arrowheads) for small graphs; for large ones a PNG
  from a built-in encoder (no libraries). `VIZ_FORMAT_AUTO` picks SVG up to 1500 vertices / 8000 edges.
  The PNG has no text and draws arrowheads only for directed graphs of <= 300 vertices.
* **Style scales with size.** Edges get fainter as there are more of them so overlaps read as density;
  vertex size follows degree; colour is by component, community (Louvain) or degree. With more than 12
  groups only the 11 largest get a colour.
* **Reading the picture.** Colour by group (`VIZ_COLOR_COMMUNITY` / `_COMPONENT`): each of the 11 largest
  groups gets a colour, a name tag at the centre of its members ("Community 3"; set `group_name`), and a
  legend entry with its size and most-connected member. Colour by degree (`VIZ_COLOR_DEGREE`) is a
  **temperature scale** from blue (few connections) to dark red (very popular) on a log scale, with a colour
  bar, tick labels and a "Most popular" list; popular vertices are drawn last, larger and with a glow. The
  most-connected people are named on the picture (`label_top`, names from `labels`). PNG text uses a built-in
  bitmap font (`src/viz_font.c`, generated from Liberation Mono, SIL OFL), so there are still no libraries.
  `legend = false` removes the panel and the community tags.
* **Matrix view** (`VIZ_LAYOUT_MATRIX`) bins the adjacency matrix into pixels, which stays readable when
  a node-link drawing would be a solid blob (block structure shows as diagonal squares).
* `graph_layout` returns just the coordinates in [0,1]^2 if you want to draw it yourself.
* The pictures `make demo` writes (about a minute at `-O2`); `clusters400.png` is
  eight dense groups joined by single bridges, `clusters2000_matrix.png` is the matrix view of a
  2000-vertex version, `grid10000.png` a 100x100 lattice, `scalefree100000.png` a 100,000-vertex
  scale-free graph.
* Try it: `make demo CFLAGS="-std=c11 -O2"` writes example images to `build/demo/`
  (`viz_demo <random|scalefree|grid|tree|clusters> <n> <out.png|out.svg> [matrix]`).

### Showcase: expected vs actual (`test/showcase.c`)

`make showcase` runs about 85 library functions plus the baseline algorithms, generators, I/O,
visualization and error codes on nine small fixed graphs and prints, for every function, the
**expected** outcome next to the **actual** one:

```
[PASS] graph_dijkstra(U, 0): distances
         expected: [0 2 3 6 8 9 12 inf]
         actual  : [0 2 3 6 8 9 12 inf]    (list = matrix)
[PASS] graph_color_dsatur(U)
         expected: proper colouring using exactly 3 colours (DSatur is optimal here)
         actual  : 3 colours    (list = matrix)
```

* The expected values were derived **independently** of the C code, by hand and with a separate
  brute-force Python script, never by calling the library.
* Every check runs twice, on a list-backed and a matrix-backed graph, and the output shows both when
  they differ. 221 checks currently pass; it exits non-zero if any fails.
* Where an answer is not unique (colourings, community splits, Hamiltonian order) the expected column
  states the property and the check verifies it.
* It also draws the test graph to `build/showcase_U_list.svg` and `build/showcase_U.png`.

### Example: a growing social network (`examples/social_network.c`)

`make social CFLAGS="-std=c11 -O2"` (about 20 s) simulates a society that grows from six friends to
20,000 people and studies it with the library; `./build/social_network 50000 7` changes the size and seed.

* **Growth** in waves, with preferential attachment, friends-of-friends and a few cross-community ties,
  on an adaptive graph (the report says whenever the library changes how it is stored).
* **Churn**: 2% of people leave (`graph_remove_vertex`), then `graph_compact` renumbers the rest.
* **Analysis**: components, degree distribution, triangles, clustering, assortativity, Louvain communities
  (checked against the planted ones), PageRank influencers, betweenness brokers, largest clique and a
  DSatur seating plan on the inner circle, k-cores, bridges, small-world separation and a concrete chain
  of friends.
* **Dynamics**: an information cascade from the most influential person versus an ordinary one.
* **Threads**: three analyses at once on a frozen graph (`graph_freeze`).
* **Compare**: `graph_compare` times list against matrix on this very network.
* **Pictures** in `build/social/`: the 60-person network with names (SVG), the network at 2000 people, the
  final network by community and by degree, and the adjacency matrix with people sorted by community
  (twelve planted communities show up as diagonal blocks).

### Error codes

`graph_status.h` defines `GraphStatus`: `GRAPH_OK`, `GRAPH_ERR_INVALID_ARGUMENT`, `_INVALID_VERTEX`,
`_DUPLICATE_EDGE`, `_EDGE_NOT_FOUND`, `_INVALID_WEIGHT`, `_NEGATIVE_WEIGHT`, `_CYCLE`, `_OUT_OF_MEMORY`,
`_UNSUPPORTED`, `_NEGATIVE_CYCLE`, `_TOO_LARGE`, `_NOT_FOUND`, `_NO_CONVERGENCE`, `_IO`,
`_NO_BACKEND_SUPPORT`, `_FROZEN`, plus `graph_status_string()`.

* Algorithm-library functions **return** a `GraphStatus` (`graph_pagerank`, `list_bellman_ford`, ...).
* Core operations and baseline algorithms keep their `bool` / `size_t` / `double` / `int` return types and
  **record** the reason in a per-thread slot. After a failure read `graph_last_error()`; reset it with
  `graph_clear_error()`. A success never clears it.

### Algorithms

`list_algorithms.h` and `matrix_algorithms.h` expose the same ~85 functions (also reachable per
representation as `list_*` / `matrix_*`, with the same contracts); the unified `graph_*` versions pick the
representation for you and feed the observer.

* **Analysis**: degree/closeness/betweenness centrality, triangles, clustering, bipartite, eccentricity,
  diameter/radius, transpose, common neighbours, Jaccard, network stats, degree histogram, BFS/DFS visitors.
* **Connectivity**: connected/strongly connected/DAG/tree/forest predicates, Tarjan SCC, articulation
  points, bridges, biconnected components, find-cycle, transitive closure and reduction.
* **Paths**: Bellman-Ford, Floyd-Warshall and Johnson (N <= 8192), A*, multi-source Dijkstra, 0-1 BFS,
  Yen k-shortest paths, path reconstruction.
* **Flow / matching**: Prim, Dinic, Edmonds-Karp, min cut, Hopcroft-Karp, Hungarian.
* **Combinatorial**: Euler trails, Hamiltonian (<= 24 vertices), greedy/Welsh-Powell/DSatur colouring,
  Bron-Kerbosch cliques, vertex cover approximation, independent set, DAG longest/critical path.
* **Network science**: PageRank, HITS, eigenvector centrality, k-core, label propagation, Louvain,
  modularity, assortativity.
* **Operations and generators**: copy, induced subgraph, union, intersection, complement, isomorphism;
  Erdos-Renyi, Barabasi-Albert, grid, path, cycle, star, complete, random DAG, random tree.
* **I/O**: Graphviz DOT, plain and named edge lists, binary save/load, GraphML/GML/JSON export.

The observer prices library algorithms in three classes (`stats.h`): `OP_ITERATIVE_RUN` (repeated edge
sweeps: PageRank, Louvain, k-core, SCC, Bellman-Ford), `OP_ALL_PAIRS_RUN` (a traversal per vertex) and
`OP_ADJACENCY_RUN` (many `has_edge` probes: triangles, cliques, flow). Rule of thumb from `make
crossover`: sweeps favour the list at every density; probe-heavy work favours the matrix once the graph
is dense.

---

---

## 11. Repository layout and branches

| Branch | Contains |
|---|---|
| `main` | Everything implemented in both branches, plus the generic code: facade, cache, migration, policy, mapping, wrappers, stub backend. |
| `Adjacency_list` | **List code only**, plus the shared interface files (`backend.h`, `graph_types.h`, `graph_status.*`, `graph_results.*`, `stats.*`). |
| `adjacency_matrix` | **Matrix code only**, plus the same shared interface files. |

`main` is made by merging the two. The representation branches deliberately lack the generic files, so a
naive merge deletes them from `main` (it happened once). After merging, run
`git diff --name-status <old main> HEAD | grep '^D'`; it must print nothing.

```
include/
  graph.h, graph_algorithms.h, graph_viz.h   the unified public API
  graph_advanced.h (generated), graph_social.h, graph_stream.h   advanced / per-user / streaming APIs
  memory_budget.h               the library-wide memory budget (default 4 GiB)
  graph_status.h, graph_results.h  shared error codes and result types
  mapping.h                     value <-> ID mapping and NamedGraph
  backend.h, stats.h, policy.h, cache.h, migrate.h      adaptive machinery
  adjacency_list.h, list_backend.h, list_algorithms.h   list representation + algorithms
  adjacency_matrix.h, matrix_backend.h, matrix_algorithms.h   matrix representation + algorithms
src/
  graph.c, graph_build.c, graph_algorithms.c + graph_advanced.c (generated), graph_compare.c, graph_crossover.c, graph_internal.h
  graph_social.c, graph_stream.c, memory_budget.c
  graph_viz.c (layout + SVG + raster), png_writer.c/.h (PNG encoder)
  adjacency_list.c, list_backend.c, list_{analysis,connectivity,paths,flow,combinatorial,network,ops,io}.c, list_util.h
  adjacency_matrix.c, matrix_backend.c, matrix_{analysis,...}.c, matrix_util.h
  policy.c, cache.c, migrate.c, stats.c, mapping.c, graph_status.c, graph_results.c
test/
  test_graph.c, test_facade.c, test_adaptive.c (+ stub_backend), test_adjacency_list.c, test_adjacency_matrix.c,
  test_matrix_backend.c, test_list_*.c, test_matrix_*.c (+ test_util_{list,matrix}.h), bench_list.c, bench_matrix.c
tools/
  gen_graph_algorithms.py       regenerates include/graph_algorithms.h and src/graph_algorithms.c
  gen_graph_advanced.py         regenerates include/graph_advanced.h and src/graph_advanced.c from list_advanced.h
  make_readme_diagrams.py       draws docs/images/architecture.svg and scale_1M.svg
  crossover_report.c            `make crossover`
  viz_demo.c                    `make demo`
  (test/showcase.c              `make showcase`: expected vs actual for every function)
  (examples/social_network.c    `make social`: a growing social network, analysed and drawn)
  examples/{quickstart,algorithms,big_graph,advanced}.c   the README's programs (`make examples`)
  examples/readme_images.c      draws docs/images/*.png and friends.svg (`make images`)
  test/test_scale.c             the scale suite: 48 cases up to 1 000 000 vertices (`scale_test`)
docs/images/                    the pictures used in this README
```

---

## 12. Build and test

```
make                 # build/libgraph.a, every test binary, build/crossover_report
make test            # build and run all suites (a couple of minutes)
make bench           # CSV timing sweeps per representation   (use CFLAGS="-std=c11 -O2")
make crossover       # where the matrix overtakes the list      (use CFLAGS="-std=c11 -O2")
make showcase        # expected vs actual for every algorithm on fixed graphs (list and matrix)
make bench-million   # THE BENCHMARK: a million users analysed and drawn, report in docs/images/million/ (about 3 minutes)
make verify          # the verification report: expected vs actual over 424 checks, graphs up to 1M vertices
make verify-quick    # the same without the million-vertex graphs
make gallery         # one picture per kind and size of graph, in docs/images/gallery/
make demo            # example pictures in build/demo/            (use CFLAGS="-std=c11 -O2")
make social          # social-network simulation, pictures in build/social/ (use CFLAGS="-std=c11 -O2")
make examples        # build and run the README's example programs
make images          # redraw docs/images/ (use CFLAGS="-std=c11 -O2")
make build/scale_test && ./build/scale_test 1000000        # scale suite, list backend up to 1M vertices
make clean
make CFLAGS="-std=c11 -Wall -Wextra -Werror -pedantic -g -O1"   # debug build for valgrind
valgrind --leak-check=full --error-exitcode=9 ./build/facade_test
```

Flags are `-std=c11 -Wall -Wextra -Werror -pedantic`; link with `-lm` (and `-pthread` for `facade_test`).
Builds must stay warning-free at `-O0` through `-O3`. `libasan` is not installed here; use valgrind.
To use the library: `cc -Iinclude my_program.c build/libgraph.a -lm`.
After changing `list_algorithms.h` re-run `python3 tools/gen_graph_algorithms.py`; after changing
`list_advanced.h` re-run `python3 tools/gen_graph_advanced.py`.

---

## 13. Conventions

- C11. Failures are reported through return values and `graph_last_error()`; the library never calls `exit()`.
- Whoever creates an object destroys it. The graph never owns user data; the mapping owns its copied keys.
- Every public function documents its parameters, return value and complexity.
- Never commit or push without being asked; pushes are fast-forward only (`git fetch` first, no force).
- `CLAUDE.local.md`, `.claude/`, `.idea/`, `TASKS.md` are local-only (`.git/info/exclude`).

## 14. Contributing: what goes on which branch, and how to add a function

### What may be committed where

| Branch | Allowed | Not allowed |
|---|---|---|
| `Adjacency_list` | `adjacency_list.*`, `list_*` headers, sources and tests, `test_util_list.h`, `bench_list.c`; the **shared interface files** (below); its own `Makefile` and `README.md` | any `matrix_*` / `adjacency_matrix.*` file; any `main`-only file (list below) |
| `adjacency_matrix` | the same for `adjacency_matrix.*` and `matrix_*` | any `list_*` / `adjacency_list.*` file; any `main`-only file |
| `main` | everything: the merged union, plus the **`main`-only files**: `graph.h`, `graph_algorithms.*` (generated), `graph_build.c`, `graph_compare.c`, `graph_crossover.c`, `graph_internal.h`, `graph_viz.*`, `png_writer.*`, `viz_*`, `cache.*`, `policy.*`, `migrate.*`, `mapping.*`, `test/stub_backend.*`, `test_graph.c`, `test_facade.c`, `test_adaptive.c`, `test_viz.c`, `showcase.c`, `tools/`, `examples/`, this README | **editing `list_*`, `matrix_*`, `adjacency_*.{c,h}` directly**: change them on their own branch and merge, so the branches never drift |

**Shared interface files** (identical on all three branches): `backend.h`, `graph_types.h`,
`graph_status.h/.c`, `graph_results.h/.c`, `stats.h/.c`. Change one on any branch, then copy the exact file
to the others (`git checkout <branch> -- <file>`, commit on each). `test_util_list.h` / `test_util_matrix.h`
are deliberately separate files so both can live on `main`.

General rules: never commit or push unless asked; commit messages end with the `Co-Authored-By` line;
push fast-forward only after `git fetch`, never force-push; `main` receives merges from the representation
branches and is never merged back into them. Each representation branch has its own README with the same
rules written for that branch.

Before committing, check the paths: `git diff --name-only --cached` (and `git diff --name-only
origin/<branch>...HEAD` for everything unpushed). Every path must be allowed on that branch.

### Adding a new algorithm (the library, both representations)

A function is part of the library only when **both** representations have it with the same signature and
behaviour, and `main` exposes it on `Graph`.

1. **On one representation branch** (say `Adjacency_list`, function `list_foo`):
   declare it in `include/list_algorithms.h` with a full contract (arguments, size of every output array,
   `GraphStatus` codes, complexity); implement it in the matching `src/list_*.c` module using
   `list_graph_next_neighbor`; test it in `test/test_list_<module>.c` against an independent oracle on random
   graphs (with removed vertices, directed and undirected) and test every error code; add a new test file to
   `TESTS` in that branch's `Makefile`. Verify: `make test`, warning-free at `-O0..-O3`, valgrind clean, and
   check that deliberately breaking the code makes a test fail. Commit on that branch.
2. **On the other branch** (`adjacency_matrix`): the same function as `matrix_foo` (usually the list version
   with `list_` -> `matrix_` and `List` -> `Matrix`), same tests, same verification, commit.
3. **On `main`:** `git merge --no-ff Adjacency_list`, then `git merge --no-ff adjacency_matrix`. If `README.md`
   conflicts keep main's (`git checkout --ours README.md && git add README.md`) and keep main's `Makefile`
   additions. Before committing the merge run `git diff --name-status <old main> HEAD | grep '^D'`; it must print
   nothing (the branches lack the generic files, so a merge can delete them; restore with `git checkout
   <old main> -- <file>`).
4. **Expose it on `Graph`:** run `python3 tools/gen_graph_algorithms.py` (it reads `list_algorithms.h`, so
   `graph_foo` and its declaration appear automatically). Add `foo` to the right workload class in the
   generator (`ITERATIVE`, `ALL_PAIRS`, `ADJACENCY`, `WEIGHTED`; default is a plain traversal) so the adaptive
   observer prices it. If it creates or combines graphs, write the wrapper by hand in `src/graph_build.c` and
   add it to `MANUAL` in the generator.
5. **Show and cross-check it:** add a row to `test/showcase.c` whose expected value comes from hand work or an
   independent script (never from the function), and add it to the list-vs-matrix cross-check in
   `test/test_facade.c`.
6. **Finish:** `make test`, update this README's algorithm list and `TASKS.md`, commit on `main`, then push
   `Adjacency_list`, `adjacency_matrix` and `main`.

### Other kinds of change

* **A facade or generic feature** (cache, policy, visualization, examples, tools): write it on `main` only.
* **A new error code, result type or workload kind:** edit the shared file, append new codes before
  `GRAPH_ERR_COUNT` (and add the name in `graph_status.c`), then copy the file to all branches.
* **A new core operation on a representation:** add it to `adjacency_<rep>.*` and `<rep>_backend.c` on both
  branches, add the slot to the shared `backend.h`, then wire it into `graph.c` on `main`.
* **A change to the generator or Makefile on `main`** does not go to the branches; each branch keeps its own
  `Makefile`.

## 15. Status

**Complete**
- Unified API, adaptive migration (opt in/out), manual conversion, result cache, `graph_compare`,
  `graph_crossover_report`.
- The full algorithm library on both representations (about 90 functions) plus the advanced set
  (`graph_adv_*`), with the same API and the shared error codes.
- Per-user social queries, streaming analyses of edge-list files, visualization of any size.
- Memory safety: a library-wide budget (4 GiB default) enforced by both backends, an unweighted mode,
  shrink-to-fit, a sparse matrix tile index, and a calculated size limit for the matrix.
- Tests: oracle-based suites for every module, list-vs-matrix cross-checks, error-code tests, a multi-thread
  test, planted-bug checks, a scale suite to 1 000 000 vertices; valgrind-clean; warning-free at -O0..-O3.

**Known limits (deliberate)**
- GraphML/GML/JSON are export-only (the binary format round-trips).
- Isomorphism and Hamiltonian search are exponential (Hamiltonian <= 24 vertices); maximal-clique output can
  be exponential; Floyd-Warshall/Johnson need N <= 8192.
- Slow beyond about 10 000 vertices: Lex-BFS (O(V^2)), colour refinement on path-like graphs, and
  `graph_adv_graph_isomorphism` (47 s at 10 000 vertices even on an asymmetric graph).
- The matrix needs V^2/8 bytes of presence bits, so it is limited to about 131 000 vertices (weighted) or
  176 000 (unweighted) under the default budget, and costs far more memory than the list on sparse graphs.
- A graph that does not fit in memory can only be analysed by the streaming functions (degrees, components);
  there is no sharded / out-of-core `Graph`.
- A normal `Graph` is not thread-safe (freeze it to share it); library algorithms are not cached.
- Cost-model constants were fitted on one machine; re-run `make crossover` to check yours.
- A backend with a NULL baseline kernel makes that algorithm fail with `GRAPH_ERR_NO_BACKEND_SUPPORT`.

**Pending ideas**
- Presence bits per 64x64 tile (lowering the matrix's V^2/8 floor) and a hybrid dense-core / sparse-tail
  representation.
- A compact read-only snapshot, and a sharded graph for networks larger than memory.
- Cache the heavier library algorithms (PageRank, centralities) per graph version.
- GraphML/GML readers; refit the cost models automatically from a crossover run.
