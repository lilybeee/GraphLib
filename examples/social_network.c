/*
 * social_network.c — a small society grows from six friends to tens of thousands
 * of people, and the graph library is used to study it.
 *
 *   make social CFLAGS="-std=c11 -O2"            # 20,000 people
 *   ./build/social_network 50000 7               # people, random seed
 *
 * What it shows
 *   1. growth    people join in waves; each makes friends by preferential
 *                attachment (popular people attract more), triadic closure
 *                (friends of friends) and a few ties across communities.
 *                Adaptive mode is ON, so the library may change how the graph is
 *                stored while it grows; the report says when that happens.
 *   2. churn     some people leave (graph_remove_vertex); IDs are compacted.
 *   3. analysis  components, degree distribution, clustering, communities
 *                (Louvain, checked against the planted ones), influence
 *                (PageRank), brokers (betweenness on the inner circle), small-
 *                world separation, bridges, k-cores, colouring, cliques.
 *   4. dynamics  an information cascade from the most influential person
 *                versus from an ordinary one.
 *   5. threads   three analyses run at once on a frozen (read-only) graph.
 *   6. compare   list versus matrix timing on this very network.
 *   7. pictures  early network (SVG with names), growing network, final network
 *                by community and by degree, and the adjacency matrix with people
 *                sorted by community.  Written to build/social/.
 */

#define _POSIX_C_SOURCE 200809L

#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "graph_algorithms.h"
#include "graph_viz.h"

#define COMMUNITIES 12
#define OUT_DIR     "build/social/"

/* ------------------------------------------------------------------ */
/*  Small utilities                                                     */
/* ------------------------------------------------------------------ */

static double now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec + (double)t.tv_nsec * 1e-9;
}

static unsigned long long rng_state = 88172645463325252ULL;
static unsigned long long rnd64(void)
{
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return rng_state;
}
static size_t rnd(size_t bound) { return (size_t)(rnd64() % bound); }
static double rnd01(void) { return (double)(rnd64() >> 11) / 9007199254740992.0; }

static void *xmalloc(size_t n)
{
    void *p = malloc(n ? n : 1);
    if (p == NULL) { fputs("out of memory\n", stderr); exit(1); }
    return p;
}

static void check(GraphStatus st, const char *what)
{
    if (st != GRAPH_OK) {
        fprintf(stderr, "%s failed: %s\n", what, graph_status_string(st));
        exit(1);
    }
}

static void heading(const char *text)
{
    printf("\n=== %s ===\n", text);
}

/** A growable int vector. */
typedef struct { int *a; size_t n, cap; } Vec;
static void vec_push(Vec *v, int x)
{
    if (v->n == v->cap) {
        v->cap = v->cap ? v->cap * 2 : 16;
        v->a = realloc(v->a, v->cap * sizeof(int));
        if (v->a == NULL) { fputs("out of memory\n", stderr); exit(1); }
    }
    v->a[v->n++] = x;
}

/* ------------------------------------------------------------------ */
/*  The society                                                         */
/* ------------------------------------------------------------------ */

typedef struct {
    Graph  *g;
    char  **names;                      /* names[id] (owned)               */
    int    *community;                  /* planted community of each id    */
    size_t  cap;                        /* capacity of names[] / community[] */
    Vec     ends[COMMUNITIES];          /* edge endpoints per community: sampling one is
                                           sampling a person proportionally to their degree */
    size_t  members[COMMUNITIES];
} Society;

static const char *SYL_A[] = { "Al", "Be", "Ca", "Da", "El", "Fa", "Gi", "Ha", "Io", "Ja", "Ka", "Li",
                               "Ma", "Ni", "Ol", "Pa", "Qu", "Ra", "Sa", "Ti", "Ur", "Vi", "Wa", "Xe", "Yu", "Za" };
static const char *SYL_B[] = { "ra", "no", "li", "mo", "ta", "ve", "su", "ki", "do", "ne" };
static const char *SYL_C[] = { "n", "s", "th", "ra", "x", "l", "m", "na", "ro", "lyn" };

static char *make_name(size_t serial)
{
    char *s = xmalloc(24);
    snprintf(s, 24, "%s%s%s%zu", SYL_A[rnd(26)], SYL_B[rnd(10)], SYL_C[rnd(10)], serial % 100);
    return s;
}

static void society_reserve(Society *s, size_t need)
{
    if (need <= s->cap) { return; }
    size_t cap = s->cap ? s->cap : 64;
    while (cap < need) { cap *= 2; }
    s->names     = realloc(s->names, cap * sizeof(char *));
    s->community = realloc(s->community, cap * sizeof(int));
    if (s->names == NULL || s->community == NULL) { fputs("out of memory\n", stderr); exit(1); }
    for (size_t i = s->cap; i < cap; ++i) { s->names[i] = NULL; s->community[i] = -1; }
    s->cap = cap;
}

static void befriend(Society *s, int a, int b)
{
    if (a == b || graph_has_edge(s->g, a, b, NULL)) { return; }
    const double strength = 1.0 + (double)rnd(5);                /* 1 = acquaintance .. 5 = best friend */
    if (!graph_add_edge(s->g, a, b, strength)) { return; }
    vec_push(&s->ends[s->community[a]], a);
    vec_push(&s->ends[s->community[b]], b);
}

/** A random person of community c, more likely if they already have many friends. */
static int popular_member(const Society *s, int c)
{
    const Vec *e = &s->ends[c];
    return e->n > 0 ? e->a[rnd(e->n)] : -1;
}

/** Someone who knows `p`: a random friend of p, or -1. */
static int friend_of(Society *s, int p)
{
    const size_t deg = graph_degree(s->g, p);
    if (deg == 0) { return -1; }
    int *nb = xmalloc(deg * sizeof(int));
    (void)graph_neighbors(s->g, p, nb, NULL, deg);
    const int f = nb[rnd(deg)];
    free(nb);
    return f;
}

static int add_person(Society *s, int community)
{
    const int id = graph_add_vertex(s->g);
    if (id < 0) { fputs("could not add a person\n", stderr); exit(1); }
    society_reserve(s, (size_t)id + 1);
    free(s->names[id]);
    s->names[id]     = make_name((size_t)id);
    s->community[id] = community;
    s->members[community]++;
    return id;
}

/** A newcomer joins a community and makes about `m` friends. */
static void join(Society *s)
{
    /* Bigger communities attract more newcomers. */
    size_t total = 0;
    for (int c = 0; c < COMMUNITIES; ++c) { total += s->members[c] + 3; }
    size_t pick = rnd(total);
    int c = 0;
    while (pick >= s->members[c] + 3) { pick -= s->members[c] + 3; ++c; }

    const int me = add_person(s, c);
    const int first = popular_member(s, c);
    if (first < 0 || first == me) {                       /* an empty community: join anyone */
        const int any = (int)rnd((size_t)me + 1);
        if (any != me && graph_is_alive(s->g, any)) { befriend(s, me, any); }
        return;
    }
    befriend(s, me, first);                               /* preferential attachment */
    const int extra = 1 + (int)rnd(3);                    /* 1..3 more friends */
    int anchor = first;
    for (int i = 0; i < extra; ++i) {
        const double r = rnd01();
        if (r < 0.60) {                                   /* triadic closure: a friend of a friend */
            const int ff = friend_of(s, anchor);
            if (ff >= 0) { befriend(s, me, ff); anchor = ff; }
        } else if (r < 0.90) {                            /* another popular member of the same community */
            const int p = popular_member(s, c);
            if (p >= 0) { befriend(s, me, p); }
        } else {                                          /* a bridge to another community */
            const int other = (c + 1 + (int)rnd(COMMUNITIES - 1)) % COMMUNITIES;
            const int p = popular_member(s, other);
            if (p >= 0 && graph_is_alive(s->g, p)) { befriend(s, me, p); }
        }
    }
}

/* ------------------------------------------------------------------ */
/*  Reporting helpers                                                   */
/* ------------------------------------------------------------------ */

static void checkpoint(Society *s, const char *label)
{
    const size_t n = graph_vertex_count(s->g), e = graph_edge_count(s->g);
    printf("  %-14s %7zu people  %8zu friendships  avg friends %5.2f  stored as %-6s (migrations so far: %zu)\n",
           label, n, e, n ? 2.0 * (double)e / (double)n : 0.0,
           graph_rep_name(graph_representation(s->g)), graph_migration_count(s->g));
}

static void bar(size_t value, size_t max, int width)
{
    const int len = max ? (int)((double)value * width / (double)max + 0.5) : 0;
    for (int i = 0; i < len; ++i) { putchar('#'); }
}

static void top_k(const double *score, const size_t n, const int *skip_dead_flag, int k, int *out)
{
    for (int i = 0; i < k; ++i) {
        int best = -1;
        for (size_t v = 0; v < n; ++v) {
            if (skip_dead_flag[v] == 0) { continue; }
            bool taken = false;
            for (int j = 0; j < i; ++j) { taken = taken || out[j] == (int)v; }
            if (!taken && (best < 0 || score[v] > score[best])) { best = (int)v; }
        }
        out[i] = best;
    }
}

static void draw(Graph *g, const char *file, VizOptions o, const char *what)
{
    char path[160];
    snprintf(path, sizeof path, OUT_DIR "%s", file);
    const double t = now();
    const GraphStatus st = graph_visualize(g, path, &o);
    printf("  picture: %-34s %s (%.1fs) - %s\n", path, st == GRAPH_OK ? "ok" : graph_status_string(st), now() - t, what);
}

/* ------------------------------------------------------------------ */
/*  Threads on a frozen graph                                           */
/* ------------------------------------------------------------------ */

typedef struct { Graph *g; double seconds; double result; const char *what; int which; } Job;

static void *worker(void *arg)
{
    Job *j = arg;
    const size_t n = graph_id_bound(j->g);
    const double t = now();
    if (j->which == 0) {
        double *pr = xmalloc(n * sizeof(double));
        check(graph_pagerank(j->g, 0.85, 1e-9, 100, pr, NULL), "pagerank (thread)");
        double best = 0;
        for (size_t v = 0; v < n; ++v) { best = pr[v] > best ? pr[v] : best; }
        j->result = best;
        free(pr);
    } else if (j->which == 1) {
        long long tri;
        check(graph_triangle_count(j->g, &tri), "triangles (thread)");
        j->result = (double)tri;
    } else {
        int *core = xmalloc(n * sizeof(int));
        check(graph_core_numbers(j->g, core), "core numbers (thread)");
        int top = 0;
        for (size_t v = 0; v < n; ++v) { top = core[v] > top ? core[v] : top; }
        j->result = top;
        free(core);
    }
    j->seconds = now() - t;
    return NULL;
}

/* ------------------------------------------------------------------ */
/*  Information cascade                                                 */
/* ------------------------------------------------------------------ */

/** Independent cascade: each newly informed person tells each friend with probability
 *  0.045 * tie strength.  Returns how many people end up informed. */
static size_t cascade(Graph *g, int seed_person, char *informed, int *queue, int *nb, double *wt, size_t maxdeg)
{
    const size_t n = graph_id_bound(g);
    memset(informed, 0, n);
    size_t head = 0, tail = 0;
    queue[tail++] = seed_person;
    informed[seed_person] = 1;
    while (head < tail) {
        const int u = queue[head++];
        const size_t d = graph_neighbors(g, u, nb, wt, maxdeg);
        for (size_t i = 0; i < d; ++i) {
            if (!informed[nb[i]] && rnd01() < 0.045 * wt[i]) {
                informed[nb[i]] = 1;
                queue[tail++] = nb[i];
            }
        }
    }
    return tail;
}

/* ------------------------------------------------------------------ */
/*  Main                                                                */
/* ------------------------------------------------------------------ */

int main(int argc, char **argv)
{
    const size_t target = argc > 1 ? (size_t)strtoul(argv[1], NULL, 10) : 20000;
    if (argc > 2) { rng_state ^= strtoull(argv[2], NULL, 10) * 0x9E3779B97F4A7C15ULL; }
    if (target < 200) { fputs("use at least 200 people\n", stderr); return 2; }
    const double t_start = now();

    Society s;
    memset(&s, 0, sizeof s);
    s.g = graph_create(false);                                    /* undirected, starts as a list, adaptive */
    if (s.g == NULL) { return 1; }

    /* ================= 1. growth ================= */
    heading("1. A society grows");
    printf("  Six friends meet. Then newcomers keep arriving; the library is free to change how\n"
           "  the graph is stored as it grows (adaptive mode is on).\n\n");
    for (int i = 0; i < 6; ++i) { (void)add_person(&s, 0); }
    for (int i = 0; i < 6; ++i) {
        for (int j = i + 1; j < 6; ++j) { if (rnd01() < 0.7) { befriend(&s, i, j); } }
    }
    checkpoint(&s, "the founders");

    const size_t waves[] = { 60, 400, 2000, target / 4, target / 2, target };
    char label[32];
    bool drew_early = false, drew_mid = false;
    for (size_t w = 0; w < sizeof waves / sizeof waves[0]; ++w) {
        while (graph_vertex_count(s.g) < waves[w]) { join(&s); }
        snprintf(label, sizeof label, "wave %zu", w + 1);
        checkpoint(&s, label);
        if (!drew_early && graph_vertex_count(s.g) >= 60) {
            drew_early = true;
            VizOptions o = graph_viz_default_options();
            o.title = "60 people: friendship network (colour = community)";
            o.labels = (const char *const *)s.names;
            o.color = VIZ_COLOR_COMMUNITY;
            o.group_name = "Community";
            o.width = o.height = 1100;
            o.show_labels = 1;
            draw(s.g, "1_early_60_people.svg", o, "the early network, with names");
        }
        if (!drew_mid && graph_vertex_count(s.g) >= 2000) {
            drew_mid = true;
            VizOptions o = graph_viz_default_options();
            o.color = VIZ_COLOR_COMMUNITY;
            o.title = "2000 people";
            o.labels = (const char *const *)s.names;
            o.group_name = "Community";
            o.label_top = 8;
            draw(s.g, "2_growing_2000_people.png", o, "the network while it is still growing");
        }
    }
    const size_t grown = graph_vertex_count(s.g);

    /* ================= 2. churn ================= */
    heading("2. People leave; IDs are compacted");
    size_t left = 0;
    const size_t want_leave = grown / 50;
    while (left < want_leave) {
        const int v = (int)rnd(graph_id_bound(s.g));
        if (graph_is_alive(s.g, v)) {
            if (graph_remove_vertex(s.g, v)) {
                s.members[s.community[v]]--;
                free(s.names[v]);
                s.names[v] = NULL;
                s.community[v] = -1;
                ++left;
            }
        }
    }
    printf("  %zu people left (2%%). The graph now has %zu people but IDs still run up to %zu\n"
           "  (the library only advises compaction once 25%% of IDs are dead: advised = %s). We compact anyway.\n",
           left, graph_vertex_count(s.g), graph_id_bound(s.g) - 1, graph_needs_compaction(s.g) ? "yes" : "no");
    {
        const size_t old_bound = graph_id_bound(s.g);
        int *remap = xmalloc(old_bound * sizeof(int));
        check(graph_compact(s.g, remap) ? GRAPH_OK : graph_last_error(), "graph_compact");
        char **names = xmalloc(old_bound * sizeof(char *));
        int *comm = xmalloc(old_bound * sizeof(int));
        for (size_t i = 0; i < old_bound; ++i) { names[i] = NULL; comm[i] = -1; }
        for (size_t old = 0; old < old_bound; ++old) {
            if (remap[old] >= 0) { names[remap[old]] = s.names[old]; comm[remap[old]] = s.community[old]; }
        }
        free(s.names); free(s.community); free(remap);
        s.names = names; s.community = comm; s.cap = old_bound;
        for (int c = 0; c < COMMUNITIES; ++c) { s.ends[c].n = 0; }   /* the endpoint lists are stale now */
        printf("  After graph_compact: %zu people, IDs 0..%zu (names and communities remapped).\n",
               graph_vertex_count(s.g), graph_id_bound(s.g) - 1);
    }
    checkpoint(&s, "final society");
    Graph *g = s.g;
    const size_t n = graph_id_bound(g);          /* == vertex count after compaction */

    /* ================= 3. analysis ================= */
    heading("3. Analysis");
    double t;

    /* connectivity */
    int *comp = xmalloc(n * sizeof(int));
    t = now();
    const size_t ncomp = graph_connected_components(g, comp);
    size_t *csize = calloc(ncomp ? ncomp : 1, sizeof(size_t));
    for (size_t v = 0; v < n; ++v) { csize[comp[v]]++; }
    size_t giant = 0;
    for (size_t c = 0; c < ncomp; ++c) { giant = csize[c] > giant ? csize[c] : giant; }
    printf("  Connectivity   %zu component(s); the largest holds %.1f%% of everyone (%.2fs)\n",
           ncomp, 100.0 * (double)giant / (double)n, now() - t);
    free(csize);

    /* degree distribution */
    size_t hist[12] = { 0 };
    size_t maxdeg = 0;
    for (size_t v = 0; v < n; ++v) {
        const size_t d = graph_degree(g, (int)v);
        maxdeg = d > maxdeg ? d : maxdeg;
        size_t b = 0;
        while (((size_t)1 << (b + 1)) <= d && b < 11) { ++b; }
        hist[b]++;
    }
    printf("  Friends per person (log-scale buckets); the most-connected person has %zu friends:\n", maxdeg);
    size_t hmax = 0;
    for (int b = 0; b < 12; ++b) { hmax = hist[b] > hmax ? hist[b] : hmax; }
    for (int b = 0; b < 12; ++b) {
        if (hist[b] == 0) { continue; }
        printf("     %4zu-%-5zu %6zu people ", (size_t)1 << b, ((size_t)1 << (b + 1)) - 1, hist[b]);
        bar(hist[b], hmax, 40);
        putchar('\n');
    }

    /* clustering, assortativity, stats */
    long long tri = 0;
    double avgcl = 0, assort = 0;
    t = now();
    check(graph_triangle_count(g, &tri), "triangle_count");
    printf("  Triangles      %lld closed friend-triples (%.2fs)\n", tri, now() - t);
    double *local = xmalloc(n * sizeof(double));
    t = now();
    check(graph_clustering_coefficient(g, local, &avgcl), "clustering");
    printf("  Clustering     %.3f: a person's friends are friends with each other this often (%.2fs)\n", avgcl, now() - t);
    free(local);
    if (graph_degree_assortativity(g, &assort) == GRAPH_OK) {
        printf("  Assortativity  %+.3f (%s)\n", assort, assort > 0.05 ? "popular people befriend popular people"
                                                       : assort < -0.05 ? "popular people befriend newcomers" : "no strong preference");
    }

    /* communities */
    int *louv = xmalloc(n * sizeof(int));
    double Q = 0;
    size_t nlouv = 0;
    t = now();
    check(graph_louvain(g, louv, &Q, &nlouv), "louvain");
    {
        /* purity: how many people sit in a found community whose majority matches their planted one */
        size_t *tab = calloc(nlouv * COMMUNITIES, sizeof(size_t));
        for (size_t v = 0; v < n; ++v) { tab[(size_t)louv[v] * COMMUNITIES + (size_t)s.community[v]]++; }
        size_t agree = 0;
        for (size_t l = 0; l < nlouv; ++l) {
            size_t best = 0;
            for (int c = 0; c < COMMUNITIES; ++c) { best = tab[l * COMMUNITIES + (size_t)c] > best ? tab[l * COMMUNITIES + (size_t)c] : best; }
            agree += best;
        }
        free(tab);
        printf("  Communities    Louvain found %zu groups (modularity %.3f); %.1f%% of people sit in a group dominated by\n"
               "                 their true community (%d were planted) (%.2fs)\n",
               nlouv, Q, 100.0 * (double)agree / (double)n, COMMUNITIES, now() - t);
    }

    /* influence */
    double *pr = xmalloc(n * sizeof(double));
    t = now();
    size_t iters = 0;
    check(graph_pagerank(g, 0.85, 1e-10, 200, pr, &iters), "pagerank");
    int *alive = xmalloc(n * sizeof(int));
    for (size_t v = 0; v < n; ++v) { alive[v] = 1; }
    int top5[5];
    top_k(pr, n, alive, 5, top5);
    printf("  Influence      PageRank converged in %zu rounds (%.2fs). The five most influential people:\n", iters, now() - t);
    for (int i = 0; i < 5; ++i) {
        printf("     %d. %-12s %5zu friends, community %2d, rank %.5f\n", i + 1, s.names[top5[i]],
               graph_degree(g, top5[i]), s.community[top5[i]], pr[top5[i]]);
    }

    /* brokers among the inner circle (betweenness is O(V*E): run it on the top people) */
    {
        const size_t circle = n < 1500 ? n : 1500;
        int *ids = xmalloc(circle * sizeof(int));
        int *all = xmalloc(n * sizeof(int));
        for (size_t v = 0; v < n; ++v) { all[v] = 1; }
        top_k(pr, n, all, (int)circle > 1500 ? 1500 : (int)circle, ids);
        Graph *core = NULL;
        check(graph_induced_subgraph(g, ids, circle, &core), "inner circle");
        double *btw = xmalloc(circle * sizeof(double));
        t = now();
        check(graph_betweenness_centrality(core, btw), "betweenness");
        int *alive2 = xmalloc(circle * sizeof(int));
        for (size_t v = 0; v < circle; ++v) { alive2[v] = 1; }
        int top3[3];
        top_k(btw, circle, alive2, 3, top3);
        printf("  Brokers        inner circle = the %zu most influential (%zu friendships among them, %.2fs). Brokers, who sit on\n"
               "                 the most shortest paths between others: %s, %s, %s\n",
               circle, graph_edge_count(core), now() - t, s.names[ids[top3[0]]], s.names[ids[top3[1]]], s.names[ids[top3[2]]]);
        int *cl = xmalloc(circle * sizeof(int));
        size_t clique = 0;
        t = now();
        if (graph_max_clique(core, cl, &clique) == GRAPH_OK) {
            printf("                 The largest group of mutual friends in the circle has %zu people (%.2fs): ", clique, now() - t);
            for (size_t i = 0; i < clique && i < 8; ++i) { printf("%s%s", i ? ", " : "", s.names[ids[cl[i]]]); }
            putchar('\n');
        }
        int *col = xmalloc(circle * sizeof(int));
        size_t tables = 0;
        if (graph_color_dsatur(core, col, &tables) == GRAPH_OK) {
            printf("                 Seating plan: %zu tables are enough so that no two friends share a table (DSatur colouring)\n", tables);
        }
        free(ids); free(all); free(btw); free(alive2); free(cl); free(col);
        graph_destroy(core);
    }

    /* cohesion */
    {
        int *core = xmalloc(n * sizeof(int));
        check(graph_core_numbers(g, core), "core numbers");
        int top = 0;
        size_t in_top = 0;
        for (size_t v = 0; v < n; ++v) { top = core[v] > top ? core[v] : top; }
        for (size_t v = 0; v < n; ++v) { in_top += core[v] == top; }
        printf("  Cohesion       the densest nucleus is the %d-core: %zu people who each have >= %d friends within it\n", top, in_top, top);
        free(core);
        bool *cut = xmalloc(n * sizeof(bool));
        size_t ncut = 0;
        check(graph_articulation_points(g, cut, &ncut), "articulation points");
        GraphEdge *br = xmalloc((graph_edge_count(g) + 1) * sizeof(GraphEdge));
        size_t nbr = 0;
        check(graph_bridges(g, br, &nbr), "bridges");
        printf("  Fragility      %zu people and %zu friendships are single points of failure (removing them disconnects someone)\n", ncut, nbr);
        free(cut); free(br);
    }

    /* small world */
    {
        const int samples = 60;
        int *order = xmalloc(n * sizeof(int)), *parent = xmalloc(n * sizeof(int)), *depth = xmalloc(n * sizeof(int));
        double sum = 0;
        size_t pairs = 0;
        int far = 0;
        size_t hist_h[16] = { 0 };
        t = now();
        for (int k = 0; k < samples; ++k) {
            const int src = (int)rnd(n);
            const size_t got = graph_bfs(g, src, order, parent);
            depth[src] = 0;
            for (size_t i = 1; i < got; ++i) {
                depth[order[i]] = depth[parent[order[i]]] + 1;
                sum += depth[order[i]];
                ++pairs;
                far = depth[order[i]] > far ? depth[order[i]] : far;
                hist_h[depth[order[i]] < 15 ? depth[order[i]] : 15]++;
            }
        }
        printf("  Small world    from %d random people, everyone reachable is on average %.2f handshakes away; the farthest seen is %d (%.2fs)\n",
               samples, pairs ? sum / (double)pairs : 0.0, far, now() - t);
        printf("                 handshakes:");
        for (int h = 1; h <= far && h < 15; ++h) { printf("  %d:%.0f%%", h, 100.0 * (double)hist_h[h] / (double)(pairs ? pairs : 1)); }
        putchar('\n');
        /* the most distant pair found: a concrete chain of friends */
        {
            const int a = (int)rnd(n);
            const size_t got = graph_bfs(g, a, order, parent);
            const int b = order[got - 1];
            int *path = xmalloc(n * sizeof(int));
            const size_t plen = graph_path_from_parent(parent, n, a, b, path, n);
            printf("                 A chain from %s to the person farthest from them, %s:\n                 ", s.names[a], s.names[b]);
            for (size_t i = 0; i < plen; ++i) { printf("%s%s", i ? " -> " : "", s.names[path[i]]); }
            putchar('\n');
            /* strongest chain: shortest path where strong ties are short (1 / strength) is not directly
             * available with this weight convention, so show the three cheapest routes by total weight */
            GraphPath routes[3];
            size_t found = 0;
            if (graph_k_shortest_paths(g, a, b, 3, routes, &found) == GRAPH_OK) {
                printf("                 The %zu lowest-total-weight routes between them have %.0f, %.0f, %.0f units of weight\n", found,
                       routes[0].cost, found > 1 ? routes[1].cost : 0.0, found > 2 ? routes[2].cost : 0.0);
                for (size_t i = 0; i < found; ++i) { graph_path_free(&routes[i]); }
            }
            free(path);
        }
        free(order); free(parent); free(depth);
    }

    /* ================= 4. dynamics ================= */
    heading("4. A rumour spreads");
    {
        char *informed = xmalloc(n);
        int *queue = xmalloc(n * sizeof(int));
        int *nb = xmalloc((maxdeg + 1) * sizeof(int));
        double *wt = xmalloc((maxdeg + 1) * sizeof(double));
        /* an ordinary person: median degree */
        int ordinary = 0;
        size_t want = 0;
        for (size_t v = 0; v < n; ++v) { want += graph_degree(g, (int)v); }
        want /= n;
        for (size_t v = 0; v < n; ++v) { if (graph_degree(g, (int)v) == want) { ordinary = (int)v; break; } }
        const int trials = 30;
        double reach_top = 0, reach_ord = 0, reach_rand = 0;
        for (int i = 0; i < trials; ++i) {
            reach_top  += (double)cascade(g, top5[0], informed, queue, nb, wt, maxdeg + 1);
            reach_ord  += (double)cascade(g, ordinary, informed, queue, nb, wt, maxdeg + 1);
            reach_rand += (double)cascade(g, (int)rnd(n), informed, queue, nb, wt, maxdeg + 1);
        }
        printf("  Each informed person tells each friend with probability 4.5%% x friendship strength. Average reach of %d runs:\n", trials);
        printf("     started by %-12s (most influential, %zu friends): %8.0f people (%.1f%%)\n", s.names[top5[0]], graph_degree(g, top5[0]),
               reach_top / trials, 100.0 * reach_top / trials / (double)n);
        printf("     started by %-12s (average person, %zu friends)  : %8.0f people (%.1f%%)\n", s.names[ordinary], graph_degree(g, ordinary),
               reach_ord / trials, 100.0 * reach_ord / trials / (double)n);
        printf("     started by a random person                                : %8.0f people (%.1f%%)\n", reach_rand / trials,
               100.0 * reach_rand / trials / (double)n);
        free(informed); free(queue); free(nb); free(wt);
    }

    /* ================= 5. threads ================= */
    heading("5. Three analyses at the same time on a frozen graph");
    {
        graph_freeze(g);
        Job jobs[3] = { { g, 0, 0, "PageRank", 0 }, { g, 0, 0, "triangle count", 1 }, { g, 0, 0, "k-cores", 2 } };
        pthread_t th[3];
        const double t0 = now();
        for (int i = 0; i < 3; ++i) { pthread_create(&th[i], NULL, worker, &jobs[i]); }
        for (int i = 0; i < 3; ++i) { pthread_join(th[i], NULL); }
        const double wall = now() - t0;
        double serial = 0;
        for (int i = 0; i < 3; ++i) { serial += jobs[i].seconds; }
        printf("  graph_freeze() makes the handle read-only and safe to share. PageRank (max rank %.5f), triangle count (%.0f) and\n"
               "  k-cores (deepest core %.0f) ran at the same time: %.2fs wall time, %.2fs of work in total.\n",
               jobs[0].result, jobs[1].result, jobs[2].result, wall, serial);
        graph_clear_error();
        const bool refused = graph_add_vertex(g) < 0;
        printf("  While frozen, adding a person is refused: %s (%s)\n", refused ? "yes" : "NO", graph_status_string(graph_last_error()));
        graph_thaw(g);
    }

    /* ================= 6. compare ================= */
    heading("6. List or matrix? The library can measure it on this network");
    {
        printf("  Current storage: %s. graph_compare() copies the graph into both representations and times the same algorithm:\n",
               graph_rep_name(graph_representation(g)));
        const GraphAlgo algos[] = { GRAPH_ALGO_BFS, GRAPH_ALGO_PAGERANK, GRAPH_ALGO_TRIANGLES, GRAPH_ALGO_CORE_NUMBERS };
        for (size_t i = 0; i < sizeof algos / sizeof algos[0]; ++i) {
            GraphCompareResult r;
            if (graph_compare(g, algos[i], 0, 2, &r)) {
                printf("     %-16s list %9.2f ms   matrix %9.2f ms   -> %s %.1fx faster; memory list %.1f MB vs matrix %.1f MB; results %s\n",
                       graph_algo_name(algos[i]), r.list_ns / 1e6, r.matrix_ns / 1e6, graph_rep_name(r.faster),
                       r.list_ns > r.matrix_ns ? r.list_ns / r.matrix_ns : r.matrix_ns / r.list_ns,
                       (double)r.list_bytes / 1e6, (double)r.matrix_bytes / 1e6, r.agreed ? "agree" : "DIFFER");
            }
        }
        printf("  A social network is sparse (%.3f%% of possible friendships exist), so the list is the natural home;\n"
               "  the adaptive layer kept it there: migrations during the whole run = %zu.\n",
               100.0 * 2.0 * (double)graph_edge_count(g) / ((double)n * (double)(n - 1)), graph_migration_count(g));
    }

    /* ================= 7. pictures ================= */
    heading("7. Pictures");
    {
        VizOptions o = graph_viz_default_options();
        o.color = VIZ_COLOR_COMMUNITY;
        o.title = "The whole society, coloured by community";
        o.labels = (const char *const *)s.names;       /* names of the most-connected people */
        o.group_name = "Community";
        o.label_top = 14;
        o.width = 2600; o.height = 1900;
        draw(g, "3_final_by_community.png", o, "colour = community found by Louvain");
        o.color = VIZ_COLOR_DEGREE;
        o.title = "Popularity: number of friends (temperature scale)";
        o.label_top = 18;
        o.width = 2600; o.height = 1900;
        draw(g, "4_final_by_degree.png", o, "red = very popular, blue = few friends");

        /* Re-number people so each community is contiguous, then draw the adjacency matrix. */
        size_t *start = calloc(nlouv + 1, sizeof(size_t));
        for (size_t v = 0; v < n; ++v) { start[louv[v] + 1]++; }
        for (size_t c = 0; c < nlouv; ++c) { start[c + 1] += start[c]; }
        int *sorted = xmalloc(n * sizeof(int));
        size_t *fill = xmalloc((nlouv + 1) * sizeof(size_t));
        memcpy(fill, start, nlouv * sizeof(size_t));
        for (size_t v = 0; v < n; ++v) { sorted[fill[louv[v]]++] = (int)v; }
        Graph *ordered = NULL;
        check(graph_induced_subgraph(g, sorted, n, &ordered), "reorder by community");
        VizOptions m = graph_viz_default_options();
        m.layout = VIZ_LAYOUT_MATRIX;
        m.width = m.height = 1400;
        draw(ordered, "5_adjacency_matrix_by_community.png", m, "who knows whom; people sorted by community, so groups form blocks");
        graph_destroy(ordered);
        free(start); free(sorted); free(fill);
    }

    heading("Done");
    printf("  %zu people, %zu friendships, %.1f seconds in total. Pictures are in " OUT_DIR "\n", n, graph_edge_count(g), now() - t_start);

    for (size_t i = 0; i < s.cap; ++i) { free(s.names[i]); }
    free(s.names); free(s.community);
    for (int c = 0; c < COMMUNITIES; ++c) { free(s.ends[c].a); }
    free(comp); free(louv); free(pr); free(alive);
    graph_destroy(g);
    return 0;
}
