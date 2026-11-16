/*
 * adjacency_list.c — Adjacency-list graph backend
 *
 * Implements all functions declared in adjacency_list.h.
 *
 * Internal representation
 * -----------------------
 *   Vertex records ("slots") are stored in fixed-size segments that are reached
 *   through a small directory:
 *
 *     segs[id >> LIST_SEG_SHIFT][id & SEG_MASK]
 *
 *   Segments are never moved or reallocated, so growing the graph costs one
 *   calloc and compaction can release whole segments.  A slot holds one
 *   contiguous block per vertex: the weights come first and the neighbour IDs
 *   follow, so a traversal that ignores weights touches 4 bytes per edge.  Each row
 *   is kept sorted by neighbour ID, which makes has_edge a binary search and
 *   appending edges in ascending order O(1).
 *
 *   Removed vertices keep their ID on a free list (every dead ID is on it, so
 *   dead count == n_free).  add_vertex pops from that list before growing.
 *
 * Decisions worth knowing about
 * -----------------------------
 *   Growth slack.  A block doubles when full, so on average about a third of the entry
 *   memory is spare room.  list_graph_shrink_to_fit() gives it back once a graph is finished.
 *
 *   Unweighted mode.  An unweighted graph (list_graph_create_unweighted) has no weight array:
 *   the block holds only neighbour IDs (4 bytes per entry instead of 12) and every edge reads
 *   back as 1.0.  Slot.w == NULL is the marker; every place that touches weights checks it.
 *
 *   Memory budget.  Each graph remembers the budget in force when it was created
 *   (memory_budget.h) and counts the bytes its blocks hold (entry_bytes).  A growth step
 *   that would pass the budget fails with GRAPH_ERR_TOO_LARGE and leaves the graph unchanged,
 *   so a too-big request is an error the program can handle, not an out-of-memory kill.
 *
 * Compile: cc -std=c11 -Wall -Wextra -Werror -pedantic adjacency_list.c -lm
 */

#include "adjacency_list.h"
#include "graph_status.h"
#include "memory_budget.h"

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef LIST_SEG_SHIFT
#define LIST_SEG_SHIFT 10
#endif
#define SEG_SIZE ((size_t)1 << LIST_SEG_SHIFT)
#define SEG_MASK (SEG_SIZE - 1)

/* ------------------------------------------------------------------ */
/*  Data structures                                                     */
/* ------------------------------------------------------------------ */

typedef struct {
    double  *w;     /* weights; also the start of the block (NULL if cap == 0) */
    int     *to;    /* neighbour IDs, stored after the weights in the block   */
    uint32_t deg;
    uint32_t cap;
    bool     alive;
} Slot;

/** The start of a slot's allocation: the weights if there are any, else the neighbour IDs. */
static inline void *slot_block(const Slot *s)
{
    return s->w != NULL ? (void *)s->w : (void *)s->to;
}

struct ListGraph {
    Slot   **segs;       /* directory of segments */
    size_t   n_segs;
    size_t   segs_cap;
    size_t   id_bound;   /* IDs 0 .. id_bound-1 exist (alive or dead) */
    int     *free_ids;   /* dead IDs available for reuse */
    size_t   n_free;
    size_t   free_cap;
    size_t   n_edges;    /* undirected edges counted once */
    size_t   neg_entries;/* stored adjacency entries with weight < 0 */
    uint64_t version;
    size_t   budget;       /* bytes this graph may use (0 = unlimited); see memory_budget.h */
    size_t   entry_bytes;  /* bytes held by the adjacency blocks */
    bool     directed;
    bool     unweighted;   /* no weight array is stored; every edge reads back as 1.0 */
};

/* ------------------------------------------------------------------ */
/*  Internal helpers                                                    */
/* ------------------------------------------------------------------ */

/* Record a failure code, then return the given value. */
#define FAIL(code, value) do { graph_set_error(code); return (value); } while (0)
/* A failed allocation is OUT_OF_MEMORY unless the budget check already said TOO_LARGE. */
#define FAIL_ALLOC(value) do {                                         \
        if (graph_last_error() != GRAPH_ERR_TOO_LARGE) {               \
            graph_set_error(GRAPH_ERR_OUT_OF_MEMORY);                  \
        }                                                              \
        return (value);                                                \
    } while (0)
#define BAD_VERTEX(graph) ((graph) == NULL ? GRAPH_ERR_INVALID_ARGUMENT : GRAPH_ERR_INVALID_VERTEX)

static inline Slot *slot_at(const ListGraph *graph, int v)
{
    return &graph->segs[(size_t)v >> LIST_SEG_SHIFT][(size_t)v & SEG_MASK];
}

static inline bool valid_vertex(const ListGraph *graph, int v)
{
    return graph != NULL && v >= 0 && (size_t)v < graph->id_bound &&
           slot_at(graph, v)->alive;
}

static inline size_t live_count(const ListGraph *graph)
{
    return graph->id_bound - graph->n_free;
}

static inline bool bit_get(const uint64_t *bits, size_t i)
{
    return (bits[i >> 6] >> (i & 63)) & 1u;
}

static inline void bit_set(uint64_t *bits, size_t i)
{
    bits[i >> 6] |= (uint64_t)1 << (i & 63);
}

/**
 * find_index - binary search for neighbour v in s (rows are sorted by ID).
 * Returns true and the position if found; otherwise false and *index is the
 * position where v would be inserted.  O(log degree).
 */
static bool find_index(const Slot *s, int v, size_t *index)
{
    size_t lo = 0, hi = s->deg;
    while (lo < hi) {
        const size_t mid = lo + (hi - lo) / 2;
        if (s->to[mid] < v) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    *index = lo;
    return lo < s->deg && s->to[lo] == v;
}

/** slot_reserve - make room for at least 'need' entries (doubling growth). */
static bool slot_reserve(ListGraph *graph, Slot *s, uint32_t need)
{
    if (need <= s->cap) {
        return true;
    }
    uint32_t cap = s->cap != 0 ? s->cap : 4;
    while (cap < need) {
        if (cap > UINT32_MAX / 2) {
            return false;
        }
        cap *= 2;
    }

    /* Bytes per stored neighbour: an ID, plus a weight unless the graph is unweighted. */
    const size_t entry = graph->unweighted ? sizeof(int) : sizeof(double) + sizeof(int);
    /* During the copy both the old and the new block exist, so the budget must cover both. */
    if (graph->budget != 0 &&
        graph->n_segs * SEG_SIZE * sizeof(Slot) + graph->entry_bytes + (size_t)cap * entry > graph->budget) {
        graph_set_error(GRAPH_ERR_TOO_LARGE);
        return false;
    }
    void *block = malloc((size_t)cap * entry);
    if (block == NULL) {
        return false;
    }
    graph->entry_bytes += (size_t)cap * entry - (size_t)s->cap * entry;
    double *w  = graph->unweighted ? NULL : (double *)block;
    int    *to = graph->unweighted ? (int *)block : (int *)((double *)block + cap);
    if (s->deg != 0) {
        if (w != NULL) {
            memcpy(w, s->w, (size_t)s->deg * sizeof(double));
        }
        memcpy(to, s->to, (size_t)s->deg * sizeof(int));
    }
    free(slot_block(s));
    s->w   = w;
    s->to  = to;
    s->cap = cap;
    return true;
}

/**
 * insert_entry - put (v, w) at position pos, shifting the tail up by one so the
 * row stays sorted.  Capacity must already be reserved.  Appending (pos == deg)
 * is O(1), which is the common case when edges are added in ascending order.
 */
static void insert_entry(ListGraph *graph, Slot *s, size_t pos, int v, double w)
{
    const size_t tail = (size_t)s->deg - pos;
    if (tail != 0) {
        memmove(s->to + pos + 1, s->to + pos, tail * sizeof(*s->to));
        if (s->w != NULL) {
            memmove(s->w + pos + 1, s->w + pos, tail * sizeof(*s->w));
        }
    }
    s->to[pos] = v;
    if (s->w != NULL) {
        s->w[pos] = w;
    }
    ++s->deg;
    if (w < 0.0) {
        ++graph->neg_entries;
    }
}

/** erase_entry - remove entry i from s, shifting the tail down (keeps the row sorted). */
static void erase_entry(ListGraph *graph, Slot *s, size_t i)
{
    if (s->w != NULL && s->w[i] < 0.0) {
        --graph->neg_entries;
    }
    const size_t tail = (size_t)s->deg - i - 1;
    if (tail != 0) {
        memmove(s->to + i, s->to + i + 1, tail * sizeof(*s->to));
        if (s->w != NULL) {
            memmove(s->w + i, s->w + i + 1, tail * sizeof(*s->w));
        }
    }
    --s->deg;
}

/** ensure_segments - make the directory cover IDs below 'bound'. */
static bool ensure_segments(ListGraph *graph, size_t bound)
{
    const size_t need = (bound + SEG_SIZE - 1) >> LIST_SEG_SHIFT;

    if (need > graph->segs_cap) {
        size_t cap = graph->segs_cap != 0 ? graph->segs_cap : 4;
        while (cap < need) {
            cap *= 2;
        }
        Slot **dir = realloc(graph->segs, cap * sizeof(*dir));
        if (dir == NULL) {
            return false;
        }
        graph->segs     = dir;
        graph->segs_cap = cap;
    }
    if (need > graph->n_segs && graph->budget != 0 &&
        need * SEG_SIZE * sizeof(Slot) + graph->entry_bytes > graph->budget) {
        graph_set_error(GRAPH_ERR_TOO_LARGE);
        return false;
    }
    while (graph->n_segs < need) {
        Slot *seg = calloc(SEG_SIZE, sizeof(*seg));
        if (seg == NULL) {
            return false;
        }
        graph->segs[graph->n_segs++] = seg;
    }
    return true;
}

/* ------------------------------------------------------------------ */
/*  Construction / destruction                                          */
/* ------------------------------------------------------------------ */

static ListGraph *create_graph(size_t vertices, bool directed, bool unweighted)
{
    if (vertices > (size_t)INT_MAX) {
        return NULL;
    }
    ListGraph *graph = calloc(1, sizeof(*graph));
    if (graph == NULL) {
        return NULL;
    }
    graph->directed   = directed;
    graph->unweighted = unweighted;
    graph->budget     = graph_memory_budget();

    if (!ensure_segments(graph, vertices)) {
        list_graph_destroy(graph); /* id_bound is still 0: frees segments only */
        if (graph_last_error() != GRAPH_ERR_TOO_LARGE) {
            graph_set_error(GRAPH_ERR_OUT_OF_MEMORY);
        }
        return NULL;
    }
    for (size_t v = 0; v < vertices; ++v) {
        slot_at(graph, (int)v)->alive = true;
    }
    graph->id_bound = vertices;
    return graph;
}

ListGraph *list_graph_create(size_t vertices, bool directed)
{
    return create_graph(vertices, directed, false);
}

ListGraph *list_graph_create_unweighted(size_t vertices, bool directed)
{
    return create_graph(vertices, directed, true);
}

bool list_graph_is_unweighted(const ListGraph *graph)
{
    return graph != NULL && graph->unweighted;
}

void list_graph_destroy(ListGraph *graph)
{
    if (graph == NULL) {
        return;
    }
    for (size_t v = 0; v < graph->id_bound; ++v) {
        free(slot_block(slot_at(graph, (int)v)));
    }
    for (size_t s = 0; s < graph->n_segs; ++s) {
        free(graph->segs[s]);
    }
    free(graph->segs);
    free(graph->free_ids);
    free(graph);
}

/* ------------------------------------------------------------------ */
/*  Queries                                                             */
/* ------------------------------------------------------------------ */

size_t list_graph_vertex_count(const ListGraph *graph)
{
    return graph == NULL ? 0 : live_count(graph);
}

size_t list_graph_id_bound(const ListGraph *graph)
{
    return graph == NULL ? 0 : graph->id_bound;
}

bool list_graph_is_alive(const ListGraph *graph, int v)
{
    return valid_vertex(graph, v);
}

size_t list_graph_edge_count(const ListGraph *graph)
{
    return graph == NULL ? 0 : graph->n_edges;
}

bool list_graph_is_directed(const ListGraph *graph)
{
    return graph != NULL && graph->directed;
}

size_t list_graph_degree(const ListGraph *graph, int v)
{
    return valid_vertex(graph, v) ? slot_at(graph, v)->deg : 0;
}

size_t list_graph_neighbors(const ListGraph *graph, int v,
                            const int **to, const double **weight)
{
    if (!valid_vertex(graph, v)) {
        if (to != NULL)     { *to = NULL; }
        if (weight != NULL) { *weight = NULL; }
        return 0;
    }
    const Slot *s = slot_at(graph, v);
    if (to != NULL)     { *to = s->to; }
    if (weight != NULL) { *weight = s->w; }
    return s->deg;
}

bool list_graph_next_neighbor(const ListGraph *graph, int v, size_t *cursor,
                              int *to, double *weight)
{
    if (!valid_vertex(graph, v) || cursor == NULL || to == NULL) {
        return false;
    }
    const Slot *s = slot_at(graph, v);
    if (*cursor >= s->deg) {
        return false;
    }
    *to = s->to[*cursor];
    if (weight != NULL) {
        *weight = s->w != NULL ? s->w[*cursor] : 1.0;
    }
    ++*cursor;
    return true;
}

bool list_graph_shrink_to_fit(ListGraph *graph)
{
    if (graph == NULL) {
        FAIL(GRAPH_ERR_INVALID_ARGUMENT, false);
    }
    const size_t entry = graph->unweighted ? sizeof(int) : sizeof(double) + sizeof(int);
    bool ok = true;
    for (size_t v = 0; v < graph->id_bound; ++v) {
        Slot *s = slot_at(graph, (int)v);
        if (!s->alive || s->cap == s->deg) {
            continue;
        }
        if (s->deg == 0) {                       /* nothing stored: give the whole block back */
            graph->entry_bytes -= (size_t)s->cap * entry;
            free(slot_block(s));
            s->w = NULL;
            s->to = NULL;
            s->cap = 0;
            continue;
        }
        void *block = malloc((size_t)s->deg * entry);
        if (block == NULL) {                     /* keep the slot as it was; it is still valid */
            ok = false;
            continue;
        }
        double *w  = graph->unweighted ? NULL : (double *)block;
        int    *to = graph->unweighted ? (int *)block : (int *)((double *)block + s->deg);
        if (w != NULL) {
            memcpy(w, s->w, (size_t)s->deg * sizeof(double));
        }
        memcpy(to, s->to, (size_t)s->deg * sizeof(int));
        free(slot_block(s));
        graph->entry_bytes -= (size_t)(s->cap - s->deg) * entry;
        s->w   = w;
        s->to  = to;
        s->cap = s->deg;
    }
    if (!ok) {
        graph_set_error(GRAPH_ERR_OUT_OF_MEMORY);
    }
    return ok;
}

uint64_t list_graph_version(const ListGraph *graph)
{
    return graph == NULL ? 0 : graph->version;
}

size_t list_graph_memory_bytes(const ListGraph *graph)
{
    if (graph == NULL) {
        return 0;
    }
    size_t bytes = sizeof(*graph) + graph->segs_cap * sizeof(*graph->segs) +
                   graph->n_segs * SEG_SIZE * sizeof(Slot) +
                   graph->free_cap * sizeof(*graph->free_ids);
    for (size_t v = 0; v < graph->id_bound; ++v) {
        bytes += (size_t)slot_at(graph, (int)v)->cap * (graph->unweighted ? sizeof(int) : sizeof(double) + sizeof(int));
    }
    return bytes;
}

/* ------------------------------------------------------------------ */
/*  Vertex modification                                                 */
/* ------------------------------------------------------------------ */

int list_graph_add_vertex(ListGraph *graph)
{
    if (graph == NULL) {
        FAIL(GRAPH_ERR_INVALID_ARGUMENT, -1);
    }
    int id;
    if (graph->n_free != 0) {
        id = graph->free_ids[--graph->n_free];
    } else {
        if (graph->id_bound >= (size_t)INT_MAX) {
            FAIL(GRAPH_ERR_TOO_LARGE, -1);
        }
        if (!ensure_segments(graph, graph->id_bound + 1)) {
            FAIL_ALLOC(-1);
        }
        id = (int)graph->id_bound++;
    }
    slot_at(graph, id)->alive = true;
    ++graph->version;
    return id;
}

bool list_graph_remove_vertex(ListGraph *graph, int v)
{
    if (!valid_vertex(graph, v)) {
        FAIL(BAD_VERTEX(graph), false);
    }

    /* Reserve the free-list slot first so a failure leaves the graph intact. */
    if (graph->n_free == graph->free_cap) {
        const size_t cap = graph->free_cap != 0 ? graph->free_cap * 2 : 16;
        int *ids = realloc(graph->free_ids, cap * sizeof(*ids));
        if (ids == NULL) {
            FAIL(GRAPH_ERR_OUT_OF_MEMORY, false);
        }
        graph->free_ids = ids;
        graph->free_cap = cap;
    }

    Slot *sv = slot_at(graph, v);

    if (!graph->directed) {
        /* Every edge is mirrored in the neighbour's list. */
        for (size_t i = 0; i < sv->deg; ++i) {
            const int w = sv->to[i];
            size_t j;
            if (w != v && find_index(slot_at(graph, w), v, &j)) {
                erase_entry(graph, slot_at(graph, w), j);
            }
        }
    } else {
        /* No reverse index: scan every list for in-edges. */
        for (size_t u = 0; u < graph->id_bound; ++u) {
            Slot *su = slot_at(graph, (int)u);
            size_t j;
            if (su->alive && u != (size_t)v && find_index(su, v, &j)) {
                erase_entry(graph, su, j);
                --graph->n_edges;
            }
        }
    }

    graph->n_edges -= sv->deg;
    for (size_t i = 0; sv->w != NULL && i < sv->deg; ++i) {
        if (sv->w[i] < 0.0) {
            --graph->neg_entries;
        }
    }
    graph->entry_bytes -= (size_t)sv->cap * (graph->unweighted ? sizeof(int) : sizeof(double) + sizeof(int));
    free(slot_block(sv));
    memset(sv, 0, sizeof(*sv)); /* alive = false */

    graph->free_ids[graph->n_free++] = v;
    ++graph->version;
    return true;
}

bool list_graph_needs_compaction(const ListGraph *graph)
{
    return graph != NULL && graph->n_free != 0 &&
           graph->n_free * 4 > graph->id_bound;
}

bool list_graph_compact(ListGraph *graph, int *remap)
{
    if (graph == NULL) {
        FAIL(GRAPH_ERR_INVALID_ARGUMENT, false);
    }
    const size_t old_bound = graph->id_bound;
    int *map = remap;
    if (map == NULL && old_bound != 0) {
        map = malloc(old_bound * sizeof(*map));
        if (map == NULL) {
            FAIL(GRAPH_ERR_OUT_OF_MEMORY, false);
        }
    }

    /* Alive vertices keep their relative order, so new ID <= old ID. */
    size_t next = 0;
    for (size_t v = 0; v < old_bound; ++v) {
        map[v] = slot_at(graph, (int)v)->alive ? (int)next++ : -1;
    }

    if (next != old_bound) {
        /* Ascending order makes the in-place move safe: the destination
         * index never exceeds the source, and everything below the source
         * has either been consumed or is dead. */
        for (size_t v = 0; v < old_bound; ++v) {
            if (map[v] >= 0 && (size_t)map[v] != v) {
                *slot_at(graph, map[v]) = *slot_at(graph, (int)v);
            }
        }
        for (size_t v = next; v < old_bound; ++v) {
            memset(slot_at(graph, (int)v), 0, sizeof(Slot));
        }
        for (size_t v = 0; v < next; ++v) {
            Slot *s = slot_at(graph, (int)v);
            for (size_t i = 0; i < s->deg; ++i) {
                s->to[i] = map[s->to[i]];
            }
        }

        const size_t keep = (next + SEG_SIZE - 1) >> LIST_SEG_SHIFT;
        for (size_t s = keep; s < graph->n_segs; ++s) {
            free(graph->segs[s]);
        }
        graph->n_segs = keep;

        free(graph->free_ids);
        graph->free_ids = NULL;
        graph->free_cap = 0;
        graph->n_free   = 0;
        graph->id_bound = next;
        ++graph->version;
    }

    if (remap == NULL) {
        free(map);
    }
    return true;
}

/* ------------------------------------------------------------------ */
/*  Edge modification                                                   */
/* ------------------------------------------------------------------ */

static bool add_edge_impl(ListGraph *graph, int u, int v, double weight,
                          bool check_duplicate)
{
    if (!valid_vertex(graph, u) || !valid_vertex(graph, v)) {
        FAIL(BAD_VERTEX(graph), false);
    }
    if (!isfinite(weight) || (graph->unweighted && weight != 1.0)) {
        FAIL(GRAPH_ERR_INVALID_WEIGHT, false);   /* an unweighted graph cannot remember another weight */
    }
    Slot *su = slot_at(graph, u);
    Slot *sv = slot_at(graph, v);

    size_t pos_u;
    if (find_index(su, v, &pos_u) && check_duplicate) {
        FAIL(GRAPH_ERR_DUPLICATE_EDGE, false);
    }
    size_t pos_v = 0;

    const bool mirror = !graph->directed && u != v;
    if (su->deg == UINT32_MAX || (mirror && sv->deg == UINT32_MAX)) {
        FAIL(GRAPH_ERR_TOO_LARGE, false);
    }
    /* Reserve both ends before touching either, so there is nothing to undo. */
    if (!slot_reserve(graph, su, su->deg + 1) ||
        (mirror && !slot_reserve(graph, sv, sv->deg + 1))) {
        FAIL_ALLOC(false);
    }

    if (mirror) {
        find_index(sv, u, &pos_v);
    }
    insert_entry(graph, su, pos_u, v, weight);
    if (mirror) {
        insert_entry(graph, sv, pos_v, u, weight);
    }
    ++graph->n_edges;
    ++graph->version;
    return true;
}

bool list_graph_add_edge(ListGraph *graph, int u, int v, double weight)
{
    return add_edge_impl(graph, u, v, weight, true);
}

bool list_graph_add_new_edge(ListGraph *graph, int u, int v, double weight)
{
    return add_edge_impl(graph, u, v, weight, false);
}

bool list_graph_remove_edge(ListGraph *graph, int u, int v)
{
    if (!valid_vertex(graph, u) || !valid_vertex(graph, v)) {
        FAIL(BAD_VERTEX(graph), false);
    }
    Slot *su = slot_at(graph, u);
    size_t i;
    if (!find_index(su, v, &i)) {
        FAIL(GRAPH_ERR_EDGE_NOT_FOUND, false);
    }
    erase_entry(graph, su, i);

    if (!graph->directed && u != v) {
        Slot *sv = slot_at(graph, v);
        size_t j;
        if (find_index(sv, u, &j)) { /* always true for a consistent graph */
            erase_entry(graph, sv, j);
        }
    }
    --graph->n_edges;
    ++graph->version;
    return true;
}

bool list_graph_has_edge(const ListGraph *graph, int u, int v, double *weight_out)
{
    if (!valid_vertex(graph, u) || !valid_vertex(graph, v)) {
        return false;
    }
    const Slot *su = slot_at(graph, u);
    size_t i;
    if (!find_index(su, v, &i)) {
        return false;
    }
    if (weight_out != NULL) {
        *weight_out = su->w != NULL ? su->w[i] : 1.0;
    }
    return true;
}

void list_graph_print(const ListGraph *graph)
{
    if (graph == NULL) {
        return;
    }
    for (size_t v = 0; v < graph->id_bound; ++v) {
        const Slot *s = slot_at(graph, (int)v);
        if (!s->alive) {
            continue;
        }
        printf("%zu:", v);
        for (size_t i = 0; i < s->deg; ++i) {
            printf(" ->[%d, %.4g]", s->to[i], s->w != NULL ? s->w[i] : 1.0);
        }
        printf("\n");
    }
}

/* ------------------------------------------------------------------ */
/*  BFS                                                                 */
/* ------------------------------------------------------------------ */

/*
 * The visit order of a BFS is exactly its queue order, so the caller's
 * 'order' array doubles as the queue and only a visited bitset is allocated.
 */
size_t list_bfs(const ListGraph *graph, int start, int *order, int *parent)
{
    if (!valid_vertex(graph, start)) {
        FAIL(BAD_VERTEX(graph), 0);
    }
    if (order == NULL) {
        FAIL(GRAPH_ERR_INVALID_ARGUMENT, 0);
    }

    const size_t n = graph->id_bound;
    uint64_t *visited = calloc((n + 63) >> 6, sizeof(*visited));
    if (visited == NULL) {
        FAIL(GRAPH_ERR_OUT_OF_MEMORY, 0);
    }
    if (parent != NULL) {
        for (size_t i = 0; i < n; ++i) {
            parent[i] = -1;
        }
    }

    size_t head = 0, tail = 0;
    order[tail++] = start;
    bit_set(visited, (size_t)start);

    while (head < tail) {
        const int u = order[head++];
        const Slot *s = slot_at(graph, u);
        for (size_t i = 0; i < s->deg; ++i) {
            const int v = s->to[i];
            if (!bit_get(visited, (size_t)v)) {
                bit_set(visited, (size_t)v);
                if (parent != NULL) {
                    parent[v] = u;
                }
                order[tail++] = v;
            }
        }
    }

    free(visited);
    return tail;
}

/* ------------------------------------------------------------------ */
/*  DFS                                                                 */
/* ------------------------------------------------------------------ */

/*
 * Each stack frame is a vertex plus a cursor into its adjacency array, so the
 * stack never holds more than one entry per vertex (a stack of pending
 * neighbours can reach one entry per edge) and no buffers are built per
 * vertex.
 */
size_t list_dfs(const ListGraph *graph, int start, int *order, int *parent)
{
    if (!valid_vertex(graph, start)) {
        FAIL(BAD_VERTEX(graph), 0);
    }
    if (order == NULL) {
        FAIL(GRAPH_ERR_INVALID_ARGUMENT, 0);
    }

    const size_t n = graph->id_bound;
    uint64_t *visited = calloc((n + 63) >> 6, sizeof(*visited));
    int      *stack   = malloc(n * 2 * sizeof(*stack));
    if (visited == NULL || stack == NULL) {
        free(visited);
        free(stack);
        FAIL(GRAPH_ERR_OUT_OF_MEMORY, 0);
    }
    uint32_t *cursor = (uint32_t *)(stack + n);

    if (parent != NULL) {
        for (size_t i = 0; i < n; ++i) {
            parent[i] = -1;
        }
    }

    size_t top = 0, count = 0;
    stack[top++] = start;
    cursor[start] = 0;
    bit_set(visited, (size_t)start);
    order[count++] = start;

    while (top > 0) {
        const int u = stack[top - 1];
        const Slot *s = slot_at(graph, u);
        if (cursor[u] < s->deg) {
            const int v = s->to[cursor[u]++];
            if (!bit_get(visited, (size_t)v)) {
                bit_set(visited, (size_t)v);
                if (parent != NULL) {
                    parent[v] = u;
                }
                order[count++] = v;
                cursor[v] = 0;
                stack[top++] = v;
            }
        } else {
            --top;
        }
    }

    free(visited);
    free(stack);
    return count;
}

/* ------------------------------------------------------------------ */
/*  Connected components (weakly connected for directed graphs)         */
/* ------------------------------------------------------------------ */

static int uf_find(int *par, int v)
{
    while (par[v] != v) {
        par[v] = par[par[v]]; /* path halving */
        v      = par[v];
    }
    return v;
}

/*
 * Undirected graphs: flood fill.  Directed graphs: a union-find over every
 * edge, which treats edges as undirected without needing a reverse index.
 */
size_t list_connected_components(const ListGraph *graph, int *component)
{
    if (graph == NULL || (graph->id_bound != 0 && component == NULL)) {
        FAIL(GRAPH_ERR_INVALID_ARGUMENT, 0);
    }

    const size_t n = graph->id_bound;
    if (n == 0) {
        return 0;
    }
    int *work = malloc(n * 2 * sizeof(*work));
    if (work == NULL) {
        FAIL(GRAPH_ERR_OUT_OF_MEMORY, 0);
    }

    size_t comp_count = 0;

    if (!graph->directed) {
        int *queue = work;
        for (size_t i = 0; i < n; ++i) {
            component[i] = -1;
        }
        for (size_t start = 0; start < n; ++start) {
            if (component[start] != -1 || !slot_at(graph, (int)start)->alive) {
                continue;
            }
            size_t head = 0, tail = 0;
            queue[tail++]    = (int)start;
            component[start] = (int)comp_count;
            while (head < tail) {
                const Slot *s = slot_at(graph, queue[head++]);
                for (size_t i = 0; i < s->deg; ++i) {
                    const int v = s->to[i];
                    if (component[v] == -1) {
                        component[v] = (int)comp_count;
                        queue[tail++] = v;
                    }
                }
            }
            ++comp_count;
        }
    } else {
        int *par   = work;
        int *label = work + n;
        for (size_t i = 0; i < n; ++i) {
            par[i]   = (int)i;
            label[i] = -1;
        }
        for (size_t u = 0; u < n; ++u) {
            const Slot *s = slot_at(graph, (int)u);
            for (size_t i = 0; i < s->deg; ++i) {
                const int ru = uf_find(par, (int)u);
                const int rv = uf_find(par, s->to[i]);
                if (ru != rv) {
                    par[ru] = rv;
                }
            }
        }
        for (size_t v = 0; v < n; ++v) {
            if (!slot_at(graph, (int)v)->alive) {
                component[v] = -1;
                continue;
            }
            const int root = uf_find(par, (int)v);
            if (label[root] == -1) {
                label[root] = (int)comp_count++;
            }
            component[v] = label[root];
        }
    }

    free(work);
    return comp_count;
}

/* ------------------------------------------------------------------ */
/*  Cycle detection                                                     */
/* ------------------------------------------------------------------ */

/*
 * Iterative DFS with colours (0 unseen, 1 on the stack, 2 finished).
 *
 *   Directed:   reaching a colour-1 vertex is a back edge, i.e. a cycle.
 *   Undirected: reaching any already-seen vertex other than the one we came
 *               from is a non-tree edge.  Parallel edges cannot exist, so
 *               skipping the parent vertex skips exactly the tree edge.
 */
bool list_has_cycle(const ListGraph *graph)
{
    if (graph == NULL) {
        FAIL(GRAPH_ERR_INVALID_ARGUMENT, false);
    }
    if (live_count(graph) == 0) {
        return false;
    }

    const size_t n = graph->id_bound;
    int           *stack  = malloc(n * 2 * sizeof(*stack));
    uint32_t      *cursor = malloc(n * sizeof(*cursor));
    unsigned char *color  = calloc(n, sizeof(*color));
    if (stack == NULL || cursor == NULL || color == NULL) {
        free(stack);
        free(cursor);
        free(color);
        FAIL(GRAPH_ERR_OUT_OF_MEMORY, false);
    }
    int *par = stack + n;

    bool found = false;
    for (size_t root = 0; root < n && !found; ++root) {
        if (color[root] != 0 || !slot_at(graph, (int)root)->alive) {
            continue;
        }
        size_t top = 0;
        stack[top++]  = (int)root;
        color[root]   = 1;
        cursor[root]  = 0;
        par[root]     = -1;

        while (top > 0 && !found) {
            const int u = stack[top - 1];
            const Slot *s = slot_at(graph, u);
            if (cursor[u] >= s->deg) {
                color[u] = 2;
                --top;
                continue;
            }
            const int w = s->to[cursor[u]++];
            if (color[w] == 0) {
                color[w]  = 1;
                cursor[w] = 0;
                par[w]    = u;
                stack[top++] = w;
            } else if (graph->directed ? color[w] == 1 : w != par[u]) {
                found = true;
            }
        }
    }

    free(stack);
    free(cursor);
    free(color);
    return found;
}

/* ------------------------------------------------------------------ */
/*  Topological sort (Kahn's BFS algorithm)                            */
/* ------------------------------------------------------------------ */

bool list_topological_sort(const ListGraph *graph, int *order)
{
    if (graph == NULL) {
        FAIL(GRAPH_ERR_INVALID_ARGUMENT, false);
    }
    if (!graph->directed) {
        FAIL(GRAPH_ERR_UNSUPPORTED, false);
    }
    const size_t live = live_count(graph);
    if (live != 0 && order == NULL) {
        FAIL(GRAPH_ERR_INVALID_ARGUMENT, false);
    }
    if (live == 0) {
        return true;
    }

    const size_t n = graph->id_bound;
    uint32_t *in_degree = calloc(n, sizeof(*in_degree));
    if (in_degree == NULL) {
        FAIL(GRAPH_ERR_OUT_OF_MEMORY, false);
    }

    for (size_t u = 0; u < n; ++u) {
        const Slot *s = slot_at(graph, (int)u);
        for (size_t i = 0; i < s->deg; ++i) {
            ++in_degree[s->to[i]];
        }
    }

    /* The output array is also the work queue. */
    size_t head = 0, tail = 0;
    for (size_t v = 0; v < n; ++v) {
        if (slot_at(graph, (int)v)->alive && in_degree[v] == 0) {
            order[tail++] = (int)v;
        }
    }
    while (head < tail) {
        const Slot *s = slot_at(graph, order[head++]);
        for (size_t i = 0; i < s->deg; ++i) {
            const int v = s->to[i];
            if (--in_degree[v] == 0) {
                order[tail++] = v;
            }
        }
    }

    free(in_degree);
    if (tail != live) {  /* fewer means a cycle blocked some vertices */
        FAIL(GRAPH_ERR_CYCLE, false);
    }
    return true;
}

/* ------------------------------------------------------------------ */
/*  Dijkstra's shortest paths                                           */
/* ------------------------------------------------------------------ */

/*
 * Indexed binary min-heap keyed by distance[].  pos[v] is v's heap index,
 * -1 if it has not been reached yet, -2 once it is settled.  Decrease-key
 * keeps the heap at most V entries, unlike lazy deletion (up to E entries).
 */
static void heap_up(int *heap, int *pos, const double *key, size_t i)
{
    const int v = heap[i];
    while (i > 0) {
        const size_t p = (i - 1) / 2;
        if (key[heap[p]] <= key[v]) {
            break;
        }
        heap[i]       = heap[p];
        pos[heap[i]]  = (int)i;
        i = p;
    }
    heap[i] = v;
    pos[v]  = (int)i;
}

static void heap_down(int *heap, int *pos, const double *key, size_t size, size_t i)
{
    const int v = heap[i];
    for (;;) {
        size_t c = 2 * i + 1;
        if (c >= size) {
            break;
        }
        if (c + 1 < size && key[heap[c + 1]] < key[heap[c]]) {
            ++c;
        }
        if (key[v] <= key[heap[c]]) {
            break;
        }
        heap[i]      = heap[c];
        pos[heap[i]] = (int)i;
        i = c;
    }
    heap[i] = v;
    pos[v]  = (int)i;
}

bool list_dijkstra(const ListGraph *graph, int source, double *distance, int *parent)
{
    if (!valid_vertex(graph, source)) {
        FAIL(BAD_VERTEX(graph), false);
    }
    if (distance == NULL) {
        FAIL(GRAPH_ERR_INVALID_ARGUMENT, false);
    }
    /* Dijkstra requires non-negative weights; the graph counts negatives. */
    if (graph->neg_entries != 0) {
        FAIL(GRAPH_ERR_NEGATIVE_WEIGHT, false);
    }

    const size_t n = graph->id_bound;
    int *heap = malloc(n * 2 * sizeof(*heap));
    if (heap == NULL) {
        FAIL(GRAPH_ERR_OUT_OF_MEMORY, false);
    }
    int *pos = heap + n;

    for (size_t v = 0; v < n; ++v) {
        distance[v] = INFINITY;
        pos[v]      = -1;
        if (parent != NULL) {
            parent[v] = -1;
        }
    }
    distance[source] = 0.0;
    heap[0]          = source;
    pos[source]      = 0;
    size_t size      = 1;

    while (size > 0) {
        const int u = heap[0];
        --size;
        if (size > 0) {
            heap[0] = heap[size];
            heap_down(heap, pos, distance, size, 0);
        }
        pos[u] = -2;

        const Slot *s = slot_at(graph, u);
        for (size_t i = 0; i < s->deg; ++i) {
            const int w = s->to[i];
            if (pos[w] == -2) {
                continue;
            }
            const double candidate = distance[u] + (s->w != NULL ? s->w[i] : 1.0);
            if (candidate < distance[w]) {
                distance[w] = candidate;
                if (parent != NULL) {
                    parent[w] = u;
                }
                if (pos[w] == -1) {
                    heap[size] = w;
                    pos[w]     = (int)size;
                    ++size;
                    heap_up(heap, pos, distance, size - 1);
                } else {
                    heap_up(heap, pos, distance, (size_t)pos[w]);
                }
            }
        }
    }

    free(heap);
    return true;
}

/* ------------------------------------------------------------------ */
/*  Kruskal's MST                                                       */
/* ------------------------------------------------------------------ */

static int compare_list_edges(const void *lhs, const void *rhs)
{
    const ListEdge *a = lhs;
    const ListEdge *b = rhs;
    if (a->weight < b->weight) { return -1; }
    if (a->weight > b->weight) { return  1; }
    if (a->u != b->u)          { return a->u < b->u ? -1 : 1; }
    if (a->v != b->v)          { return a->v < b->v ? -1 : 1; }
    return 0;
}

double list_mst_kruskal(const ListGraph *graph, ListEdge *edges, size_t *edge_count)
{
    if (edge_count == NULL) {
        FAIL(GRAPH_ERR_INVALID_ARGUMENT, INFINITY);
    }
    *edge_count = 0;

    if (graph == NULL) {
        FAIL(GRAPH_ERR_INVALID_ARGUMENT, INFINITY);
    }
    if (graph->directed) {
        FAIL(GRAPH_ERR_UNSUPPORTED, INFINITY);
    }
    const size_t live = live_count(graph);
    if (live > 1 && edges == NULL) {
        FAIL(GRAPH_ERR_INVALID_ARGUMENT, INFINITY);
    }
    if (live <= 1) {
        return 0.0;
    }

    const size_t n = graph->id_bound;
    const size_t e = graph->n_edges;

    /*
     * Undirected edges are stored twice (u->v and v->u); collect only the
     * copy with u < v.  Self-loops never join two components, so skipping
     * them (v == u) is also correct.
     */
    ListEdge *candidates = malloc((e != 0 ? e : 1) * sizeof(*candidates));
    int      *par        = malloc(n * sizeof(*par));
    if (candidates == NULL || par == NULL) {
        free(candidates);
        free(par);
        FAIL(GRAPH_ERR_OUT_OF_MEMORY, INFINITY);
    }

    size_t cnt = 0;
    for (size_t u = 0; u < n; ++u) {
        par[u] = (int)u;
        const Slot *s = slot_at(graph, (int)u);
        for (size_t i = 0; i < s->deg; ++i) {
            if ((size_t)s->to[i] > u) {
                candidates[cnt++] = (ListEdge){ (int)u, s->to[i], s->w != NULL ? s->w[i] : 1.0 };
            }
        }
    }

    qsort(candidates, cnt, sizeof(*candidates), compare_list_edges);

    double total = 0.0;
    for (size_t i = 0; i < cnt && *edge_count + 1 < live; ++i) {
        const int ra = uf_find(par, candidates[i].u);
        const int rb = uf_find(par, candidates[i].v);
        if (ra != rb) {
            par[ra] = rb;
            edges[(*edge_count)++] = candidates[i];
            total += candidates[i].weight;
        }
    }

    free(candidates);
    free(par);

    /* A spanning tree over the alive vertices has exactly live-1 edges. */
    if (*edge_count + 1 != live) {
        FAIL(GRAPH_ERR_NOT_FOUND, INFINITY);   /* graph is disconnected */
    }
    return total;
}
