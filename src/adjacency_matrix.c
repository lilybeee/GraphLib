/*
 * adjacency_matrix.c — Adjacency-matrix graph backend
 *
 * Implements everything declared in adjacency_matrix.h.
 *
 * Internal representation
 * -----------------------
 *   Presence bits.  Rows are grouped in strips of 64; strips[row >> 6] points
 *   to 64 rows of `dim` words each, so
 *
 *     word (row, j) = strips[row >> 6][(row & 63) * dim + j]
 *
 *   holds the edges from `row` to columns j*64 .. j*64+63.  A row is one
 *   contiguous run of memory, which is what makes the scans below fast: BFS
 *   ANDs a row's words with the complement of a visited bitset and gets the
 *   newly discovered neighbours 64 at a time.  The visited bitsets used by the
 *   algorithms have the same layout (word j covers columns j*64 .. j*64+63).
 *   A strip is NULL until an edge leaves one of its rows; NULL reads as all
 *   zeros.
 *
 *   Weights.  A weight tile is 64 x 64 doubles and exists only where an edge has landed in
 *   that block; the tiles are found through a sparse hash index keyed by (tile row, tile
 *   column), so the index costs memory in proportion to the tiles that exist -- the dense
 *   (V/64)^2 pointer directory it replaced was 2 GB at a million vertices before a single
 *   edge was added.  Tiles never move, so growing the matrix copies only the bit strips.  A
 *   weight is only read where its presence bit is set, so tiles are never initialised.  An
 *   unweighted graph (matrix_graph_create_unweighted) has no tiles at all.
 *
 *   Size limit.  Once every row strip exists the presence bits take V^2/8 bytes, so
 *   matrix_max_vertices() bounds V by the memory budget (about 131 000 weighted, 176 000
 *   unweighted under 4 GiB) and a larger matrix is refused outright.
 *
 *   Vertex bookkeeping (alive bits, out-degrees, the free list of dead IDs) is
 *   the same as in adjacency_list.c.
 *
 * Compile: cc -std=c11 -Wall -Wextra -Werror -pedantic adjacency_matrix.c -lm
 */

#include "adjacency_matrix.h"
#include "graph_status.h"
#include "memory_budget.h"

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TILE        64
#define TILE_SHIFT  6
#define TILE_MASK   63

typedef struct Tile {
    double weights[TILE * TILE];   /* weights[r * 64 + c]; valid iff the bit is set */
} Tile;

/* Sparse index of the weight tiles: open addressing, key = (tile row, tile column) packed
 * into one word (+1 so that 0 means "empty slot").  Memory is proportional to the number of
 * tiles that exist, not to (V/64)^2 as the dense directory it replaces was. */
typedef struct {
    uint64_t *key;
    Tile    **val;
    size_t    cap;        /* a power of two, or 0 before the first tile */
    size_t    count;
} TileMap;

struct MatrixGraph {
    uint64_t **strips;    /* dim pointers; NULL = every row of the strip is empty */
    TileMap   tiles;      /* weight tiles that exist; empty for an unweighted graph */
    size_t    dim;        /* words per row; the ID capacity is dim * 64 */
    uint64_t *alive;      /* dim words; bit v set iff v is an alive vertex */
    uint32_t *degree;     /* dim * 64 out-degrees */
    size_t    id_bound;   /* IDs 0 .. id_bound-1 exist (alive or dead) */
    int      *free_ids;   /* dead IDs available for reuse */
    size_t    n_free;
    size_t    free_cap;
    size_t    edges;      /* undirected edges counted once */
    size_t    neg_entries;/* stored cells with weight < 0 */
    size_t    n_strips;   /* bit strips currently allocated */
    size_t    budget;     /* bytes this graph may use (0 = unlimited); see memory_budget.h */
    uint64_t  version;
    bool      directed;
    bool      unweighted; /* no weight tiles: every edge reads back as 1.0 */
};

/** A failed allocation is OUT_OF_MEMORY unless the budget check already said TOO_LARGE. */
#define FAIL_ALLOC(value) do {                                         \
        if (graph_last_error() != GRAPH_ERR_TOO_LARGE) {               \
            graph_set_error(GRAPH_ERR_OUT_OF_MEMORY);                  \
        }                                                              \
        return (value);                                                \
    } while (0)

/** Would growing by 'extra' bytes exceed the graph's memory budget?  Records TOO_LARGE if so. */
static bool over_budget(const MatrixGraph *graph, size_t extra);

/* ------------------------------------------------------------------ */
/*  Bit helpers                                                         */
/* ------------------------------------------------------------------ */

static inline int ctz64(uint64_t x) {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_ctzll(x);
#else
    int n = 0;
    while ((x & 1) == 0) {
        x >>= 1;
        ++n;
    }
    return n;
#endif
}

static inline bool bit_get(const uint64_t *bits, size_t i) {
    return (bits[i >> 6] >> (i & 63)) & 1u;
}

static inline void bit_set(uint64_t *bits, size_t i) {
    bits[i >> 6] |= (uint64_t)1 << (i & 63);
}

/* ------------------------------------------------------------------ */
/*  Internal helpers                                                    */
/* ------------------------------------------------------------------ */

static inline size_t strip_words(const MatrixGraph *graph) {
    return TILE * graph->dim;
}

/** The `dim` presence words of a row, or NULL if the row's strip is empty. */
static inline uint64_t *row_bits(const MatrixGraph *graph, size_t row) {
    uint64_t *strip = graph->strips[row >> TILE_SHIFT];
    return strip == NULL ? NULL : strip + (row & TILE_MASK) * graph->dim;
}

static inline uint64_t tile_key(size_t tile_row, size_t tile_col) {
    return (((uint64_t)tile_row << 32) | (uint64_t)tile_col) + 1;
}

static inline size_t tile_hash(uint64_t key, size_t mask) {
    key *= 0x9E3779B97F4A7C15ULL;
    return (size_t)(key >> 32) & mask;
}

/** The tile at (tile_row, tile_col), or NULL if no edge has landed in that block. */
static Tile *tile_at(const MatrixGraph *graph, size_t tile_row, size_t tile_col) {
    const TileMap *m = &graph->tiles;
    if (m->cap == 0) {
        return NULL;
    }
    const uint64_t key = tile_key(tile_row, tile_col);
    for (size_t i = tile_hash(key, m->cap - 1);; i = (i + 1) & (m->cap - 1)) {
        if (m->key[i] == key) {
            return m->val[i];
        }
        if (m->key[i] == 0) {
            return NULL;
        }
    }
}

/** Adds a tile to the index (the key must not be present).  False if out of memory. */
static bool tilemap_insert(TileMap *m, uint64_t key, Tile *tile) {
    if ((m->count + 1) * 2 > m->cap) {                       /* keep the load under one half */
        const size_t cap = m->cap != 0 ? m->cap * 2 : 64;
        uint64_t *nkey = calloc(cap, sizeof(*nkey));
        Tile    **nval = malloc(cap * sizeof(*nval));
        if (nkey == NULL || nval == NULL) {
            free(nkey);
            free(nval);
            return false;
        }
        for (size_t i = 0; i < m->cap; ++i) {
            if (m->key[i] != 0) {
                size_t j = tile_hash(m->key[i], cap - 1);
                while (nkey[j] != 0) {
                    j = (j + 1) & (cap - 1);
                }
                nkey[j] = m->key[i];
                nval[j] = m->val[i];
            }
        }
        free(m->key);
        free(m->val);
        m->key = nkey;
        m->val = nval;
        m->cap = cap;
    }
    size_t i = tile_hash(key, m->cap - 1);
    while (m->key[i] != 0) {
        i = (i + 1) & (m->cap - 1);
    }
    m->key[i] = key;
    m->val[i] = tile;
    ++m->count;
    return true;
}

/** The 64 weights of one row inside a tile.  An unweighted graph has no tiles: all ones. */
static const double ONES[TILE] = {
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

static inline const double *row_weights(const MatrixGraph *graph, size_t tile_row, size_t tile_col,
                                        size_t row_in_tile) {
    if (graph->unweighted) {
        return ONES;
    }
    return tile_at(graph, tile_row, tile_col)->weights + row_in_tile * TILE;
}

/* Record a failure code, then return the given value. */
#define FAIL(code, value) do { graph_set_error(code); return (value); } while (0)

static bool valid_vertex(const MatrixGraph *graph, int v) {
    return graph != NULL && v >= 0 && (size_t)v < graph->id_bound &&
           bit_get(graph->alive, (size_t)v);
}

static inline size_t live_count(const MatrixGraph *graph) {
    return graph->id_bound - graph->n_free;
}

static bool has_cell(const MatrixGraph *graph, size_t row, size_t col) {
    const uint64_t *bits = row_bits(graph, row);
    return bits != NULL && bit_get(bits, col);
}

/** Weight of a cell whose presence bit is set. */
static inline double weight_at(const MatrixGraph *graph, size_t row, size_t col) {
    if (graph->unweighted) {
        return 1.0;
    }
    const Tile *tile = tile_at(graph, row >> TILE_SHIFT, col >> TILE_SHIFT);
    return tile->weights[(row & TILE_MASK) * TILE + (col & TILE_MASK)];
}

/** Makes sure the strip and the weight tile that cover (row, col) exist. */
static bool ensure_cell_storage(MatrixGraph *graph, size_t row, size_t col) {
    uint64_t **strip = &graph->strips[row >> TILE_SHIFT];
    if (*strip == NULL) {
        if (over_budget(graph, strip_words(graph) * sizeof(**strip))) {
            return false;
        }
        *strip = calloc(strip_words(graph), sizeof(**strip));
        if (*strip == NULL) {
            return false;
        }
        ++graph->n_strips;
    }
    if (graph->unweighted) {
        return true;                          /* the presence bit is all there is */
    }
    if (tile_at(graph, row >> TILE_SHIFT, col >> TILE_SHIFT) == NULL) {
        /* the tile plus its index slot (key + pointer, at most doubled by the load factor) */
        if (over_budget(graph, sizeof(Tile) + 4 * sizeof(uint64_t))) {
            return false;
        }
        Tile *tile = malloc(sizeof(*tile));
        if (tile == NULL) {
            return false;
        }
        if (!tilemap_insert(&graph->tiles, tile_key(row >> TILE_SHIFT, col >> TILE_SHIFT), tile)) {
            free(tile);
            return false;
        }
    }
    return true;
}

/** Sets a cell whose storage already exists. */
static void set_cell(MatrixGraph *graph, size_t row, size_t col, double weight) {
    bit_set(row_bits(graph, row), col);
    if (!graph->unweighted) {
        tile_at(graph, row >> TILE_SHIFT, col >> TILE_SHIFT)
            ->weights[(row & TILE_MASK) * TILE + (col & TILE_MASK)] = weight;
    }
}

/** Clears a cell that is known to be set. */
static void clear_cell(MatrixGraph *graph, size_t row, size_t col) {
    row_bits(graph, row)[col >> 6] &= ~((uint64_t)1 << (col & 63));
}

/**
 * next_set - smallest column >= from that is set in row u, or SIZE_MAX.
 * If 'exclude' is non-NULL, columns whose bit is set there are skipped (a
 * visited bitset in the row layout).
 */
static size_t next_set(const MatrixGraph *graph, size_t u, size_t from,
                       const uint64_t *exclude) {
    const uint64_t *bits = row_bits(graph, u);
    if (bits == NULL) {
        return SIZE_MAX;
    }
    size_t word_index = from >> 6;
    uint64_t mask = ~(uint64_t)0 << (from & 63);

    for (; word_index < graph->dim; ++word_index, mask = ~(uint64_t)0) {
        uint64_t word = bits[word_index] & mask;
        if (exclude != NULL) {
            word &= ~exclude[word_index];
        }
        if (word != 0) {
            return (word_index << 6) + (size_t)ctz64(word);
        }
    }
    return SIZE_MAX;
}

/**
 * grow_capacity - widen the matrix to new_dim words per row.  The weight tiles
 * are not touched; the pointer directory and the (small) bit strips are
 * re-laid-out.  Leaves the graph unchanged on failure.
 */
static bool grow_capacity(MatrixGraph *graph, size_t new_dim) {
    if (new_dim <= graph->dim) {
        return true;
    }
    if (new_dim > SIZE_MAX / new_dim || new_dim > SIZE_MAX / TILE ||
        new_dim * TILE > (size_t)INT_MAX ||
        new_dim > SIZE_MAX / (TILE * TILE * sizeof(uint64_t))) {
        return false;
    }
    /* The calculated limit: every row strip allocated is V^2/8 bytes, so a matrix over this many
     * vertices could not be filled within the budget however few edges it starts with. */
    if (new_dim * TILE > matrix_max_vertices(graph->budget, !graph->unweighted)) {
        graph_set_error(GRAPH_ERR_TOO_LARGE);
        return false;
    }

    uint64_t **strips = calloc(new_dim, sizeof(*strips));
    uint64_t *alive = calloc(new_dim, sizeof(*alive));
    uint32_t *degree = calloc(new_dim * TILE, sizeof(*degree));
    bool ok = strips != NULL && alive != NULL && degree != NULL;

    /* Re-lay-out every allocated strip at the new row width. */
    for (size_t s = 0; ok && s < graph->dim; ++s) {
        if (graph->strips[s] == NULL) {
            continue;
        }
        strips[s] = calloc(TILE * new_dim, sizeof(**strips));
        ok = strips[s] != NULL;
        for (size_t r = 0; ok && r < TILE; ++r) {
            memcpy(strips[s] + r * new_dim, graph->strips[s] + r * graph->dim,
                   graph->dim * sizeof(**strips));
        }
    }
    if (!ok) {
        for (size_t s = 0; strips != NULL && s < new_dim; ++s) {
            free(strips[s]);
        }
        free(strips);
        free(alive);
        free(degree);
        return false;
    }

    if (graph->dim != 0) {
        memcpy(alive, graph->alive, graph->dim * sizeof(*alive));
        memcpy(degree, graph->degree, graph->dim * TILE * sizeof(*degree));
    }
    for (size_t s = 0; s < graph->dim; ++s) {
        free(graph->strips[s]);
    }
    free(graph->strips);
    free(graph->alive);
    free(graph->degree);
    graph->strips = strips;
    graph->alive = alive;
    graph->degree = degree;
    graph->dim = new_dim;
    return true;
}

/** Most vertices a matrix may hold under @budget bytes (0 = unlimited). */
size_t matrix_max_vertices(size_t budget, bool weighted) {
    if (budget == 0) {
        return (size_t)INT_MAX;
    }
    /* Presence bits alone are V^2/8 bytes once every strip exists.  Weighted graphs keep half the
     * budget back for weight tiles, unweighted ones the tenth for the per-vertex arrays. */
    const double share = weighted ? 0.5 : 0.9;
    const double n = floor(sqrt(8.0 * share * (double)budget));
    return n >= (double)INT_MAX ? (size_t)INT_MAX : (size_t)n;
}

static void release_storage(MatrixGraph *graph) {
    for (size_t i = 0; i < graph->tiles.cap; ++i) {
        if (graph->tiles.key[i] != 0) {
            free(graph->tiles.val[i]);
        }
    }
    free(graph->tiles.key);
    free(graph->tiles.val);
    for (size_t s = 0; s < graph->dim; ++s) {
        free(graph->strips[s]);
    }
    free(graph->strips);
    free(graph->alive);
    free(graph->degree);
    free(graph->free_ids);
}

/* ------------------------------------------------------------------ */
/*  Construction / destruction                                          */
/* ------------------------------------------------------------------ */

static MatrixGraph *create_graph(size_t vertices, bool directed, bool unweighted) {
    if (vertices > (size_t)INT_MAX) {
        return NULL;
    }
    MatrixGraph *graph = calloc(1, sizeof(*graph));
    if (graph == NULL) {
        return NULL;
    }
    graph->directed   = directed;
    graph->unweighted = unweighted;
    graph->budget     = graph_memory_budget();

    if (vertices != 0 && !grow_capacity(graph, (vertices + TILE - 1) >> TILE_SHIFT)) {
        free(graph);
        if (graph_last_error() != GRAPH_ERR_TOO_LARGE) {
            graph_set_error(GRAPH_ERR_OUT_OF_MEMORY);
        }
        return NULL;
    }
    for (size_t v = 0; v < vertices; ++v) {
        bit_set(graph->alive, v);
    }
    graph->id_bound = vertices;
    return graph;
}

MatrixGraph *matrix_graph_create(size_t vertices, bool directed) {
    return create_graph(vertices, directed, false);
}

MatrixGraph *matrix_graph_create_unweighted(size_t vertices, bool directed) {
    return create_graph(vertices, directed, true);
}

bool matrix_graph_is_unweighted(const MatrixGraph *graph) {
    return graph != NULL && graph->unweighted;
}

void matrix_graph_destroy(MatrixGraph *graph) {
    if (graph == NULL) {
        return;
    }
    release_storage(graph);
    free(graph);
}

/* ------------------------------------------------------------------ */
/*  Queries                                                             */
/* ------------------------------------------------------------------ */

size_t matrix_graph_vertex_count(const MatrixGraph *graph) {
    return graph == NULL ? 0 : live_count(graph);
}

size_t matrix_graph_id_bound(const MatrixGraph *graph) {
    return graph == NULL ? 0 : graph->id_bound;
}

bool matrix_graph_is_alive(const MatrixGraph *graph, int v) {
    return valid_vertex(graph, v);
}

size_t matrix_graph_edge_count(const MatrixGraph *graph) {
    return graph == NULL ? 0 : graph->edges;
}

bool matrix_graph_is_directed(const MatrixGraph *graph) {
    return graph != NULL && graph->directed;
}

size_t matrix_graph_degree(const MatrixGraph *graph, int v) {
    return valid_vertex(graph, v) ? graph->degree[v] : 0;
}

bool matrix_graph_next_neighbor(const MatrixGraph *graph, int v, size_t *cursor,
                                int *to, double *weight) {
    if (!valid_vertex(graph, v) || cursor == NULL || to == NULL) {
        return false;
    }
    const size_t column = next_set(graph, (size_t)v, *cursor, NULL);
    if (column == SIZE_MAX) {
        return false;
    }
    *to = (int)column;
    if (weight != NULL) {
        *weight = weight_at(graph, (size_t)v, column);
    }
    *cursor = column + 1;
    return true;
}

uint64_t matrix_graph_version(const MatrixGraph *graph) {
    return graph == NULL ? 0 : graph->version;
}

size_t matrix_graph_memory_bytes(const MatrixGraph *graph) {
    if (graph == NULL) {
        return 0;
    }
    return sizeof(*graph) + graph->tiles.cap * (sizeof(uint64_t) + sizeof(Tile *)) +
           graph->dim * sizeof(uint64_t *) + graph->dim * sizeof(uint64_t) +
           graph->dim * TILE * sizeof(uint32_t) + graph->free_cap * sizeof(int) +
           graph->n_strips * strip_words(graph) * sizeof(uint64_t) +
           graph->tiles.count * sizeof(Tile);
}

static bool over_budget(const MatrixGraph *graph, size_t extra) {
    if (graph->budget != 0 && matrix_graph_memory_bytes(graph) + extra > graph->budget) {
        graph_set_error(GRAPH_ERR_TOO_LARGE);
        return true;
    }
    return false;
}

/* ------------------------------------------------------------------ */
/*  Vertex modification                                                 */
/* ------------------------------------------------------------------ */

int matrix_graph_add_vertex(MatrixGraph *graph) {
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
        if (graph->id_bound == graph->dim * TILE) {
            /* Double, but never past the size limit: grow to the limit itself if that is nearer. */
            size_t want = graph->dim != 0 ? graph->dim * 2 : 1;
            const size_t limit_dim = matrix_max_vertices(graph->budget, !graph->unweighted) / TILE;
            if (want > limit_dim) {
                want = limit_dim;
            }
            if (want <= graph->dim) {
                FAIL(GRAPH_ERR_TOO_LARGE, -1);       /* already at the limit */
            }
            if (!grow_capacity(graph, want)) {
                FAIL_ALLOC(-1);
            }
        }
        id = (int)graph->id_bound++;
    }
    bit_set(graph->alive, (size_t)id);
    ++graph->version;
    return id;
}

bool matrix_graph_remove_vertex(MatrixGraph *graph, int v) {
    if (!valid_vertex(graph, v)) {
        FAIL(graph == NULL ? GRAPH_ERR_INVALID_ARGUMENT : GRAPH_ERR_INVALID_VERTEX, false);
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

    const size_t vertex = (size_t)v;
    const size_t vertex_word = vertex >> 6;
    const uint64_t vertex_bit = (uint64_t)1 << (vertex & 63);
    uint64_t *bits = row_bits(graph, vertex);

    /* Row v: every edge leaving v.  In an undirected graph each of them has a
     * mirror cell in column v of the neighbour's row, cleared here too. */
    if (bits != NULL) {
        for (size_t j = 0; j < graph->dim; ++j) {
            uint64_t word = bits[j];
            while (word != 0) {
                const size_t column = (j << 6) + (size_t)ctz64(word);
                word &= word - 1;
                const bool negative = weight_at(graph, vertex, column) < 0.0;

                if (!graph->directed && column != vertex) {
                    clear_cell(graph, column, vertex);
                    --graph->degree[column];
                    if (negative) {
                        --graph->neg_entries;   /* the mirror cell */
                    }
                }
                if (negative) {
                    --graph->neg_entries;
                }
                --graph->edges;
            }
            bits[j] = 0;
        }
    }
    graph->degree[vertex] = 0;

    /* Column v: every edge entering v.  Undirected graphs already handled it
     * through the mirrors; a directed self-loop was in the row. */
    if (graph->directed) {
        for (size_t s = 0; s < graph->dim; ++s) {
            uint64_t *strip = graph->strips[s];
            if (strip == NULL) {
                continue;
            }
            for (size_t r = 0; r < TILE; ++r) {
                uint64_t *word = &strip[r * graph->dim + vertex_word];
                if ((*word & vertex_bit) == 0) {
                    continue;
                }
                const size_t source = (s << TILE_SHIFT) + r;
                *word &= ~vertex_bit;
                --graph->degree[source];
                --graph->edges;
                if (!graph->unweighted &&
                    tile_at(graph, s, vertex_word)->weights[r * TILE + (vertex & TILE_MASK)] < 0.0) {
                    --graph->neg_entries;
                }
            }
        }
    }

    graph->alive[vertex_word] &= ~vertex_bit;
    graph->free_ids[graph->n_free++] = v;
    ++graph->version;
    return true;
}

bool matrix_graph_needs_compaction(const MatrixGraph *graph) {
    return graph != NULL && graph->n_free != 0 && graph->n_free * 4 > graph->id_bound;
}

bool matrix_graph_compact(MatrixGraph *graph, int *remap) {
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
        map[v] = bit_get(graph->alive, v) ? (int)next++ : -1;
    }

    bool ok = true;
    if (next != old_bound) {
        /* Build the compacted matrix beside the old one, then swap. */
        MatrixGraph *fresh = graph->unweighted ? matrix_graph_create_unweighted(next, graph->directed)
                                               : matrix_graph_create(next, graph->directed);
        ok = fresh != NULL;
        for (size_t u = 0; ok && u < old_bound; ++u) {
            if (map[u] < 0) {
                continue;
            }
            const size_t new_row = (size_t)map[u];
            size_t cursor = 0;
            int to;
            double weight;
            while (ok && matrix_graph_next_neighbor(graph, (int)u, &cursor, &to, &weight)) {
                const size_t new_col = (size_t)map[to];
                ok = ensure_cell_storage(fresh, new_row, new_col);
                if (ok) {
                    set_cell(fresh, new_row, new_col, weight);
                }
            }
            if (ok) {
                fresh->degree[new_row] = graph->degree[u];
            }
        }
        if (ok) {
            fresh->edges = graph->edges;
            fresh->neg_entries = graph->neg_entries;
            fresh->version = graph->version + 1;
            release_storage(graph);
            *graph = *fresh;
            free(fresh);
        } else {
            matrix_graph_destroy(fresh);
        }
    }

    if (remap == NULL) {
        free(map);
    }
    return ok;
}

/* ------------------------------------------------------------------ */
/*  Edge modification                                                   */
/* ------------------------------------------------------------------ */

bool matrix_graph_add_edge(MatrixGraph *graph, int source, int destination, double weight) {
    if (!valid_vertex(graph, source) || !valid_vertex(graph, destination)) {
        FAIL(graph == NULL ? GRAPH_ERR_INVALID_ARGUMENT : GRAPH_ERR_INVALID_VERTEX, false);
    }
    if (!isfinite(weight) || (graph->unweighted && weight != 1.0)) {   /* nowhere to keep another weight */
        FAIL(GRAPH_ERR_INVALID_WEIGHT, false);
    }
    const size_t u = (size_t)source;
    const size_t v = (size_t)destination;
    if (has_cell(graph, u, v)) {
        FAIL(GRAPH_ERR_DUPLICATE_EDGE, false);
    }

    /* Allocate the storage for both cells before touching either, so there is
     * nothing to undo if the second allocation fails. */
    const bool mirror = !graph->directed && u != v;
    if (!ensure_cell_storage(graph, u, v) || (mirror && !ensure_cell_storage(graph, v, u))) {
        FAIL_ALLOC(false);
    }

    set_cell(graph, u, v, weight);
    ++graph->degree[u];
    if (mirror) {
        set_cell(graph, v, u, weight);
        ++graph->degree[v];
    }
    if (weight < 0.0) {
        graph->neg_entries += mirror ? 2 : 1;
    }
    ++graph->edges;
    ++graph->version;
    return true;
}

bool matrix_graph_remove_edge(MatrixGraph *graph, int source, int destination) {
    if (!valid_vertex(graph, source) || !valid_vertex(graph, destination)) {
        FAIL(graph == NULL ? GRAPH_ERR_INVALID_ARGUMENT : GRAPH_ERR_INVALID_VERTEX, false);
    }
    const size_t u = (size_t)source;
    const size_t v = (size_t)destination;
    if (!has_cell(graph, u, v)) {
        FAIL(GRAPH_ERR_EDGE_NOT_FOUND, false);
    }

    const bool negative = weight_at(graph, u, v) < 0.0;
    const bool mirror = !graph->directed && u != v;
    clear_cell(graph, u, v);
    --graph->degree[u];
    if (mirror) {
        clear_cell(graph, v, u);
        --graph->degree[v];
    }
    if (negative) {
        graph->neg_entries -= mirror ? 2 : 1;
    }
    --graph->edges;
    ++graph->version;
    return true;
}

bool matrix_graph_has_edge(const MatrixGraph *graph, int source, int destination,
                           double *weight_out) {
    if (!valid_vertex(graph, source) || !valid_vertex(graph, destination)) {
        return false;
    }
    if (!has_cell(graph, (size_t)source, (size_t)destination)) {
        return false;
    }
    if (weight_out != NULL) {
        *weight_out = weight_at(graph, (size_t)source, (size_t)destination);
    }
    return true;
}

void matrix_graph_print(const MatrixGraph *graph) {
    if (graph == NULL) {
        return;
    }
    for (size_t row = 0; row < graph->id_bound; ++row) {
        if (!bit_get(graph->alive, row)) {
            continue;
        }
        bool first = true;
        for (size_t column = 0; column < graph->id_bound; ++column) {
            if (!bit_get(graph->alive, column)) {
                continue;
            }
            if (!first) {
                putchar('\t');
            }
            first = false;
            if (has_cell(graph, row, column)) {
                printf("%g", weight_at(graph, row, column));
            } else {
                printf(".");
            }
        }
        printf("\n");
    }
}

/* ------------------------------------------------------------------ */
/*  BFS                                                                 */
/* ------------------------------------------------------------------ */

/*
 * The visit order of a BFS is its queue order, so the caller's 'order' array
 * doubles as the queue.  Expanding a vertex ANDs each word of its row with the
 * complement of the visited word: whatever survives is exactly the set of newly
 * discovered neighbours, 64 columns per operation.
 */
size_t matrix_bfs(const MatrixGraph *graph, int start, int *order, int *parent) {
    if (!valid_vertex(graph, start)) {
        FAIL(graph == NULL ? GRAPH_ERR_INVALID_ARGUMENT : GRAPH_ERR_INVALID_VERTEX, 0);
    }
    if (order == NULL) {
        FAIL(GRAPH_ERR_INVALID_ARGUMENT, 0);
    }

    uint64_t *visited = calloc(graph->dim, sizeof(*visited));
    if (visited == NULL) {
        FAIL(GRAPH_ERR_OUT_OF_MEMORY, 0);
    }
    if (parent != NULL) {
        for (size_t i = 0; i < graph->id_bound; ++i) {
            parent[i] = -1;
        }
    }

    size_t head = 0;
    size_t tail = 0;
    order[tail++] = start;
    bit_set(visited, (size_t)start);

    while (head < tail) {
        const int vertex = order[head++];
        const uint64_t *bits = row_bits(graph, (size_t)vertex);
        if (bits == NULL) {
            continue;
        }
        for (size_t j = 0; j < graph->dim; ++j) {
            uint64_t fresh = bits[j] & ~visited[j];
            if (fresh == 0) {
                continue;
            }
            visited[j] |= fresh;
            while (fresh != 0) {
                const int neighbor = (int)((j << 6) + (size_t)ctz64(fresh));
                fresh &= fresh - 1;
                if (parent != NULL) {
                    parent[neighbor] = vertex;
                }
                order[tail++] = neighbor;
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
 * Each stack frame is a vertex plus a cursor (the next column to examine in
 * its row).  A vertex is marked when it is first reached, so parent[] is a real
 * DFS tree.  next_set skips already-visited columns 64 at a time, so a finished
 * neighbourhood costs O(V / 64), not O(degree).
 */
size_t matrix_dfs(const MatrixGraph *graph, int start, int *order, int *parent) {
    if (!valid_vertex(graph, start)) {
        FAIL(graph == NULL ? GRAPH_ERR_INVALID_ARGUMENT : GRAPH_ERR_INVALID_VERTEX, 0);
    }
    if (order == NULL) {
        FAIL(GRAPH_ERR_INVALID_ARGUMENT, 0);
    }

    const size_t n = graph->id_bound;
    uint64_t *visited = calloc(graph->dim, sizeof(*visited));
    int *stack = malloc(n * 2 * sizeof(*stack));
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

    size_t top = 0;
    size_t count = 0;
    stack[top++] = start;
    cursor[start] = 0;
    bit_set(visited, (size_t)start);
    order[count++] = start;

    while (top > 0) {
        const int vertex = stack[top - 1];
        const size_t column = next_set(graph, (size_t)vertex, cursor[vertex], visited);
        if (column == SIZE_MAX) {
            --top;
            continue;
        }
        cursor[vertex] = (uint32_t)(column + 1);
        bit_set(visited, column);
        if (parent != NULL) {
            parent[column] = vertex;
        }
        order[count++] = (int)column;
        cursor[column] = 0;
        stack[top++] = (int)column;
    }

    free(visited);
    free(stack);
    return count;
}

/* ------------------------------------------------------------------ */
/*  Connected components (weakly connected for directed graphs)         */
/* ------------------------------------------------------------------ */

static int uf_find(int *par, int v) {
    while (par[v] != v) {
        par[v] = par[par[v]];   /* path halving */
        v = par[v];
    }
    return v;
}

/*
 * Undirected graphs: flood fill with the same word-at-a-time expansion as BFS.
 * Directed graphs: a union-find over every stored edge, which treats edges as
 * undirected without having to look at columns.
 */
size_t matrix_connected_components(const MatrixGraph *graph, int *component) {
    if (graph == NULL || (graph->id_bound != 0 && component == NULL)) {
        FAIL(GRAPH_ERR_INVALID_ARGUMENT, 0);
    }
    const size_t n = graph->id_bound;
    if (n == 0) {
        return 0;
    }

    size_t comp_count = 0;

    if (!graph->directed) {
        uint64_t *visited = calloc(graph->dim, sizeof(*visited));
        int *queue = malloc(n * sizeof(*queue));
        if (visited == NULL || queue == NULL) {
            free(visited);
            free(queue);
            FAIL(GRAPH_ERR_OUT_OF_MEMORY, 0);
        }
        for (size_t i = 0; i < n; ++i) {
            component[i] = -1;
        }
        for (size_t start = 0; start < n; ++start) {
            if (!bit_get(graph->alive, start) || bit_get(visited, start)) {
                continue;
            }
            size_t head = 0;
            size_t tail = 0;
            queue[tail++] = (int)start;
            bit_set(visited, start);
            component[start] = (int)comp_count;
            while (head < tail) {
                const uint64_t *bits = row_bits(graph, (size_t)queue[head++]);
                if (bits == NULL) {
                    continue;
                }
                for (size_t j = 0; j < graph->dim; ++j) {
                    uint64_t fresh = bits[j] & ~visited[j];
                    if (fresh == 0) {
                        continue;
                    }
                    visited[j] |= fresh;
                    while (fresh != 0) {
                        const int neighbor = (int)((j << 6) + (size_t)ctz64(fresh));
                        fresh &= fresh - 1;
                        component[neighbor] = (int)comp_count;
                        queue[tail++] = neighbor;
                    }
                }
            }
            ++comp_count;
        }
        free(visited);
        free(queue);
        return comp_count;
    }

    int *work = malloc(n * 2 * sizeof(*work));
    if (work == NULL) {
        FAIL(GRAPH_ERR_OUT_OF_MEMORY, 0);
    }
    int *par = work;
    int *label = work + n;
    for (size_t i = 0; i < n; ++i) {
        par[i] = (int)i;
        label[i] = -1;
    }
    for (size_t u = 0; u < n; ++u) {
        const uint64_t *bits = row_bits(graph, u);
        if (bits == NULL) {
            continue;
        }
        for (size_t j = 0; j < graph->dim; ++j) {
            uint64_t word = bits[j];
            while (word != 0) {
                const int column = (int)((j << 6) + (size_t)ctz64(word));
                word &= word - 1;
                const int ru = uf_find(par, (int)u);
                const int rv = uf_find(par, column);
                if (ru != rv) {
                    par[ru] = rv;
                }
            }
        }
    }
    for (size_t v = 0; v < n; ++v) {
        if (!bit_get(graph->alive, v)) {
            component[v] = -1;
            continue;
        }
        const int root = uf_find(par, (int)v);
        if (label[root] == -1) {
            label[root] = (int)comp_count++;
        }
        component[v] = label[root];
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
bool matrix_has_cycle(const MatrixGraph *graph) {
    if (graph == NULL) {
        FAIL(GRAPH_ERR_INVALID_ARGUMENT, false);
    }
    if (live_count(graph) == 0) {
        return false;
    }

    const size_t n = graph->id_bound;
    int *stack = malloc(n * 2 * sizeof(*stack));
    uint32_t *cursor = malloc(n * sizeof(*cursor));
    unsigned char *color = calloc(n, sizeof(*color));
    if (stack == NULL || cursor == NULL || color == NULL) {
        free(stack);
        free(cursor);
        free(color);
        FAIL(GRAPH_ERR_OUT_OF_MEMORY, false);
    }
    int *par = stack + n;

    bool found = false;
    for (size_t root = 0; root < n && !found; ++root) {
        if (color[root] != 0 || !bit_get(graph->alive, root)) {
            continue;
        }
        size_t top = 0;
        stack[top++] = (int)root;
        color[root] = 1;
        cursor[root] = 0;
        par[root] = -1;

        while (top > 0 && !found) {
            const int vertex = stack[top - 1];
            const size_t column = next_set(graph, (size_t)vertex, cursor[vertex], NULL);
            if (column == SIZE_MAX) {
                color[vertex] = 2;
                --top;
                continue;
            }
            cursor[vertex] = (uint32_t)(column + 1);
            if (color[column] == 0) {
                color[column] = 1;
                cursor[column] = 0;
                par[column] = vertex;
                stack[top++] = (int)column;
            } else if (graph->directed ? color[column] == 1 : (int)column != par[vertex]) {
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

bool matrix_topological_sort(const MatrixGraph *graph, int *order) {
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

    /* In-degrees: one pass over every set bit. */
    for (size_t u = 0; u < n; ++u) {
        const uint64_t *bits = row_bits(graph, u);
        if (bits == NULL) {
            continue;
        }
        for (size_t j = 0; j < graph->dim; ++j) {
            uint64_t word = bits[j];
            while (word != 0) {
                ++in_degree[(j << 6) + (size_t)ctz64(word)];
                word &= word - 1;
            }
        }
    }

    /* The output array is also the work queue. */
    size_t head = 0;
    size_t tail = 0;
    for (size_t v = 0; v < n; ++v) {
        if (bit_get(graph->alive, v) && in_degree[v] == 0) {
            order[tail++] = (int)v;
        }
    }
    while (head < tail) {
        const uint64_t *bits = row_bits(graph, (size_t)order[head++]);
        if (bits == NULL) {
            continue;
        }
        for (size_t j = 0; j < graph->dim; ++j) {
            uint64_t word = bits[j];
            while (word != 0) {
                const size_t column = (j << 6) + (size_t)ctz64(word);
                word &= word - 1;
                if (--in_degree[column] == 0) {
                    order[tail++] = (int)column;
                }
            }
        }
    }

    free(in_degree);
    if (tail != live) {    /* fewer means a cycle blocked some vertices */
        FAIL(GRAPH_ERR_CYCLE, false);
    }
    return true;
}

/* ------------------------------------------------------------------ */
/*  Dijkstra's shortest paths                                           */
/* ------------------------------------------------------------------ */

/*
 * Indexed binary min-heap keyed by distance[].  pos[v] is v's heap index,
 * -1 if it has not been reached yet, -2 once it is settled.
 */
static void heap_up(int *heap, int *pos, const double *key, size_t i) {
    const int v = heap[i];
    while (i > 0) {
        const size_t p = (i - 1) / 2;
        if (key[heap[p]] <= key[v]) {
            break;
        }
        heap[i] = heap[p];
        pos[heap[i]] = (int)i;
        i = p;
    }
    heap[i] = v;
    pos[v] = (int)i;
}

static void heap_down(int *heap, int *pos, const double *key, size_t size, size_t i) {
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
        heap[i] = heap[c];
        pos[heap[i]] = (int)i;
        i = c;
    }
    heap[i] = v;
    pos[v] = (int)i;
}

bool matrix_dijkstra(const MatrixGraph *graph, int source, double *distance, int *parent) {
    if (!valid_vertex(graph, source)) {
        FAIL(graph == NULL ? GRAPH_ERR_INVALID_ARGUMENT : GRAPH_ERR_INVALID_VERTEX, false);
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
        pos[v] = -1;
        if (parent != NULL) {
            parent[v] = -1;
        }
    }
    distance[source] = 0.0;
    heap[0] = source;
    pos[source] = 0;
    size_t size = 1;

    while (size > 0) {
        const int vertex = heap[0];
        --size;
        if (size > 0) {
            heap[0] = heap[size];
            heap_down(heap, pos, distance, size, 0);
        }
        pos[vertex] = -2;

        const uint64_t *bits = row_bits(graph, (size_t)vertex);
        if (bits == NULL) {
            continue;
        }
        const size_t tile_row = (size_t)vertex >> TILE_SHIFT;
        const double *weights_row = NULL;   /* this row's weights in the current tile */
        for (size_t j = 0; j < graph->dim; ++j) {
            uint64_t word = bits[j];
            if (word == 0) {
                continue;
            }
            weights_row = row_weights(graph, tile_row, j, (size_t)vertex & TILE_MASK);
            while (word != 0) {
                const unsigned offset = (unsigned)ctz64(word);
                word &= word - 1;
                const size_t column = (j << 6) + offset;
                if (pos[column] == -2) {
                    continue;
                }
                const double candidate = distance[vertex] + weights_row[offset];
                if (candidate < distance[column]) {
                    distance[column] = candidate;
                    if (parent != NULL) {
                        parent[column] = vertex;
                    }
                    if (pos[column] == -1) {
                        heap[size] = (int)column;
                        pos[column] = (int)size;
                        ++size;
                        heap_up(heap, pos, distance, size - 1);
                    } else {
                        heap_up(heap, pos, distance, (size_t)pos[column]);
                    }
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

static int compare_edges(const void *left, const void *right) {
    const MatrixEdge *a = left;
    const MatrixEdge *b = right;
    if (a->weight < b->weight) {
        return -1;
    }
    if (a->weight > b->weight) {
        return 1;
    }
    if (a->u != b->u) {
        return a->u < b->u ? -1 : 1;
    }
    if (a->v != b->v) {
        return a->v < b->v ? -1 : 1;
    }
    return 0;
}

double matrix_mst_kruskal(const MatrixGraph *graph, MatrixEdge *edges, size_t *edge_count) {
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
    MatrixEdge *candidates = malloc((graph->edges != 0 ? graph->edges : 1) * sizeof(*candidates));
    int *par = malloc(n * sizeof(*par));
    if (candidates == NULL || par == NULL) {
        free(candidates);
        free(par);
        FAIL(GRAPH_ERR_OUT_OF_MEMORY, INFINITY);
    }
    for (size_t v = 0; v < n; ++v) {
        par[v] = (int)v;
    }

    /*
     * Undirected edges are stored in both directions; collect only the cells
     * above the diagonal.  Self-loops never join two components, so skipping
     * them is also correct.
     */
    size_t count = 0;
    for (size_t u = 0; u < n; ++u) {
        const uint64_t *bits = row_bits(graph, u);
        if (bits == NULL) {
            continue;
        }
        const size_t first_word = u >> 6;
        for (size_t j = first_word; j < graph->dim; ++j) {
            uint64_t word = bits[j];
            if (j == first_word) {
                word = ((u & 63) == 63) ? 0 : (word & (~(uint64_t)0 << ((u & 63) + 1)));
            }
            if (word == 0) {
                continue;
            }
            const double *weights_row = row_weights(graph, u >> TILE_SHIFT, j, u & TILE_MASK);
            while (word != 0) {
                const unsigned offset = (unsigned)ctz64(word);
                word &= word - 1;
                candidates[count++] =
                    (MatrixEdge){(int)u, (int)((j << 6) + offset), weights_row[offset]};
            }
        }
    }
    qsort(candidates, count, sizeof(*candidates), compare_edges);

    double total = 0.0;
    for (size_t i = 0; i < count && *edge_count + 1 < live; ++i) {
        const int left = uf_find(par, candidates[i].u);
        const int right = uf_find(par, candidates[i].v);
        if (left != right) {
            par[left] = right;
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
