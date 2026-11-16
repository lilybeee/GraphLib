/*
 * graph_shard.c — the sharding engine; see graph_shard.h for what it is for and why it is built
 * this way.  Read that header first: the decisions are explained there.
 *
 * File layout of one shard (little-endian, native integers; the files are private scratch data
 * that live only as long as the handle, so no portability is needed):
 *
 *     uint64  offsets[count + 1]     row r of the shard starts at neighbours[offsets[r]]
 *     int32   neighbours[entries]    the target of every stored arc, rows back to back
 *
 * where count = the number of users in the shard's ID range [bound[i], bound[i+1]).
 *
 * Building uses three passes over the edge list, each holding only O(users) integers:
 *   1. degrees   -> decides where to cut the ID range into shards;
 *   2. bucket    -> writes every arc (u, v) to the raw file of the shard that owns u;
 *   3. finish    -> turns each raw file into a CSR shard with a counting sort (one shard in
 *                   memory at a time) and deletes the raw file.
 */

#define _POSIX_C_SOURCE 200809L

#include "graph_shard.h"

#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "graph_algorithms.h"
#include "graph_stream.h"
#include "list_backend.h"
#include "memory_budget.h"

#define MAX_SHARDS      256     /* raw files are open together while bucketing; stay under fd limits */
#define MIN_SHARD_ENTRIES 1024
#define BYTES_PER_USER  24      /* what the algorithms need per user: degree + PageRank's two arrays */

typedef struct {
    int                shard;   /* which shard is held, or -1 for a free slot */
    uint64_t          *offsets;
    int32_t           *nbrs;
    unsigned           pins;    /* >0: in use, must not be evicted */
    unsigned long long stamp;   /* LRU clock value of the last use */
} Slot;

struct ShardedGraph {
    char    *dir;
    bool     owns_dir;
    bool     directed;
    size_t   vertices;
    size_t   edges;             /* distinct edges once repeats are merged */
    size_t   loops;             /* self-loops among them (stored once) */
    size_t   nshards;
    int     *bound;             /* nshards + 1 user-ID boundaries */
    size_t  *entries;           /* stored arcs per shard */
    uint32_t *degree;           /* out-degree of every user: the one per-user array kept for good */
    size_t   budget;
    size_t   shard_bytes_max;
    size_t   cache_cap;
    Slot    *slots;
    unsigned long long clock;
    size_t   loads;
    uint64_t bytes_read;
};

ShardConfig shard_default_config(void)
{
    ShardConfig c = {0, NULL, false};
    return c;
}

/* ------------------------------------------------------------------ */
/*  Small helpers                                                       */
/* ------------------------------------------------------------------ */

static void shard_path(const ShardedGraph *sg, size_t i, bool raw, char *out, size_t cap)
{
    snprintf(out, cap, "%s/%s_%06zu.%s", sg->dir, raw ? "raw" : "shard", i, raw ? "tmp" : "bin");
}

/** Which shard owns user v: binary search on the boundaries. */
static size_t shard_of(const ShardedGraph *sg, int v)
{
    size_t lo = 0, hi = sg->nshards;            /* find the last i with bound[i] <= v */
    while (hi - lo > 1) {
        const size_t mid = (lo + hi) / 2;
        if (sg->bound[mid] <= v) lo = mid; else hi = mid;
    }
    return lo;
}

static size_t shard_bytes(size_t users, size_t entries)
{
    return (users + 1) * sizeof(uint64_t) + entries * sizeof(int32_t);
}

static size_t effective_budget(const ShardConfig *cfg)
{
    size_t b = cfg != NULL && cfg->memory_budget != 0 ? cfg->memory_budget : graph_memory_budget();
    return b != 0 ? b : GRAPH_DEFAULT_MEMORY_BUDGET;     /* "unlimited" still needs a shard size */
}

void shard_destroy(ShardedGraph *sg)
{
    if (sg == NULL) {
        return;
    }
    char path[1024];
    for (size_t i = 0; i < sg->nshards; i++) {
        shard_path(sg, i, false, path, sizeof path);
        remove(path);
        shard_path(sg, i, true, path, sizeof path);
        remove(path);                                 /* only exists if a build was interrupted */
    }
    if (sg->owns_dir && sg->dir != NULL) {
        rmdir(sg->dir);
    }
    if (sg->slots != NULL) {
        for (size_t i = 0; i < sg->cache_cap; i++) {
            free(sg->slots[i].offsets);
            free(sg->slots[i].nbrs);
        }
    }
    free(sg->slots);
    free(sg->degree);
    free(sg->bound);
    free(sg->entries);
    free(sg->dir);
    free(sg);
}

/* ------------------------------------------------------------------ */
/*  Building                                                            */
/* ------------------------------------------------------------------ */

typedef struct {
    ShardedGraph *sg;
    uint32_t     *deg;       /* pass 1 */
    FILE        **raw;       /* pass 2: one open raw file per shard */
    bool          failed;
} BuildCtx;

static GraphStatus count_edge(void *ctx, int u, int v)
{
    BuildCtx *b = ctx;
    if ((size_t)u >= b->sg->vertices || (size_t)v >= b->sg->vertices) {
        return GRAPH_ERR_IO;                          /* the file changed between passes */
    }
    b->deg[u]++;
    if (!b->sg->directed && u != v) {
        b->deg[v]++;
    }
    return GRAPH_OK;
}

static bool write_arc(BuildCtx *b, int from, int to)
{
    const int32_t rec[2] = {from, to};
    return fwrite(rec, sizeof rec, 1, b->raw[shard_of(b->sg, from)]) == 1;
}

static GraphStatus bucket_edge(void *ctx, int u, int v)
{
    BuildCtx *b = ctx;
    if ((size_t)u >= b->sg->vertices || (size_t)v >= b->sg->vertices) {
        return GRAPH_ERR_IO;
    }
    if (!write_arc(b, u, v) || (!b->sg->directed && u != v && !write_arc(b, v, u))) {
        return GRAPH_ERR_IO;                          /* disk full or unwritable */
    }
    return GRAPH_OK;
}

/** Cuts [0, n) into shards of about `target` arcs each.  Returns the boundaries (count in *k). */
static int *cut_ranges(const uint32_t *deg, size_t n, size_t target, size_t *k, size_t *max_entries)
{
    size_t cap = 16, count = 0, in_shard = 0;
    int *bound = malloc(cap * sizeof *bound);
    if (bound == NULL) {
        return NULL;
    }
    bound[count++] = 0;
    *max_entries   = 0;
    for (size_t v = 0; v < n; v++) {
        if (in_shard > 0 && in_shard + deg[v] > target) {   /* close the shard before this user */
            if (in_shard > *max_entries) *max_entries = in_shard;
            if (count == cap) {
                int *bigger = realloc(bound, cap * 2 * sizeof *bound);
                if (bigger == NULL) { free(bound); return NULL; }
                bound = bigger;
                cap *= 2;
            }
            bound[count++] = (int)v;
            in_shard = 0;
        }
        in_shard += deg[v];
    }
    if (in_shard > *max_entries) *max_entries = in_shard;
    if (count == cap) {
        int *bigger = realloc(bound, (cap + 1) * sizeof *bound);
        if (bigger == NULL) { free(bound); return NULL; }
        bound = bigger;
    }
    bound[count++] = (int)n;
    *k = count - 1;
    return bound;
}

static int compare_int32(const void *a, const void *b)
{
    const int32_t x = *(const int32_t *)a, y = *(const int32_t *)b;
    return (x > y) - (x < y);
}

static GraphStatus finish_shard(ShardedGraph *sg, size_t i)
{
    char raw_path[1024], out_path[1024];
    shard_path(sg, i, true, raw_path, sizeof raw_path);
    shard_path(sg, i, false, out_path, sizeof out_path);
    const size_t lo = (size_t)sg->bound[i], users = (size_t)sg->bound[i + 1] - lo, E = sg->entries[i];

    GraphStatus st = GRAPH_ERR_OUT_OF_MEMORY;
    int32_t  *rec  = malloc((E > 0 ? E : 1) * 2 * sizeof *rec);
    int32_t  *nbrs = malloc((E > 0 ? E : 1) * sizeof *nbrs);
    uint64_t *off  = calloc(users + 2, sizeof *off);
    FILE *in = NULL, *out = NULL;
    if (rec == NULL || nbrs == NULL || off == NULL) {
        goto done;
    }
    st = GRAPH_ERR_IO;
    in = fopen(raw_path, "rb");
    if (in == NULL || (E > 0 && fread(rec, 2 * sizeof *rec, E, in) != E)) {
        goto done;
    }
    /* Counting sort by source: off[r+1] counts row r, prefix sums give row starts. */
    for (size_t k = 0; k < E; k++) {
        off[(size_t)rec[2 * k] - lo + 1]++;
    }
    for (size_t r = 0; r < users; r++) {
        off[r + 1] += off[r];
    }
    uint64_t *cursor = malloc((users + 1) * sizeof *cursor);
    if (cursor == NULL) {
        st = GRAPH_ERR_OUT_OF_MEMORY;
        goto done;
    }
    memcpy(cursor, off, (users + 1) * sizeof *cursor);
    for (size_t k = 0; k < E; k++) {
        nbrs[cursor[(size_t)rec[2 * k] - lo]++] = rec[2 * k + 1];
    }
    free(cursor);

    /* Merge repeated edges, as graph_add_edge does for an in-memory graph: sort each row, drop
     * equal neighbours, close the gaps.  `kept` becomes the shard's real number of arcs. */
    uint64_t kept = 0;
    uint64_t *final_off = malloc((users + 1) * sizeof *final_off);
    if (final_off == NULL) {
        st = GRAPH_ERR_OUT_OF_MEMORY;
        goto done;
    }
    for (size_t r = 0; r < users; r++) {
        final_off[r] = kept;
        int32_t *row = nbrs + off[r];
        const size_t len = (size_t)(off[r + 1] - off[r]);
        qsort(row, len, sizeof *row, compare_int32);
        for (size_t k = 0; k < len; k++) {
            if (k == 0 || row[k] != row[k - 1]) {
                nbrs[kept++] = row[k];
            }
        }
        sg->degree[lo + r] = (uint32_t)(kept - final_off[r]);
        for (uint64_t k = final_off[r]; k < kept; k++) {
            sg->loops += nbrs[k] == (int32_t)(lo + r);
        }
    }
    final_off[users] = kept;
    sg->entries[i]   = (size_t)kept;
    out = fopen(out_path, "wb");
    if (out == NULL || fwrite(final_off, sizeof *final_off, users + 1, out) != users + 1 ||
        (kept > 0 && fwrite(nbrs, sizeof *nbrs, (size_t)kept, out) != (size_t)kept)) {
        free(final_off);
        goto done;
    }
    free(final_off);
    st = GRAPH_OK;
done:
    if (in != NULL)  fclose(in);
    if (out != NULL && fclose(out) != 0 && st == GRAPH_OK) st = GRAPH_ERR_IO;
    remove(raw_path);
    free(rec); free(nbrs); free(off);
    return st;
}

GraphStatus shard_build(FILE *in, const ShardConfig *config, ShardedGraph **out)
{
    if (in == NULL || out == NULL) {
        return graph_set_error(GRAPH_ERR_INVALID_ARGUMENT);
    }
    *out = NULL;
    rewind(in);
    GraphStreamInfo info;
    GraphStatus st = graph_stream_scan(in, &info);
    if (st != GRAPH_OK) {
        return st;
    }
    const size_t n = info.vertices, budget = effective_budget(config);
    if (n >= INT_MAX) {
        return graph_set_error(GRAPH_ERR_TOO_LARGE);
    }
    /* The per-user arrays are the one thing that must always fit. */
    if (n * BYTES_PER_USER > budget / 2) {
        return graph_set_error(GRAPH_ERR_TOO_LARGE);
    }

    ShardedGraph *sg = calloc(1, sizeof *sg);
    if (sg == NULL) {
        return graph_set_error(GRAPH_ERR_OUT_OF_MEMORY);
    }
    sg->vertices = n;
    sg->edges    = info.edges;
    sg->directed = info.directed;
    sg->budget   = budget;
    if (config != NULL && config->dir != NULL) {
        sg->dir = strdup(config->dir);
    } else {
        char tmpl[] = "/tmp/graph_shards_XXXXXX";
        sg->dir      = mkdtemp(tmpl) != NULL ? strdup(tmpl) : NULL;
        sg->owns_dir = sg->dir != NULL;
    }
    sg->degree = calloc(n > 0 ? n : 1, sizeof *sg->degree);
    if (sg->dir == NULL || sg->degree == NULL) {
        st = sg->dir == NULL ? GRAPH_ERR_IO : GRAPH_ERR_OUT_OF_MEMORY;
        goto fail;
    }

    /* Pass 1: degrees. */
    BuildCtx b = {sg, sg->degree, NULL, false};
    rewind(in);
    if ((st = graph_stream_edges(in, count_edge, &b)) != GRAPH_OK) {
        goto fail;
    }
    size_t total = 0;
    for (size_t v = 0; v < n; v++) {
        total += sg->degree[v];
        if ((size_t)sg->degree[v] * 12 > budget / 2) {      /* one user's arcs must fit as a shard */
            st = GRAPH_ERR_TOO_LARGE;
            goto fail;
        }
    }

    /* Cut the ID range.  A shard is built from 12 bytes per arc (8 raw + 4 final), and the cache
     * should hold a few of them, so aim for a quarter of half the budget; fewer, larger shards
     * if that would make more than MAX_SHARDS files. */
    size_t target = (budget / 2 / 4) / 12;
    if (target < MIN_SHARD_ENTRIES) target = MIN_SHARD_ENTRIES;
    size_t max_entries = 0;
    sg->bound = cut_ranges(sg->degree, n, target, &sg->nshards, &max_entries);
    if (sg->bound != NULL && sg->nshards > MAX_SHARDS) {
        free(sg->bound);
        target    = (total + MAX_SHARDS - 1) / MAX_SHARDS;
        sg->bound = cut_ranges(sg->degree, n, target > 0 ? target : 1, &sg->nshards, &max_entries);
    }
    if (sg->bound == NULL) {
        st = GRAPH_ERR_OUT_OF_MEMORY;
        goto fail;
    }
    if (max_entries * 12 > budget / 2) {
        st = GRAPH_ERR_TOO_LARGE;
        goto fail;
    }
    sg->entries = calloc(sg->nshards > 0 ? sg->nshards : 1, sizeof *sg->entries);
    if (sg->entries == NULL) {
        st = GRAPH_ERR_OUT_OF_MEMORY;
        goto fail;
    }
    for (size_t i = 0; i < sg->nshards; i++) {
        size_t users = (size_t)(sg->bound[i + 1] - sg->bound[i]);
        for (int v = sg->bound[i]; v < sg->bound[i + 1]; v++) sg->entries[i] += sg->degree[v];
        if (shard_bytes(users, sg->entries[i]) > sg->shard_bytes_max) {
            sg->shard_bytes_max = shard_bytes(users, sg->entries[i]);
        }
    }

    /* Pass 2: bucket the arcs into one raw file per shard. */
    FILE **raw   = calloc(sg->nshards > 0 ? sg->nshards : 1, sizeof *raw);
    char **bufs  = calloc(sg->nshards > 0 ? sg->nshards : 1, sizeof *bufs);
    size_t bufsz = sg->nshards > 0 ? budget / 8 / sg->nshards : 4096;
    bufsz        = bufsz < 4096 ? 4096 : (bufsz > 262144 ? 262144 : bufsz);
    if (raw == NULL || bufs == NULL) {
        free(raw); free(bufs);
        st = GRAPH_ERR_OUT_OF_MEMORY;
        goto fail;
    }
    b.raw = raw;
    st = GRAPH_OK;
    for (size_t i = 0; i < sg->nshards && st == GRAPH_OK; i++) {
        char path[1024];
        shard_path(sg, i, true, path, sizeof path);
        bufs[i] = malloc(bufsz);
        raw[i]  = fopen(path, "wb");
        if (bufs[i] == NULL || raw[i] == NULL) {
            st = raw[i] == NULL ? GRAPH_ERR_IO : GRAPH_ERR_OUT_OF_MEMORY;
        } else {
            setvbuf(raw[i], bufs[i], _IOFBF, bufsz);
        }
    }
    if (st == GRAPH_OK) {
        rewind(in);
        st = graph_stream_edges(in, bucket_edge, &b);
    }
    for (size_t i = 0; i < sg->nshards; i++) {
        if (raw[i] != NULL && fclose(raw[i]) != 0 && st == GRAPH_OK) st = GRAPH_ERR_IO;
        free(bufs[i]);
    }
    free(raw); free(bufs);
    if (st != GRAPH_OK) {
        goto fail;
    }

    /* Pass 3: raw -> CSR, one shard in memory at a time. */
    for (size_t i = 0; i < sg->nshards; i++) {
        if ((st = finish_shard(sg, i)) != GRAPH_OK) {
            goto fail;
        }
    }

    /* Repeats are gone now: the real edge count (an undirected edge is stored as two arcs, a
     * self-loop as one) and the real size of the biggest shard. */
    size_t arcs = 0;
    sg->shard_bytes_max = 0;
    for (size_t i = 0; i < sg->nshards; i++) {
        arcs += sg->entries[i];
        const size_t users = (size_t)(sg->bound[i + 1] - sg->bound[i]);
        if (shard_bytes(users, sg->entries[i]) > sg->shard_bytes_max) {
            sg->shard_bytes_max = shard_bytes(users, sg->entries[i]);
        }
    }
    sg->edges = sg->directed ? arcs : (arcs + sg->loops) / 2;

    /* The cache: as many shards as half the budget allows, but at least two so a sweep can
     * move on while the previous shard is still resident, and never more than exist. */
    sg->cache_cap = sg->shard_bytes_max > 0 ? (budget / 2) / sg->shard_bytes_max : 1;
    if (sg->cache_cap < 2) sg->cache_cap = 2;
    if (sg->cache_cap > sg->nshards) sg->cache_cap = sg->nshards;
    sg->slots = calloc(sg->cache_cap > 0 ? sg->cache_cap : 1, sizeof *sg->slots);
    if (sg->slots == NULL) {
        st = GRAPH_ERR_OUT_OF_MEMORY;
        goto fail;
    }
    for (size_t i = 0; i < sg->cache_cap; i++) {
        sg->slots[i].shard = -1;
    }
    *out = sg;
    return GRAPH_OK;

fail:
    shard_destroy(sg);
    return graph_set_error(st);
}

GraphStatus graph_open_edge_list(FILE *in, const ShardConfig *config, Graph **graph, ShardedGraph **sharded)
{
    if (in == NULL || graph == NULL || sharded == NULL) {
        return graph_set_error(GRAPH_ERR_INVALID_ARGUMENT);
    }
    *graph   = NULL;
    *sharded = NULL;
    rewind(in);
    GraphStreamInfo info;
    GraphStatus st = graph_stream_scan(in, &info);
    if (st != GRAPH_OK) {
        return st;
    }
    /* Would the ordinary in-memory list fit?  Ask the list's own cost model. */
    const GraphShape shape = {info.vertices, info.vertices, info.edges, info.directed, false};
    const size_t need = list_backend()->estimate_memory(&shape);
    if (!(config != NULL && config->force) && need <= effective_budget(config)) {
        rewind(in);
        return graph_read_edge_list(in, info.directed, GRAPH_LIST, graph);
    }
    return shard_build(in, config, sharded);
}

/* ------------------------------------------------------------------ */
/*  The handle                                                          */
/* ------------------------------------------------------------------ */

size_t shard_vertex_count(const ShardedGraph *sg) { return sg == NULL ? 0 : sg->vertices; }
size_t shard_edge_count(const ShardedGraph *sg)   { return sg == NULL ? 0 : sg->edges; }
bool   shard_is_directed(const ShardedGraph *sg)  { return sg != NULL && sg->directed; }

size_t shard_degree(const ShardedGraph *sg, int user)
{
    return sg != NULL && user >= 0 && (size_t)user < sg->vertices ? sg->degree[user] : 0;
}

ShardStats shard_stats(const ShardedGraph *sg)
{
    ShardStats s = {0, 0, 0, 0, 0};
    if (sg != NULL) {
        s.shards = sg->nshards;
        s.cache_capacity = sg->cache_cap;
        s.shard_bytes_max = sg->shard_bytes_max;
        s.loads = sg->loads;
        s.bytes_read = sg->bytes_read;
    }
    return s;
}

/* ------------------------------------------------------------------ */
/*  The cache: pin / unpin                                              */
/* ------------------------------------------------------------------ */

/** Makes shard i resident and pinned (so it cannot be evicted until unpinned). */
static GraphStatus pin(ShardedGraph *sg, size_t i, Slot **out)
{
    Slot *hit = NULL, *victim = NULL;
    for (size_t k = 0; k < sg->cache_cap; k++) {
        Slot *s = &sg->slots[k];
        if (s->shard == (int)i) {
            hit = s;
            break;
        }
        if (s->pins == 0 && (victim == NULL || s->shard == -1 ||
                             (victim->shard != -1 && s->stamp < victim->stamp))) {
            victim = s;
        }
    }
    if (hit == NULL) {
        if (victim == NULL) {
            return GRAPH_ERR_INVALID_ARGUMENT;        /* everything pinned: a sweep pins one at a time */
        }
        free(victim->offsets);                        /* read-only data: evicting is just freeing */
        free(victim->nbrs);
        victim->offsets = NULL;
        victim->nbrs    = NULL;
        victim->shard   = -1;

        char path[1024];
        shard_path(sg, i, false, path, sizeof path);
        const size_t users = (size_t)(sg->bound[i + 1] - sg->bound[i]), E = sg->entries[i];
        victim->offsets = malloc((users + 1) * sizeof *victim->offsets);
        victim->nbrs    = malloc((E > 0 ? E : 1) * sizeof *victim->nbrs);
        FILE *f         = fopen(path, "rb");
        if (victim->offsets == NULL || victim->nbrs == NULL || f == NULL ||
            fread(victim->offsets, sizeof *victim->offsets, users + 1, f) != users + 1 ||
            (E > 0 && fread(victim->nbrs, sizeof *victim->nbrs, E, f) != E)) {
            const GraphStatus why = f == NULL ? GRAPH_ERR_IO
                                    : (victim->offsets == NULL || victim->nbrs == NULL)
                                          ? GRAPH_ERR_OUT_OF_MEMORY : GRAPH_ERR_IO;
            if (f != NULL) fclose(f);
            free(victim->offsets); free(victim->nbrs);
            victim->offsets = NULL; victim->nbrs = NULL;
            return why;
        }
        fclose(f);
        victim->shard = (int)i;
        sg->loads++;
        sg->bytes_read += shard_bytes(users, E);
        hit = victim;
    }
    hit->pins++;
    hit->stamp = ++sg->clock;
    *out = hit;
    return GRAPH_OK;
}

static void unpin(Slot *s)
{
    if (s->pins > 0) {
        s->pins--;
    }
}

/* ------------------------------------------------------------------ */
/*  Algorithms                                                          */
/* ------------------------------------------------------------------ */

GraphStatus shard_degrees(const ShardedGraph *sg, size_t *out_degree)
{
    if (sg == NULL || (sg->vertices > 0 && out_degree == NULL)) {
        return graph_set_error(GRAPH_ERR_INVALID_ARGUMENT);
    }
    for (size_t v = 0; v < sg->vertices; v++) {
        out_degree[v] = sg->degree[v];                 /* counted while the shards were cut */
    }
    return GRAPH_OK;
}

static int find_root(int *parent, int x)
{
    while (parent[x] != x) {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    return x;
}

GraphStatus shard_components(ShardedGraph *sg, int *component, size_t *count)
{
    if (sg == NULL || (sg->vertices > 0 && component == NULL)) {
        return graph_set_error(GRAPH_ERR_INVALID_ARGUMENT);
    }
    const size_t n = sg->vertices;
    int *parent = malloc((n > 0 ? n : 1) * sizeof *parent), *size = malloc((n > 0 ? n : 1) * sizeof *size);
    if (parent == NULL || size == NULL) {
        free(parent); free(size);
        return graph_set_error(GRAPH_ERR_OUT_OF_MEMORY);
    }
    for (size_t v = 0; v < n; v++) { parent[v] = (int)v; size[v] = 1; }

    GraphStatus st = GRAPH_OK;
    for (size_t i = 0; i < sg->nshards && st == GRAPH_OK; i++) {          /* one sweep */
        Slot *s;
        if ((st = pin(sg, i, &s)) != GRAPH_OK) break;
        const int lo = sg->bound[i];
        for (int u = lo; u < sg->bound[i + 1]; u++) {
            for (uint64_t k = s->offsets[u - lo]; k < s->offsets[u - lo + 1]; k++) {
                int a = find_root(parent, u), b = find_root(parent, s->nbrs[k]);
                if (a != b) {
                    if (size[a] < size[b]) { int t = a; a = b; b = t; }
                    parent[b] = a;
                    size[a] += size[b];
                }
            }
        }
        unpin(s);
    }
    if (st == GRAPH_OK) {
        for (size_t v = 0; v < n; v++) size[v] = -1;
        size_t labels = 0;
        for (size_t v = 0; v < n; v++) {
            const int r = find_root(parent, (int)v);
            if (size[r] < 0) size[r] = (int)labels++;
            component[v] = size[r];
        }
        if (count != NULL) *count = labels;
    }
    free(parent); free(size);
    return st == GRAPH_OK ? GRAPH_OK : graph_set_error(st);
}

GraphStatus shard_pagerank(ShardedGraph *sg, double damping, double tolerance, size_t max_iterations,
                           double *rank, size_t *iterations)
{
    if (sg == NULL || (sg->vertices > 0 && rank == NULL) || damping < 0.0 || damping >= 1.0 ||
        tolerance < 0.0 || max_iterations == 0) {
        return graph_set_error(GRAPH_ERR_INVALID_ARGUMENT);
    }
    const size_t n = sg->vertices;
    if (iterations != NULL) *iterations = 0;
    if (n == 0) {
        return GRAPH_OK;
    }
    double *next = malloc(n * sizeof *next);
    if (next == NULL) {
        return graph_set_error(GRAPH_ERR_OUT_OF_MEMORY);
    }
    for (size_t v = 0; v < n; v++) rank[v] = 1.0 / (double)n;

    GraphStatus st = GRAPH_OK;
    bool converged = false;
    size_t round = 0;
    while (round < max_iterations && st == GRAPH_OK && !converged) {
        double dangling = 0.0;
        for (size_t v = 0; v < n; v++) {
            next[v] = 0.0;
            if (sg->degree[v] == 0) dangling += rank[v];       /* no friends: spreads its rank evenly */
        }
        /* One sweep per round.  An LRU cache smaller than the graph would miss on every shard of a
         * sweep that always runs 0..k-1 (the oldest is evicted just before it is needed again), so
         * alternate the direction: each round then starts on the shards the last one finished with,
         * which are still resident. */
        for (size_t step = 0; step < sg->nshards && st == GRAPH_OK; step++) {
            const size_t i = (round % 2 == 0) ? step : sg->nshards - 1 - step;
            Slot *s;
            if ((st = pin(sg, i, &s)) != GRAPH_OK) break;
            const int lo = sg->bound[i];
            for (int u = lo; u < sg->bound[i + 1]; u++) {
                const uint64_t from = s->offsets[u - lo], to = s->offsets[u - lo + 1];
                if (to == from) continue;
                const double share = rank[u] / (double)(to - from);
                for (uint64_t k = from; k < to; k++) next[s->nbrs[k]] += share;
            }
            unpin(s);
        }
        if (st != GRAPH_OK) break;
        double diff = 0.0;
        for (size_t v = 0; v < n; v++) {
            const double r = (1.0 - damping) / (double)n + damping * (next[v] + dangling / (double)n);
            diff += fabs(r - rank[v]);
            rank[v] = r;
        }
        round++;
        converged = diff < tolerance;
    }
    if (iterations != NULL) *iterations = round;
    free(next);
    if (st != GRAPH_OK) {
        return graph_set_error(st);
    }
    return converged ? GRAPH_OK : graph_set_error(GRAPH_ERR_NO_CONVERGENCE);
}

GraphStatus shard_bfs(ShardedGraph *sg, int source, size_t max_hops, int *level, size_t *reached)
{
    if (sg == NULL || (sg->vertices > 0 && level == NULL)) {
        return graph_set_error(GRAPH_ERR_INVALID_ARGUMENT);
    }
    const size_t n = sg->vertices;
    if (source < 0 || (size_t)source >= n) {
        return graph_set_error(GRAPH_ERR_INVALID_VERTEX);
    }
    unsigned char *flag = calloc(sg->nshards, 1), *flag_next = calloc(sg->nshards, 1);
    if (flag == NULL || flag_next == NULL) {
        free(flag); free(flag_next);
        return graph_set_error(GRAPH_ERR_OUT_OF_MEMORY);
    }
    for (size_t v = 0; v < n; v++) level[v] = -1;
    level[source] = 0;
    size_t count = 1;
    flag[shard_of(sg, source)] = 1;

    GraphStatus st = GRAPH_OK;
    for (size_t cur = 0; cur < max_hops && st == GRAPH_OK; cur++) {
        bool any = false;
        memset(flag_next, 0, sg->nshards);
        for (size_t i = 0; i < sg->nshards && st == GRAPH_OK; i++) {
            if (!flag[i]) continue;                     /* no user of this level lives here: skip the read */
            any = true;
            Slot *s;
            if ((st = pin(sg, i, &s)) != GRAPH_OK) break;
            const int lo = sg->bound[i];
            for (int u = lo; u < sg->bound[i + 1]; u++) {
                if (level[u] != (int)cur) continue;
                for (uint64_t k = s->offsets[u - lo]; k < s->offsets[u - lo + 1]; k++) {
                    const int v = s->nbrs[k];
                    if (level[v] < 0) {
                        level[v] = (int)cur + 1;
                        count++;
                        flag_next[shard_of(sg, v)] = 1;
                    }
                }
            }
            unpin(s);
        }
        if (!any) break;
        unsigned char *t = flag; flag = flag_next; flag_next = t;
    }
    free(flag); free(flag_next);
    if (reached != NULL) *reached = count;
    return st == GRAPH_OK ? GRAPH_OK : graph_set_error(st);
}
