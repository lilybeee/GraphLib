/*
 * cache.c — versioned LRU cache.  See cache.h.
 *
 * Layout: a chained hash table for lookup, and a doubly linked list threaded
 * through the same entries for recency (head = most recently used).
 */

#include "cache.h"

#include <stdlib.h>
#include <string.h>

typedef struct Entry {
    CacheKey      key;
    uint64_t      hash;
    void         *data;
    size_t        size;
    struct Entry *prev;   /* towards the head (more recent) */
    struct Entry *next;   /* towards the tail (less recent) */
    struct Entry *chain;  /* next entry in the same hash bucket */
} Entry;

struct Cache {
    Entry   **buckets;
    size_t    n_buckets;  /* power of two */
    Entry    *head;
    Entry    *tail;
    size_t    count;
    size_t    bytes;
    size_t    max_entries;
    size_t    max_bytes;
    uint64_t  version;
    bool      have_version;
    CacheStats stats;
};

#define INITIAL_BUCKETS 16

/* ------------------------------------------------------------------ */
/*  Hashing                                                             */
/* ------------------------------------------------------------------ */

static uint64_t mix64(uint64_t x)
{
    x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27; x *= 0x94d049bb133111ebULL;
    x ^= x >> 31;
    return x;
}

static uint64_t hash_key(CacheKey k)
{
    uint64_t h = mix64(k.algo);
    h = mix64(h ^ (uint64_t)(uint32_t)k.a);
    h = mix64(h ^ ((uint64_t)(uint32_t)k.b << 1));
    return h;
}

static bool key_equal(CacheKey x, CacheKey y)
{
    return x.algo == y.algo && x.a == y.a && x.b == y.b;
}

/* ------------------------------------------------------------------ */
/*  List and table plumbing                                             */
/* ------------------------------------------------------------------ */

static void lru_unlink(Cache *c, Entry *e)
{
    if (e->prev != NULL) { e->prev->next = e->next; } else { c->head = e->next; }
    if (e->next != NULL) { e->next->prev = e->prev; } else { c->tail = e->prev; }
    e->prev = e->next = NULL;
}

static void lru_push_front(Cache *c, Entry *e)
{
    e->prev = NULL;
    e->next = c->head;
    if (c->head != NULL) { c->head->prev = e; } else { c->tail = e; }
    c->head = e;
}

static Entry *table_find(const Cache *c, CacheKey key, uint64_t hash)
{
    for (Entry *e = c->buckets[hash & (c->n_buckets - 1)]; e != NULL; e = e->chain) {
        if (e->hash == hash && key_equal(e->key, key)) {
            return e;
        }
    }
    return NULL;
}

static void table_unlink(Cache *c, Entry *e)
{
    Entry **pp = &c->buckets[e->hash & (c->n_buckets - 1)];
    while (*pp != e) {
        pp = &(*pp)->chain;
    }
    *pp = e->chain;
}

/** Doubles the bucket array once the load factor passes 1. Failure is harmless. */
static void table_grow(Cache *c)
{
    const size_t n = c->n_buckets * 2;
    Entry **buckets = calloc(n, sizeof(*buckets));
    if (buckets == NULL) {
        return;
    }
    for (Entry *e = c->head; e != NULL; e = e->next) {
        Entry **slot = &buckets[e->hash & (n - 1)];
        e->chain = *slot;
        *slot = e;
    }
    free(c->buckets);
    c->buckets   = buckets;
    c->n_buckets = n;
}

/** Remove e from the cache and free it. */
static void entry_drop(Cache *c, Entry *e)
{
    lru_unlink(c, e);
    table_unlink(c, e);
    c->bytes -= e->size;
    --c->count;
    free(e->data);
    free(e);
}

/** Drops everything if the graph changed since the cache last looked. */
static void sync_version(Cache *c, uint64_t version)
{
    if (c->have_version && c->version == version) {
        return;
    }
    while (c->head != NULL) {
        entry_drop(c, c->head);
        ++c->stats.invalidations;
    }
    c->version      = version;
    c->have_version = true;
}

/* ------------------------------------------------------------------ */
/*  Public API                                                          */
/* ------------------------------------------------------------------ */

Cache *cache_create(size_t max_entries, size_t max_bytes)
{
    Cache *c = calloc(1, sizeof(*c));
    if (c == NULL) {
        return NULL;
    }
    c->buckets = calloc(INITIAL_BUCKETS, sizeof(*c->buckets));
    if (c->buckets == NULL) {
        free(c);
        return NULL;
    }
    c->n_buckets  = INITIAL_BUCKETS;
    c->max_entries = max_entries;
    c->max_bytes   = max_bytes;
    return c;
}

void cache_destroy(Cache *cache)
{
    if (cache == NULL) {
        return;
    }
    while (cache->head != NULL) {
        entry_drop(cache, cache->head);
    }
    free(cache->buckets);
    free(cache);
}

const void *cache_lookup(Cache *cache, CacheKey key, uint64_t graph_version,
                         size_t *size)
{
    if (cache == NULL) {
        return NULL;
    }
    sync_version(cache, graph_version);

    Entry *e = table_find(cache, key, hash_key(key));
    if (e == NULL) {
        ++cache->stats.misses;
        return NULL;
    }
    ++cache->stats.hits;
    if (e != cache->head) {
        lru_unlink(cache, e);
        lru_push_front(cache, e);
    }
    if (size != NULL) {
        *size = e->size;
    }
    return e->data;
}

bool cache_store(Cache *cache, CacheKey key, uint64_t graph_version,
                 const void *data, size_t size)
{
    if (cache == NULL || cache->max_entries == 0 || size > cache->max_bytes ||
        (data == NULL && size != 0)) {
        return false;
    }
    sync_version(cache, graph_version);

    const uint64_t hash = hash_key(key);
    Entry *e = table_find(cache, key, hash);

    void *copy = malloc(size != 0 ? size : 1);
    if (copy == NULL) {
        return false;
    }
    if (size != 0) {
        memcpy(copy, data, size);
    }

    if (e != NULL) {                       /* replace in place */
        free(e->data);
        cache->bytes = cache->bytes - e->size + size;
        e->data = copy;
        e->size = size;
        lru_unlink(cache, e);
        lru_push_front(cache, e);
    } else {
        e = calloc(1, sizeof(*e));
        if (e == NULL) {
            free(copy);
            return false;
        }
        e->key  = key;
        e->hash = hash;
        e->data = copy;
        e->size = size;
        Entry **slot = &cache->buckets[hash & (cache->n_buckets - 1)];
        e->chain = *slot;
        *slot = e;
        lru_push_front(cache, e);
        cache->bytes += size;
        if (++cache->count > cache->n_buckets) {
            table_grow(cache);
        }
    }
    ++cache->stats.stores;

    /* e is the head, so it is evicted last; size <= max_bytes guarantees the
     * loop stops before reaching it. */
    while ((cache->count > cache->max_entries || cache->bytes > cache->max_bytes) &&
           cache->tail != e) {
        entry_drop(cache, cache->tail);
        ++cache->stats.evictions;
    }
    return true;
}

void cache_clear(Cache *cache)
{
    if (cache == NULL) {
        return;
    }
    while (cache->head != NULL) {
        entry_drop(cache, cache->head);
        ++cache->stats.invalidations;
    }
}

CacheStats cache_stats(const Cache *cache)
{
    CacheStats s;
    memset(&s, 0, sizeof(s));
    if (cache != NULL) {
        s         = cache->stats;
        s.entries = cache->count;
        s.bytes   = cache->bytes;
    }
    return s;
}
