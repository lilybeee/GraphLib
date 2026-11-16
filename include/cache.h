/*
 * cache.h — versioned LRU cache for algorithm results.
 *
 * Results are keyed by (algorithm, up to two integer arguments) and stored as
 * opaque byte blobs.  The cache belongs to one graph and is told that graph's
 * version on every call (list_graph_version / GraphBackend.version): the first
 * call that sees a different version drops every entry, so a stale result can
 * never be returned and no per-mutation bookkeeping is needed.
 *
 * Results do not depend on the representation, so a migration (which keeps
 * the graph's contents) does not need to clear the cache as long as the caller
 * carries the version counter over.
 *
 * Eviction is least-recently-used, bounded by entry count and by total payload
 * bytes.  Lookup, store and eviction are O(1) expected.
 */

#ifndef CACHE_H
#define CACHE_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

typedef enum {
    CACHE_BFS,
    CACHE_DFS,
    CACHE_COMPONENTS,
    CACHE_CYCLE,
    CACHE_TOPO,
    CACHE_DIJKSTRA,
    CACHE_MST,
    CACHE_USER = 64     /* first value free for the caller's own keys */
} CacheAlgo;

typedef struct {
    uint32_t algo;      /* CacheAlgo or a CACHE_USER-based value          */
    int      a;         /* e.g. the source vertex; 0 if unused            */
    int      b;         /* e.g. the target vertex; 0 if unused            */
} CacheKey;

typedef struct {
    uint64_t hits;
    uint64_t misses;
    uint64_t stores;
    uint64_t evictions;     /* removed to stay within the limits          */
    uint64_t invalidations; /* removed because the graph changed          */
    size_t   entries;
    size_t   bytes;         /* payload bytes currently held               */
} CacheStats;

typedef struct Cache Cache;

/**
 * cache_create - make an empty cache.
 * @max_entries: most entries held at once (0 disables caching).
 * @max_bytes:   most payload bytes held at once; a single result larger than
 *               this is never stored.
 * Returns NULL on allocation failure.
 */
Cache *cache_create(size_t max_entries, size_t max_bytes);

/** Frees the cache and every stored result.  Safe with NULL. */
void cache_destroy(Cache *cache);

/**
 * cache_lookup - fetch a stored result.
 * Returns a pointer to the stored bytes (and sets *size if non-NULL), or NULL
 * on a miss.  The pointer stays valid until the next call on this cache.
 * A hit makes the entry the most recently used.
 */
const void *cache_lookup(Cache *cache, CacheKey key, uint64_t graph_version,
                         size_t *size);

/**
 * cache_store - copy a result into the cache, evicting the least recently
 * used entries if needed.  Replaces an existing entry with the same key.
 * Returns false if the result was not stored (too large, cache disabled, or
 * allocation failure); the cache is still consistent in that case.
 */
bool cache_store(Cache *cache, CacheKey key, uint64_t graph_version,
                 const void *data, size_t size);

/** Drops every entry (counted as invalidations). */
void cache_clear(Cache *cache);

CacheStats cache_stats(const Cache *cache);

#endif /* CACHE_H */
