/*
 * mapping.h — bijective value ↔ integer-ID mapping
 *
 * Maps arbitrary user data (strings, ints, pointers, …) to dense integer
 * IDs 0..n-1, and back.  The mapping is bijective:
 *
 *   value  →  ID   via  mapping_get_or_create / mapping_find
 *   ID     →  value via  mapping_data
 *
 * The caller supplies a KeyOps struct that teaches the map how to hash,
 * compare, copy, and free keys.  Two ready-made KeyOps are provided for
 * the common cases of NUL-terminated strings and plain int values.
 *
 * Internally the map uses a hash table (djb2, open addressing, 50% max
 * load factor) for O(1) average value→ID lookup, and a dynamic array
 * for O(1) ID→value lookup.
 *
 * Thread safety: none.  Callers must serialise if shared.
 */

#ifndef MAPPING_H
#define MAPPING_H

#include <stdbool.h>
#include <stddef.h>

/* ------------------------------------------------------------------ */
/*  KeyOps — caller-supplied type descriptor                           */
/* ------------------------------------------------------------------ */

/**
 * KeyOps - describes how to work with a particular key type.
 *
 * @hash:     Hash the key.  Must be consistent with @equals.
 * @equals:   Return true iff two keys are equal.
 * @copy:     Allocate and return a deep copy of the key, or NULL to store
 *            the caller's pointer as-is (no ownership transfer).
 * @free_key: Release a copied key.  NULL = mapping never frees keys.
 *
 * If @copy is non-NULL the mapping takes ownership of copied keys and
 * frees them with @free_key on mapping_destroy().
 */
typedef struct KeyOps {
    size_t (*hash)    (const void *key);
    bool   (*equals)  (const void *a, const void *b);
    void  *(*copy)    (const void *key);   /* NULL = store pointer as-is */
    void   (*free_key)(void *key);         /* NULL = never free          */
} KeyOps;

/* ---- Ready-made KeyOps ---- */

/**
 * KEYOPS_STRING - for NUL-terminated char* keys.
 * Keys are deep-copied (strdup) and freed on destroy.
 * Hash: djb2.
 */
extern const KeyOps KEYOPS_STRING;

/**
 * KEYOPS_INT - for int* keys (passed as (const void *)&my_int).
 * Keys are copied (sizeof(int) malloc) and freed on destroy.
 */
extern const KeyOps KEYOPS_INT;

/* ------------------------------------------------------------------ */
/*  Mapping                                                             */
/* ------------------------------------------------------------------ */

typedef struct Mapping Mapping;

/**
 * mapping_create - allocate a new empty mapping.
 * @ops: key-type descriptor (copied by value into the mapping).
 * Returns NULL on allocation failure.
 */
Mapping    *mapping_create(KeyOps ops);

/**
 * mapping_destroy - free all memory owned by the mapping, including
 * any copied keys.  Safe to call with NULL.
 */
void        mapping_destroy(Mapping *m);

/**
 * mapping_get_or_create - look up key; create a new ID if absent.
 *
 * If the key is already in the map, returns its existing ID.
 * Otherwise copies the key (via KeyOps.copy), assigns the next
 * sequential ID (equal to the count before insertion), and returns it.
 *
 * Returns -1 on allocation failure.
 * Complexity: O(1) amortised.
 *
 * Invariant after n distinct successful calls:
 *   IDs are exactly 0 .. n-1  (dense).
 */
int         mapping_get_or_create(Mapping *m, const void *key);

/**
 * mapping_find - look up key; return -1 if absent.
 * Never creates anything.
 * Complexity: O(1) amortised.
 */
int         mapping_find(const Mapping *m, const void *key);

/**
 * mapping_assign - like mapping_get_or_create, but the caller chooses the ID.
 * @id must be unused and at most mapping_id_bound().  Returns @id, or -1 if the
 * key already exists, the ID is taken, or allocation fails.  This is how
 * NamedGraph keeps the mapping in step with the IDs a Graph hands out
 * (a Graph reuses the IDs of removed vertices).
 */
int         mapping_assign(Mapping *m, const void *key, int id);

/**
 * mapping_remove - forget the key stored under @id and free it.  The ID becomes
 * unused (mapping_data(id) returns NULL).  False if @id is not in use.
 * Complexity: O(1) amortised.
 */
bool        mapping_remove(Mapping *m, int id);

/**
 * mapping_compact - renumber after Graph compaction.  @remap has @n >=
 * mapping_id_bound() entries (remap[old] = new ID, or -1 for unused IDs), as
 * produced by graph_compact.  Surviving IDs must only move down.  After this,
 * the IDs in use are exactly 0 .. mapping_count()-1.  False on a bad @remap
 * (mapping unchanged) or allocation failure.
 */
bool        mapping_compact(Mapping *m, const int *remap, size_t n);

/** One past the largest ID ever handed out (removed IDs inside it are unused). */
size_t      mapping_id_bound(const Mapping *m);

/**
 * mapping_data - return the key associated with @id, or NULL if out of range.
 * The returned pointer is owned by the mapping; do not free it.
 * Complexity: O(1).
 */
const void *mapping_data(const Mapping *m, int id);

/**
 * mapping_count - return the number of distinct keys in the map (removed ones excluded).
 * Complexity: O(1).
 */
size_t      mapping_count(const Mapping *m);

/* ---- Convenience wrapper: graph + mapping together ---- */
typedef struct Graph Graph;

typedef struct NamedGraph {
    Graph   *g;
    Mapping *map;
} NamedGraph;

NamedGraph *named_graph_create(bool directed, KeyOps ops);
void        named_graph_destroy(NamedGraph *ng);
int         named_add_vertex(NamedGraph *ng, const void *key);
bool        named_add_edge(NamedGraph *ng, const void *u, const void *v, double w);
const void *named_vertex_data(const NamedGraph *ng, int id);

/**
 * named_remove_vertex - remove the vertex called @key and all its edges.  The
 * vertex's ID is recycled by the next named_add_vertex.  False (with
 * graph_last_error() = GRAPH_ERR_NOT_FOUND) if no such key exists.
 */
bool        named_remove_vertex(NamedGraph *ng, const void *key);

/** Remove the edge between the vertices called @u and @v. */
bool        named_remove_edge(NamedGraph *ng, const void *u, const void *v);

/**
 * named_compact - compact the underlying Graph and renumber the mapping to
 * match.  Names stay attached to the right vertices; IDs you saved from
 * named_vertex_id() before the call are invalid afterwards.
 */
bool        named_compact(NamedGraph *ng);
int         named_vertex_id(const NamedGraph *ng, const void *key);

#endif /* MAPPING_H */
