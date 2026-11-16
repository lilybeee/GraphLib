/*
 * mapping.c — bijective value ↔ integer-ID mapping
 *
 * Implementation details
 * ----------------------
 *  ID → value : a dynamic array (data[]).  data[id] = copied key.
 *  value → ID : an open-addressing hash table (slots[]).
 *               Each slot stores either EMPTY, DELETED, or a (key, id) pair.
 *               Load factor is kept below 50 % before any lookup so that
 *               the worst-case probe sequence stays short.
 *
 * Hash function: djb2 (strings) / FNV-mix (raw bytes for ints).
 *                User supplies their own via KeyOps.hash.
 *
 * Growth: the table doubles when load > 50 %, rehashing all live entries.
 */

#include "mapping.h"

#include <stdlib.h>
#include <string.h>
#include <stdint.h>

/* ------------------------------------------------------------------ */
/*  Ready-made KeyOps implementations                                   */
/* ------------------------------------------------------------------ */

/* ---- string ---- */

static size_t str_hash(const void *key)
{
    /* djb2 */
    const unsigned char *p = key;
    size_t h = 5381;
    while (*p) {
        h = ((h << 5) + h) ^ (size_t)*p;
        ++p;
    }
    return h;
}

static bool str_equals(const void *a, const void *b)
{
    return strcmp((const char *)a, (const char *)b) == 0;
}

static void *str_copy(const void *key)
{
    const char *s = key;
    size_t len = strlen(s) + 1;
    char *copy = malloc(len);
    if (copy) memcpy(copy, s, len);
    return copy;
}

const KeyOps KEYOPS_STRING = {
    .hash     = str_hash,
    .equals   = str_equals,
    .copy     = str_copy,
    .free_key = free,
};

/* ---- int ---- */

static size_t int_hash(const void *key)
{
    /* FNV-1a inspired mix */
    unsigned int v = *(const int *)key;
    size_t h = (size_t)v;
    h ^= h >> 16;
    h *= 0x45d9f3b;
    h ^= h >> 16;
    return h;
}

static bool int_equals(const void *a, const void *b)
{
    return *(const int *)a == *(const int *)b;
}

static void *int_copy(const void *key)
{
    int *copy = malloc(sizeof(int));
    if (copy) *copy = *(const int *)key;
    return copy;
}

const KeyOps KEYOPS_INT = {
    .hash     = int_hash,
    .equals   = int_equals,
    .copy     = int_copy,
    .free_key = free,
};

/* ------------------------------------------------------------------ */
/*  Hash-table slot                                                     */
/* ------------------------------------------------------------------ */

typedef enum {
    SLOT_EMPTY   = 0,
    SLOT_LIVE    = 1,
    SLOT_DELETED = 2,
} SlotState;

typedef struct {
    SlotState  state;
    const void *key;   /* owned by the data[] array, NOT by the slot */
    int        id;
} Slot;

/* ------------------------------------------------------------------ */
/*  Mapping struct                                                      */
/* ------------------------------------------------------------------ */

struct Mapping {
    KeyOps  ops;

    /* ID → value */
    void  **data;      /* data[id] = copied key (owned here); NULL = removed ID */
    size_t  data_cap;
    size_t  bound;     /* IDs 0 .. bound-1 have been used (some may be removed) */
    size_t  count;     /* number of live entries */

    /* value → ID (open addressing) */
    Slot   *slots;
    size_t  slot_cap;  /* always a power of two */
    size_t  live;      /* live + deleted ≤ slot_cap */
};

/* ------------------------------------------------------------------ */
/*  Internal helpers                                                    */
/* ------------------------------------------------------------------ */


/* Find the slot for key in a table of given capacity.
 * Returns the index of: a LIVE slot with matching key (hit),
 *                       the first EMPTY slot encountered (miss),
 *                       or the first DELETED slot if no LIVE hit found (tombstone). */
static size_t probe(const Mapping *m, const void *key)
{
    size_t mask = m->slot_cap - 1;
    size_t idx  = m->ops.hash(key) & mask;
    size_t first_deleted = m->slot_cap; /* sentinel */

    for (size_t i = 0; i < m->slot_cap; ++i) {
        const Slot *s = &m->slots[idx];
        if (s->state == SLOT_EMPTY) {
            return (first_deleted < m->slot_cap) ? first_deleted : idx;
        }
        if (s->state == SLOT_LIVE && m->ops.equals(s->key, key)) {
            return idx;
        }
        if (s->state == SLOT_DELETED && first_deleted == m->slot_cap) {
            first_deleted = idx;
        }
        idx = (idx + 1) & mask;
    }
    /* Table is full of DELETED — return first_deleted */
    return first_deleted;
}

/* Rebuild the hash table with a new capacity (power of two). */
static bool rehash(Mapping *m, size_t new_cap)
{
    Slot *new_slots = calloc(new_cap, sizeof(*new_slots));
    if (new_slots == NULL) return false;

    size_t mask = new_cap - 1;
    for (size_t i = 0; i < m->slot_cap; ++i) {
        if (m->slots[i].state != SLOT_LIVE) continue;
        size_t idx = m->ops.hash(m->slots[i].key) & mask;
        while (new_slots[idx].state == SLOT_LIVE) {
            idx = (idx + 1) & mask;
        }
        new_slots[idx] = (Slot){ SLOT_LIVE, m->slots[i].key, m->slots[i].id };
    }

    free(m->slots);
    m->slots    = new_slots;
    m->slot_cap = new_cap;
    m->live     = m->count;     /* tombstones are gone after a rebuild */
    return true;
}

/* Grow the data[] array to at least @needed capacity. */
static bool grow_data(Mapping *m, size_t needed)
{
    if (needed <= m->data_cap) return true;
    size_t new_cap = m->data_cap == 0 ? 8 : m->data_cap * 2;
    while (new_cap < needed) new_cap *= 2;
    void **new_data = realloc(m->data, new_cap * sizeof(*new_data));
    if (new_data == NULL) return false;
    m->data     = new_data;
    m->data_cap = new_cap;
    return true;
}

/* ------------------------------------------------------------------ */
/*  Public API                                                          */
/* ------------------------------------------------------------------ */

Mapping *mapping_create(KeyOps ops)
{
    Mapping *m = malloc(sizeof(*m));
    if (m == NULL) return NULL;

    m->ops      = ops;
    m->data     = NULL;
    m->data_cap = 0;
    m->count    = 0;
    m->bound    = 0;
    m->slot_cap = 16; /* initial hash-table capacity */
    m->live     = 0;

    m->slots = calloc(m->slot_cap, sizeof(*m->slots));
    if (m->slots == NULL) {
        free(m);
        return NULL;
    }
    return m;
}

void mapping_destroy(Mapping *m)
{
    if (m == NULL) return;
    if (m->ops.free_key != NULL) {
        for (size_t i = 0; i < m->bound; ++i) {
            if (m->data[i] != NULL) {
                m->ops.free_key(m->data[i]);
            }
        }
    }
    free(m->data);
    free(m->slots);
    free(m);
}

/* Insert key under a specific, currently unused ID (id <= bound).  Returns id or -1. */
static int insert_with_id(Mapping *m, const void *key, int id)
{
    /* Grow hash table before lookup so probe() always finds a free slot. */
    if (m->live * 2 >= m->slot_cap) {
        if (!rehash(m, m->slot_cap * 2)) return -1;
    }

    size_t idx = probe(m, key);
    if (m->slots[idx].state == SLOT_LIVE) {
        return -1;      /* key already present */
    }

    /* New entry — copy the key */
    void *copied = m->ops.copy ? m->ops.copy(key) : (void *)key;
    if (m->ops.copy && copied == NULL) return -1;

    if (!grow_data(m, (size_t)id + 1)) {
        if (m->ops.free_key) m->ops.free_key(copied);
        return -1;
    }
    m->data[id] = copied;
    if ((size_t)id == m->bound) ++m->bound;
    ++m->count;

    /* Insert into hash table */
    bool was_deleted = (m->slots[idx].state == SLOT_DELETED);
    m->slots[idx] = (Slot){ SLOT_LIVE, copied, id };
    if (!was_deleted) ++m->live;

    return id;
}

int mapping_get_or_create(Mapping *m, const void *key)
{
    if (m == NULL || key == NULL) return -1;
    const int existing = mapping_find(m, key);
    if (existing >= 0) return existing;
    return insert_with_id(m, key, (int)m->bound);
}

int mapping_assign(Mapping *m, const void *key, int id)
{
    if (m == NULL || key == NULL || id < 0 || (size_t)id > m->bound) return -1;
    if ((size_t)id < m->bound && m->data[id] != NULL) return -1;   /* ID in use */
    return insert_with_id(m, key, id);
}

bool mapping_remove(Mapping *m, int id)
{
    if (m == NULL || id < 0 || (size_t)id >= m->bound || m->data[id] == NULL) return false;
    const size_t idx = probe(m, m->data[id]);
    if (m->slots[idx].state != SLOT_LIVE || m->slots[idx].id != id) return false;
    m->slots[idx].state = SLOT_DELETED;       /* tombstone: keeps probe chains intact */
    m->slots[idx].key   = NULL;
    if (m->ops.free_key) m->ops.free_key(m->data[id]);
    m->data[id] = NULL;
    --m->count;
    return true;
}

bool mapping_compact(Mapping *m, const int *remap, size_t n)
{
    if (m == NULL || (n != 0 && remap == NULL) || n < m->bound) return false;
    size_t next = 0;
    for (size_t old = 0; old < m->bound; ++old) {
        if (m->data[old] == NULL) continue;
        const int to = remap[old];
        if (to < 0 || (size_t)to > old) return false;     /* surviving IDs only move down */
        m->data[to] = m->data[old];
        if (to != (int)old) m->data[old] = NULL;
        ++next;
    }
    m->bound = next;
    /* The keys did not move in memory, but their IDs changed: rebuild the table. */
    Slot *fresh = calloc(m->slot_cap, sizeof(*fresh));
    if (fresh == NULL) return false;      /* data[] is already compact; see header note */
    const size_t mask = m->slot_cap - 1;
    for (size_t id = 0; id < m->bound; ++id) {
        size_t idx = m->ops.hash(m->data[id]) & mask;
        while (fresh[idx].state == SLOT_LIVE) idx = (idx + 1) & mask;
        fresh[idx] = (Slot){ SLOT_LIVE, m->data[id], (int)id };
    }
    free(m->slots);
    m->slots = fresh;
    m->live  = m->count;
    return true;
}

int mapping_find(const Mapping *m, const void *key)
{
    if (m == NULL || key == NULL || m->count == 0) return -1;
    size_t idx = probe(m, key);
    if (m->slots[idx].state == SLOT_LIVE) {
        return m->slots[idx].id;
    }
    return -1;
}

const void *mapping_data(const Mapping *m, int id)
{
    if (m == NULL || id < 0 || (size_t)id >= m->bound) return NULL;
    return m->data[id];
}

size_t mapping_count(const Mapping *m)
{
    return m ? m->count : 0;
}

size_t mapping_id_bound(const Mapping *m)
{
    return m ? m->bound : 0;
}

/* ------------------------------------------------------------------ */
/*  NamedGraph convenience wrapper                                      */
/* ------------------------------------------------------------------ */

#include "graph.h"
#include "graph_status.h"

NamedGraph *named_graph_create(bool directed, KeyOps ops)
{
    NamedGraph *ng = malloc(sizeof(*ng));
    if (!ng) return NULL;
    ng->g = graph_create(directed);
    ng->map = mapping_create(ops);
    if (!ng->g || !ng->map) {
        if (ng->g) graph_destroy(ng->g);
        if (ng->map) mapping_destroy(ng->map);
        free(ng);
        return NULL;
    }
    return ng;
}

void named_graph_destroy(NamedGraph *ng)
{
    if (!ng) return;
    graph_destroy(ng->g);
    mapping_destroy(ng->map);
    free(ng);
}

int named_add_vertex(NamedGraph *ng, const void *key)
{
    if (!ng || !key) return -1;
    int id = mapping_find(ng->map, key);
    if (id >= 0) return id;
    /* The graph picks the ID (it may reuse a removed one); the mapping follows. */
    id = graph_add_vertex(ng->g);
    if (id < 0) return -1;
    if (mapping_assign(ng->map, key, id) != id) {
        (void)graph_remove_vertex(ng->g, id);       /* roll back */
        return -1;
    }
    return id;
}

bool named_remove_vertex(NamedGraph *ng, const void *key)
{
    if (!ng || !key) {
        graph_set_error(GRAPH_ERR_INVALID_ARGUMENT);
        return false;
    }
    const int id = mapping_find(ng->map, key);
    if (id < 0) {
        graph_set_error(GRAPH_ERR_NOT_FOUND);
        return false;
    }
    if (!graph_remove_vertex(ng->g, id)) {
        return false;                               /* graph_last_error() says why */
    }
    return mapping_remove(ng->map, id);
}

bool named_remove_edge(NamedGraph *ng, const void *u, const void *v)
{
    if (!ng || !u || !v) {
        graph_set_error(GRAPH_ERR_INVALID_ARGUMENT);
        return false;
    }
    const int uid = mapping_find(ng->map, u);
    const int vid = mapping_find(ng->map, v);
    if (uid < 0 || vid < 0) {
        graph_set_error(GRAPH_ERR_NOT_FOUND);
        return false;
    }
    return graph_remove_edge(ng->g, uid, vid);
}

bool named_compact(NamedGraph *ng)
{
    if (!ng) {
        graph_set_error(GRAPH_ERR_INVALID_ARGUMENT);
        return false;
    }
    const size_t bound = graph_id_bound(ng->g);
    int *remap = malloc((bound > 0 ? bound : 1) * sizeof(*remap));
    if (remap == NULL) {
        graph_set_error(GRAPH_ERR_OUT_OF_MEMORY);
        return false;
    }
    bool ok = graph_compact(ng->g, remap) && mapping_compact(ng->map, remap, bound);
    free(remap);
    return ok;
}

bool named_add_edge(NamedGraph *ng, const void *u, const void *v, double w)
{
    if (!ng || !u || !v) return false;
    int uid = named_add_vertex(ng, u);
    int vid = named_add_vertex(ng, v);
    if (uid < 0 || vid < 0) return false;
    return graph_add_edge(ng->g, uid, vid, w);
}

const void *named_vertex_data(const NamedGraph *ng, int id)
{
    if (!ng) return NULL;
    return mapping_data(ng->map, id);
}

int named_vertex_id(const NamedGraph *ng, const void *key)
{
    if (!ng) return -1;
    return mapping_find(ng->map, key);
}

