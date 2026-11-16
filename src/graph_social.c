/*
 * graph_social.c — see graph_social.h.
 *
 * C note: the "visited" sets are small open-addressing hash tables keyed by user
 * ID.  A flat array indexed by ID would be simpler, but it would cost O(V) to
 * allocate and clear on every query, which is exactly what these queries avoid.
 */

#include "graph_social.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/*  A tiny int -> (int, int) hash map                                   */
/* ------------------------------------------------------------------ */

typedef struct {
    int    *key;     /* -1 = empty slot */
    int    *a;       /* first payload: parent, or a count */
    int    *b;       /* second payload: hop distance */
    size_t  cap;     /* power of two */
    size_t  n;
} Map;

static bool map_init(Map *m, size_t want)
{
    size_t cap = 16;
    while (cap < want * 2) {
        cap *= 2;
    }
    m->key = malloc(cap * sizeof *m->key);
    m->a   = malloc(cap * sizeof *m->a);
    m->b   = malloc(cap * sizeof *m->b);
    m->cap = cap;
    m->n   = 0;
    if (m->key == NULL || m->a == NULL || m->b == NULL) {
        free(m->key); free(m->a); free(m->b);
        return false;
    }
    memset(m->key, 0xff, cap * sizeof *m->key);
    return true;
}

static void map_free(Map *m)
{
    free(m->key); free(m->a); free(m->b);
}

static size_t map_slot(const Map *m, int key)
{
    size_t i = ((uint32_t)key * 2654435761u) & (m->cap - 1);
    while (m->key[i] != -1 && m->key[i] != key) {
        i = (i + 1) & (m->cap - 1);
    }
    return i;
}

static bool map_find(const Map *m, int key, size_t *slot)
{
    const size_t i = map_slot(m, key);
    *slot = i;
    return m->key[i] == key;
}

static bool map_put(Map *m, int key, int a, int b)
{
    if ((m->n + 1) * 2 > m->cap) {                  /* keep the load factor under one half */
        Map bigger;
        if (!map_init(&bigger, m->cap)) {
            return false;
        }
        for (size_t i = 0; i < m->cap; i++) {
            if (m->key[i] != -1) {
                const size_t j = map_slot(&bigger, m->key[i]);
                bigger.key[j] = m->key[i];
                bigger.a[j]   = m->a[i];
                bigger.b[j]   = m->b[i];
                bigger.n++;
            }
        }
        map_free(m);
        *m = bigger;
    }
    const size_t i = map_slot(m, key);
    if (m->key[i] == -1) {
        m->key[i] = key;
        m->n++;
    }
    m->a[i] = a;
    m->b[i] = b;
    return true;
}

/* ------------------------------------------------------------------ */
/*  Reading one user's friends                                          */
/* ------------------------------------------------------------------ */

typedef struct {
    int *buf;
    size_t cap;
} Nbrs;

/* The friends of v, into n->buf.  Returns how many, or SIZE_MAX if out of memory. */
static size_t friends_of(Graph *g, int v, Nbrs *n)
{
    const size_t d = graph_degree(g, v);
    if (d > n->cap) {
        size_t cap = n->cap != 0 ? n->cap : 64;
        while (cap < d) {
            cap *= 2;
        }
        int *bigger = realloc(n->buf, cap * sizeof *bigger);
        if (bigger == NULL) {
            return SIZE_MAX;
        }
        n->buf = bigger;
        n->cap = cap;
    }
    return d == 0 ? 0 : graph_neighbors(g, v, n->buf, NULL, d);
}

/* ------------------------------------------------------------------ */
/*  Everyone within k hops                                              */
/* ------------------------------------------------------------------ */

GraphStatus graph_social_within_hops(Graph *g, int source, size_t max_hops, size_t max_results,
                                     int *users, int *hops, size_t *count, bool *truncated)
{
    if (g == NULL || count == NULL) {
        return graph_set_error(GRAPH_ERR_INVALID_ARGUMENT);
    }
    *count = 0;
    if (truncated != NULL) {
        *truncated = false;
    }
    if (!graph_is_alive(g, source)) {
        return graph_set_error(GRAPH_ERR_INVALID_VERTEX);
    }

    Map seen;
    Nbrs nb = {0};
    size_t queue_cap = 64, head = 0, tail = 0;
    int *queue = malloc(queue_cap * sizeof *queue);
    GraphStatus st = GRAPH_OK;
    if (queue == NULL || !map_init(&seen, 64)) {
        free(queue);
        return graph_set_error(GRAPH_ERR_OUT_OF_MEMORY);
    }
    map_put(&seen, source, -1, 0);
    queue[tail++] = source;

    while (head < tail && st == GRAPH_OK) {
        const int u = queue[head++];
        size_t slot;
        map_find(&seen, u, &slot);
        const int du = seen.b[slot];
        if ((size_t)du >= max_hops) {
            continue;                                  /* do not look past the requested depth */
        }
        const size_t d = friends_of(g, u, &nb);
        if (d == SIZE_MAX) {
            st = GRAPH_ERR_OUT_OF_MEMORY;
            break;
        }
        for (size_t i = 0; i < d; i++) {
            const int w = nb.buf[i];
            if (map_find(&seen, w, &slot)) {
                continue;
            }
            if (*count == max_results) {               /* full: keep the nearest, say so */
                if (truncated != NULL) {
                    *truncated = true;
                }
                head = tail;                           /* nothing nearer can still come: BFS order */
                break;
            }
            if (!map_put(&seen, w, u, du + 1)) {
                st = GRAPH_ERR_OUT_OF_MEMORY;
                break;
            }
            if (users != NULL) users[*count] = w;
            if (hops != NULL)  hops[*count]  = du + 1;
            (*count)++;
            if (tail == queue_cap) {
                int *bigger = realloc(queue, queue_cap * 2 * sizeof *bigger);
                if (bigger == NULL) {
                    st = GRAPH_ERR_OUT_OF_MEMORY;
                    break;
                }
                queue = bigger;
                queue_cap *= 2;
            }
            queue[tail++] = w;
        }
    }
    free(queue);
    free(nb.buf);
    map_free(&seen);
    return st == GRAPH_OK ? GRAPH_OK : graph_set_error(st);
}

/* ------------------------------------------------------------------ */
/*  People you may know                                                 */
/* ------------------------------------------------------------------ */

static int suggestion_order(const void *x, const void *y)
{
    const GraphSuggestion *a = x, *b = y;
    if (a->mutual != b->mutual) return a->mutual > b->mutual ? -1 : 1;
    return (a->user > b->user) - (a->user < b->user);
}

GraphStatus graph_social_suggest_friends(Graph *g, int user, size_t k, GraphSuggestion *out,
                                         size_t *count)
{
    if (g == NULL || count == NULL || (k > 0 && out == NULL)) {
        return graph_set_error(GRAPH_ERR_INVALID_ARGUMENT);
    }
    *count = 0;
    if (!graph_is_alive(g, user)) {
        return graph_set_error(GRAPH_ERR_INVALID_VERTEX);
    }
    const size_t deg = graph_degree(g, user);
    if (deg == 0 || k == 0) {
        return GRAPH_OK;                               /* nobody to be introduced through */
    }

    Nbrs mine = {0}, theirs = {0};
    const size_t d0 = friends_of(g, user, &mine);
    Map tally, friends;
    GraphStatus st = GRAPH_OK;
    if (d0 == SIZE_MAX || !map_init(&friends, deg) || !map_init(&tally, deg * 4)) {
        free(mine.buf);
        return graph_set_error(GRAPH_ERR_OUT_OF_MEMORY);
    }
    for (size_t i = 0; i < d0; i++) {
        map_put(&friends, mine.buf[i], 0, 0);
    }
    map_put(&friends, user, 0, 0);                     /* never suggest the user to themself */

    for (size_t i = 0; i < d0 && st == GRAPH_OK; i++) {
        const size_t d = friends_of(g, mine.buf[i], &theirs);
        if (d == SIZE_MAX) {
            st = GRAPH_ERR_OUT_OF_MEMORY;
            break;
        }
        for (size_t j = 0; j < d; j++) {
            const int w = theirs.buf[j];
            size_t slot;
            if (map_find(&friends, w, &slot)) {
                continue;                              /* already a friend (or the user) */
            }
            if (map_find(&tally, w, &slot)) {
                tally.a[slot]++;
            } else if (!map_put(&tally, w, 1, 0)) {
                st = GRAPH_ERR_OUT_OF_MEMORY;
                break;
            }
        }
    }
    if (st == GRAPH_OK) {
        GraphSuggestion *all = malloc((tally.n > 0 ? tally.n : 1) * sizeof *all);
        if (all == NULL) {
            st = GRAPH_ERR_OUT_OF_MEMORY;
        } else {
            size_t m = 0;
            for (size_t i = 0; i < tally.cap; i++) {
                if (tally.key[i] != -1) {
                    all[m].user   = tally.key[i];
                    all[m].mutual = (size_t)tally.a[i];
                    m++;
                }
            }
            qsort(all, m, sizeof *all, suggestion_order);
            *count = m < k ? m : k;
            memcpy(out, all, *count * sizeof *out);
            free(all);
        }
    }
    free(mine.buf);
    free(theirs.buf);
    map_free(&tally);
    map_free(&friends);
    return st == GRAPH_OK ? GRAPH_OK : graph_set_error(st);
}

/* ------------------------------------------------------------------ */
/*  Shortest chain of friends, searching from both ends                 */
/* ------------------------------------------------------------------ */

/* One side of the search: who it has reached (parent, hops) and the current level. */
typedef struct {
    Map    seen;
    int   *level;
    size_t level_n;
    int    depth;
} Side;

static bool side_init(Side *s, int root)
{
    s->level   = malloc(sizeof *s->level);
    s->level_n = 1;
    s->depth   = 0;
    if (s->level == NULL || !map_init(&s->seen, 64)) {
        free(s->level);
        return false;
    }
    s->level[0] = root;
    return map_put(&s->seen, root, -1, 0);
}

static void side_free(Side *s)
{
    free(s->level);
    map_free(&s->seen);
}

GraphStatus graph_social_path(Graph *g, int source, int target, size_t max_hops, GraphPath *out)
{
    if (g == NULL || out == NULL) {
        return graph_set_error(GRAPH_ERR_INVALID_ARGUMENT);
    }
    *out = (GraphPath){0};
    if (!graph_is_alive(g, source) || !graph_is_alive(g, target)) {
        return graph_set_error(GRAPH_ERR_INVALID_VERTEX);
    }
    if (graph_is_directed(g)) {
        return graph_set_error(GRAPH_ERR_UNSUPPORTED);   /* needs in-neighbours, which a directed graph lacks */
    }
    if (source == target) {
        out->vertices = malloc(sizeof *out->vertices);
        if (out->vertices == NULL) {
            return graph_set_error(GRAPH_ERR_OUT_OF_MEMORY);
        }
        out->vertices[0] = source;
        out->length      = 1;
        out->cost        = 0;
        return GRAPH_OK;
    }

    Side fwd = {0}, bwd = {0};
    Nbrs nb  = {0};
    GraphStatus st = GRAPH_ERR_NOT_FOUND;
    int meet_a = -1, meet_b = -1;                       /* an edge a-b joining the two searches */
    if (!side_init(&fwd, source) || !side_init(&bwd, target)) {
        side_free(&fwd);
        return graph_set_error(GRAPH_ERR_OUT_OF_MEMORY);
    }

    while (fwd.level_n > 0 && bwd.level_n > 0 && (size_t)(fwd.depth + bwd.depth) < max_hops) {
        /* Expand the side with the smaller frontier by one whole level. */
        Side *mine  = fwd.level_n <= bwd.level_n ? &fwd : &bwd;
        Side *other = mine == &fwd ? &bwd : &fwd;
        int  *next   = NULL;
        size_t next_n = 0, next_cap = 0;
        int best = -1;                                  /* shortest total over this level */
        for (size_t i = 0; i < mine->level_n && st != GRAPH_ERR_OUT_OF_MEMORY; i++) {
            const int u = mine->level[i];
            const size_t d = friends_of(g, u, &nb);
            if (d == SIZE_MAX) {
                st = GRAPH_ERR_OUT_OF_MEMORY;
                break;
            }
            for (size_t j = 0; j < d; j++) {
                const int w = nb.buf[j];
                size_t slot;
                if (map_find(&other->seen, w, &slot)) {          /* the two searches touch */
                    const int total = mine->depth + 1 + other->seen.b[slot];
                    if (best < 0 || total < best) {
                        best = total;
                        meet_a = mine == &fwd ? u : w;           /* a is on the source side */
                        meet_b = mine == &fwd ? w : u;
                    }
                }
                if (map_find(&mine->seen, w, &slot)) {
                    continue;
                }
                if (!map_put(&mine->seen, w, u, mine->depth + 1)) {
                    st = GRAPH_ERR_OUT_OF_MEMORY;
                    break;
                }
                if (next_n == next_cap) {
                    next_cap = next_cap != 0 ? next_cap * 2 : 64;
                    int *bigger = realloc(next, next_cap * sizeof *bigger);
                    if (bigger == NULL) {
                        st = GRAPH_ERR_OUT_OF_MEMORY;
                        break;
                    }
                    next = bigger;
                }
                next[next_n++] = w;
            }
        }
        free(mine->level);
        mine->level   = next;
        mine->level_n = next_n;
        mine->depth++;
        if (st == GRAPH_ERR_OUT_OF_MEMORY) {
            break;
        }
        if (best >= 0) {
            if ((size_t)best > max_hops) {
                break;                                    /* connected, but farther than allowed */
            }
            st = GRAPH_OK;
            break;
        }
    }

    if (st == GRAPH_OK) {
        /* source ... meet_a, meet_b ... target */
        size_t na = 1, nbk = 1;
        size_t slot;
        for (int v = meet_a; map_find(&fwd.seen, v, &slot) && fwd.seen.a[slot] != -1; v = fwd.seen.a[slot]) na++;
        for (int v = meet_b; map_find(&bwd.seen, v, &slot) && bwd.seen.a[slot] != -1; v = bwd.seen.a[slot]) nbk++;
        int *path = malloc((na + nbk) * sizeof *path);
        if (path == NULL) {
            st = GRAPH_ERR_OUT_OF_MEMORY;
        } else {
            size_t pos = na;
            int v = meet_a;
            while (pos > 0) {
                path[--pos] = v;
                map_find(&fwd.seen, v, &slot);
                v = fwd.seen.a[slot];
            }
            pos = na;
            v   = meet_b;
            while (pos < na + nbk) {
                path[pos++] = v;
                map_find(&bwd.seen, v, &slot);
                v = bwd.seen.a[slot];
            }
            out->vertices = path;
            out->length   = na + nbk;
            out->cost     = (double)(na + nbk - 1);
        }
    }
    free(nb.buf);
    side_free(&fwd);
    side_free(&bwd);
    return st == GRAPH_OK ? GRAPH_OK : graph_set_error(st);
}
