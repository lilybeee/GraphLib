/*
 * migrate.c — rebuild a graph in a different representation.  See migrate.h.
 *
 * The copy goes through the public backend interface (iterate the source's edges, add them to
 * the destination) rather than any representation-specific shortcut, so it works between ANY
 * pair of backends, including a test stub.  The old graph is never modified: if the copy
 * fails part-way the half-built destination is destroyed and the caller keeps the original,
 * which is what makes both manual conversion and automatic migration safe.  Vertex IDs are
 * preserved exactly (dead IDs included) so nothing a caller holds becomes invalid, and a
 * weight-free graph stays weight-free where the target supports it.
 */

#include "migrate.h"

void *graph_clone(const GraphBackend *from, const void *src, const GraphBackend *to)
{
    if (from == NULL || src == NULL || to == NULL) {
        return NULL;
    }

    const size_t bound    = from->id_bound(src);
    const bool   directed = from->is_directed(src);

    /* create() makes IDs 0 .. bound-1 alive; remove the ones that are dead
     * in the source so the ID space matches exactly. */
    /* Keep a weight-free graph weight-free where the target can; otherwise its 1.0 weights are
     * simply stored. */
    const bool plain = from->is_unweighted != NULL && from->is_unweighted(src) &&
                       to->create_unweighted != NULL;
    void *dst = plain ? to->create_unweighted(bound, directed) : to->create(bound, directed);
    if (dst == NULL) {
        return NULL;
    }
    for (size_t id = 0; id < bound; ++id) {
        if (!from->is_alive(src, (int)id) && !to->remove_vertex(dst, (int)id)) {
            goto fail;
        }
    }

    bool (*add)(void *, int, int, double) =
        to->add_new_edge != NULL ? to->add_new_edge : to->add_edge;

    for (size_t u = 0; u < bound; ++u) {
        if (!from->is_alive(src, (int)u)) {
            continue;
        }
        GraphIter it;
        int       v;
        double    w;
        from->iter_begin(src, (int)u, &it);
        while (from->iter_next(&it, &v, &w)) {
            /* An undirected edge shows up in both endpoints' lists (a
             * self-loop once); add_edge mirrors it, so take it only once. */
            if (!directed && v < (int)u) {
                continue;
            }
            if (!add(dst, (int)u, v, w)) {
                goto fail;
            }
        }
    }

    /* Defensive: a backend that silently dropped something is worse than a
     * failed migration. */
    if (to->vertex_count(dst) != from->vertex_count(src) ||
        to->edge_count(dst)   != from->edge_count(src)) {
        goto fail;
    }
    return dst;

fail:
    to->destroy(dst);
    return NULL;
}
