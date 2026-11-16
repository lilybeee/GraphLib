/*
 * graph_internal.h — private definition of Graph, shared by graph.c and
 * graph_compare.c.  Not part of the public API.
 */

#ifndef GRAPH_INTERNAL_H
#define GRAPH_INTERNAL_H

#include "backend.h"
#include "cache.h"
#include "graph.h"
#include "policy.h"

struct Graph {
    const GraphBackend *be;        /* the representation in use            */
    void               *impl;      /* that representation's own handle     */
    Observer           *obs;
    Cache              *cache;     /* NULL when caching is disabled        */
    bool                adaptive;
    bool                frozen;    /* read-only, shareable between threads (graph_freeze) */
};

/** Record an operation and, in adaptive mode, let the observer migrate (no-op when frozen). */
void graph_after_op(Graph *g, OpKind op);

/** Wrap an existing backend handle in a new Graph that takes ownership of it.
 *  Settings (adaptive, cache sizes) are copied from @like, or the defaults if NULL.
 *  Returns NULL (and destroys @impl) on allocation failure. */
Graph *graph_adopt(const GraphBackend *be, void *impl, const Graph *like);

/** The backend that implements @rep. */
const GraphBackend *graph_backend_for(GraphRep rep);

#endif /* GRAPH_INTERNAL_H */
