/*
 * migrate.h — rebuild a graph in a different representation.
 */

#ifndef MIGRATE_H
#define MIGRATE_H

#include "backend.h"

/**
 * graph_clone - build a copy of the graph 'src' (stored in 'from') in the
 * representation 'to', and return the new handle (NULL on failure; 'src' is
 * never modified).
 *
 * The copy has the same directedness, the same ID space, the same set of
 * alive vertices and the same weighted edges, so every vertex ID the caller
 * holds stays valid.  Dead IDs are recreated as dead (and are available for
 * reuse in the copy).
 *
 * Complexity: O(V + E) plus whatever 'to' charges per call; creating the dead
 * IDs costs one remove_vertex each, so a graph with many dead IDs migrates
 * faster after list_graph_compact (which the caller can do beforehand when it
 * is prepared to apply the remap).
 */
void *graph_clone(const GraphBackend *from, const void *src, const GraphBackend *to);

#endif /* MIGRATE_H */
