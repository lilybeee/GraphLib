/*
 * graph_social.h — queries a social application asks about ONE user or ONE pair
 * of users: who is within two hops, who might I know, how are we connected.
 *
 * Unlike the whole-graph algorithms in graph_algorithms.h, the cost of these
 * depends on the neighbourhood they touch, not on the number of users: asking
 * for the friends of friends of one person in a graph of a million users takes
 * about as long as in a graph of a thousand (as long as the people involved
 * have similar numbers of friends).  They keep their working state in small hash
 * tables sized to what they visit, so there is no per-query O(V) allocation.
 *
 * All of them work on either representation, only read the graph, and return a
 * GraphStatus.  "Friends" of a user are the out-neighbours (the people they
 * follow, in a directed graph).  Users with no connections are fine: they simply
 * have no neighbourhood.
 */

#ifndef GRAPH_SOCIAL_H
#define GRAPH_SOCIAL_H

#include <stdbool.h>
#include <stddef.h>

#include "graph.h"
#include "graph_results.h"
#include "graph_status.h"

/**
 * graph_social_within_hops - everyone reachable from @source in at most
 * @max_hops steps, nearest first (ties in discovery order), excluding @source.
 *
 * @users / @hops (optional) receive up to @max_results entries; *count is how
 * many were written.  If more people qualify than fit, the nearest ones are
 * kept and *truncated (optional) is set.  Cost O(edges out of the visited people).
 * Errors: INVALID_VERTEX for a dead @source.
 */
GraphStatus graph_social_within_hops(Graph *g, int source, size_t max_hops, size_t max_results,
                                     int *users, int *hops, size_t *count, bool *truncated);

/** A candidate friend and how many friends they share with the user. */
typedef struct {
    int    user;
    size_t mutual;
} GraphSuggestion;

/**
 * graph_social_suggest_friends - "people you may know": the (up to) @k people who
 * are not @user and not already their friends, ranked by how many friends they
 * share with @user (most first, ties by lower ID).  *count receives how many
 * were written to @out.  Cost O(sum of the friends' degrees).
 */
GraphStatus graph_social_suggest_friends(Graph *g, int user, size_t k, GraphSuggestion *out,
                                         size_t *count);

/**
 * graph_social_path - a shortest chain of friends from @source to @target by
 * bidirectional breadth-first search, giving up beyond @max_hops
 * (SIZE_MAX = no limit).  @out->vertices lists the users from source to target
 * and @out->cost is the number of hops; free it with graph_path_free().
 * Undirected graphs only (UNSUPPORTED otherwise).  GRAPH_ERR_NOT_FOUND if they
 * are not connected within @max_hops.  Searching from both ends visits roughly
 * the square root of what a one-sided search would.
 */
GraphStatus graph_social_path(Graph *g, int source, int target, size_t max_hops, GraphPath *out);

#endif /* GRAPH_SOCIAL_H */
