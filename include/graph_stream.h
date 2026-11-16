/*
 * graph_stream.h — analyses that read an edge-list FILE and never build the graph.
 *
 * For a network too big to hold as a graph, these make one or two passes over the
 * file and keep only a few integers per user (memory O(V), nothing per edge), so
 * a billion friendships can be summarised on a machine that could not store them.
 * They are the "beyond RAM" tools: what they can answer is limited to questions
 * that can be settled while streaming (how many users and links, how many friends
 * each user has, which users are connected to which).  Anything that needs the
 * whole graph at once (shortest paths, centrality, communities) still needs it in
 * memory, or a sample of it.
 *
 * File format: the one graph_write_edge_list writes and graph_read_edge_list reads —
 * one "u v [weight]" per line, '#' comment lines, with optional "# directed",
 * "# undirected" and "# vertices N" headers.  Each function reads from the current
 * position to the end, so rewind() between calls.
 */

#ifndef GRAPH_STREAM_H
#define GRAPH_STREAM_H

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>

#include "graph_status.h"

typedef struct {
    size_t vertices;   /* largest user ID + 1, or the "# vertices" header if larger */
    size_t edges;      /* edge lines read (an undirected edge is one line)           */
    bool   directed;   /* from the header; false if the file does not say            */
} GraphStreamInfo;

/** Pass 1: learn the size.  O(1) memory.  GRAPH_ERR_IO on a malformed line. */
GraphStatus graph_stream_scan(FILE *in, GraphStreamInfo *info);

/**
 * Degrees of every user.  @n = info.vertices from graph_stream_scan.  out_degree[v]
 * (n entries) counts the edges leaving v — for an undirected file, every edge touching
 * v.  @in_degree (optional, n entries) counts arrows arriving (directed files only).
 * Both arrays are zeroed first.  GRAPH_ERR_INVALID_VERTEX if a line names a user >= n.
 */
GraphStatus graph_stream_degrees(FILE *in, size_t n, bool directed, size_t *out_degree,
                                 size_t *in_degree);

/**
 * Connected components (directions ignored) by union-find over the streamed edges.
 * @component (n entries) gets labels 0.. in order of each component's first member;
 * *count (optional) the number of components.  Users that appear in no edge are
 * components of their own.  Memory: two int arrays of n entries plus @component.
 */
GraphStatus graph_stream_components(FILE *in, size_t n, int *component, size_t *count);

/** Called once per edge line, in file order.  Return anything but GRAPH_OK to stop the walk. */
typedef GraphStatus (*GraphStreamEdgeFn)(void *ctx, int u, int v);

/**
 * Walks the file and hands every edge to @fn, one at a time, holding nothing in memory.  This
 * is the building block the sharding engine (graph_shard.h) uses to bucket a file's edges
 * without loading it.  GRAPH_ERR_IO on a malformed line; @fn's own error is returned as is.
 */
GraphStatus graph_stream_edges(FILE *in, GraphStreamEdgeFn fn, void *ctx);

#endif /* GRAPH_STREAM_H */
