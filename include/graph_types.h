/*
 * graph_types.h — small types shared by every representation.
 */

#ifndef GRAPH_TYPES_H
#define GRAPH_TYPES_H

/* A weighted edge; used by MST results and anywhere an edge is passed by value. */
typedef struct {
    int    u;
    int    v;
    double weight;
} GraphEdge;

#endif /* GRAPH_TYPES_H */
