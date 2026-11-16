/*
 * graph_status.h — error codes shared by every representation.
 *
 * Two ways a failure is reported:
 *
 *   1. Functions that return a GraphStatus (the algorithm library, e.g.
 *      matrix_bellman_ford) return the code directly: GRAPH_OK or a GRAPH_ERR_*.
 *   2. Functions that already return bool / size_t / double / int (the core
 *      graph operations and the baseline algorithms) keep their signature, and
 *      on failure also record a code you can read with graph_last_error().
 *
 * graph_last_error() is per thread (C11 _Thread_local), is only written when
 * something fails, and is never cleared by a success; call graph_clear_error()
 * before an operation if you want to be sure the value you read belongs to it.
 */

#ifndef GRAPH_STATUS_H
#define GRAPH_STATUS_H

typedef enum {
    GRAPH_OK = 0,
    GRAPH_ERR_INVALID_ARGUMENT,   /* NULL pointer, bad parameter                    */
    GRAPH_ERR_INVALID_VERTEX,     /* vertex ID out of range or removed              */
    GRAPH_ERR_DUPLICATE_EDGE,     /* the edge is already present                    */
    GRAPH_ERR_EDGE_NOT_FOUND,     /* the edge to remove / query does not exist      */
    GRAPH_ERR_INVALID_WEIGHT,     /* weight is NaN or infinite                      */
    GRAPH_ERR_NEGATIVE_WEIGHT,    /* algorithm needs non-negative weights           */
    GRAPH_ERR_CYCLE,              /* graph has a cycle but the algorithm needs a DAG */
    GRAPH_ERR_OUT_OF_MEMORY,
    GRAPH_ERR_UNSUPPORTED,        /* e.g. directed graph given to undirected-only   */
    GRAPH_ERR_NEGATIVE_CYCLE,
    GRAPH_ERR_TOO_LARGE,          /* would need more memory / time than allowed     */
    GRAPH_ERR_NOT_FOUND,          /* no path / tour / solution exists               */
    GRAPH_ERR_NO_CONVERGENCE,     /* iterative method hit its iteration limit       */
    GRAPH_ERR_IO,                 /* file read or write failed, or bad input        */
    GRAPH_ERR_NO_BACKEND_SUPPORT, /* the active representation lacks this operation */
    GRAPH_ERR_FROZEN,             /* the graph handle is frozen (read-only)         */
    GRAPH_ERR_COUNT               /* number of codes; not an error                  */
} GraphStatus;

/** Human-readable name of a code ("out of memory", ...); never NULL. */
const char *graph_status_string(GraphStatus status);

/** The most recent failure recorded on this thread (GRAPH_OK if none). */
GraphStatus graph_last_error(void);

/** Resets the recorded failure to GRAPH_OK. */
void graph_clear_error(void);

/** Records a failure and returns it, so callers can write `return graph_fail(c);`. */
GraphStatus graph_set_error(GraphStatus status);

#endif /* GRAPH_STATUS_H */
