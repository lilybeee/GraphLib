/*
 * graph_status.c — error names and the "last error" slot.
 *
 * Functions that return a GraphStatus hand the code back directly; functions that return a
 * bool, size_t or double cannot, so they record the reason here instead.  The slot is
 * _Thread_local so that two threads using separate graphs never overwrite each other's reason,
 * and a success deliberately does NOT clear it (clear it yourself with graph_clear_error):
 * that way the reason survives the cleanup calls that often follow a failure.
 */

#include "graph_status.h"

static _Thread_local GraphStatus last_error = GRAPH_OK;

static const char *const names[GRAPH_ERR_COUNT] = {
    [GRAPH_OK]                     = "ok",
    [GRAPH_ERR_INVALID_ARGUMENT]   = "invalid argument",
    [GRAPH_ERR_INVALID_VERTEX]     = "invalid or removed vertex",
    [GRAPH_ERR_DUPLICATE_EDGE]     = "duplicate edge",
    [GRAPH_ERR_EDGE_NOT_FOUND]     = "edge not found",
    [GRAPH_ERR_INVALID_WEIGHT]     = "invalid weight (not finite, or not 1.0 in an unweighted graph)",
    [GRAPH_ERR_NEGATIVE_WEIGHT]    = "negative weight not allowed",
    [GRAPH_ERR_CYCLE]              = "graph contains a cycle",
    [GRAPH_ERR_OUT_OF_MEMORY]      = "out of memory",
    [GRAPH_ERR_UNSUPPORTED]        = "unsupported for this graph",
    [GRAPH_ERR_NEGATIVE_CYCLE]     = "negative cycle",
    [GRAPH_ERR_TOO_LARGE]          = "input too large",
    [GRAPH_ERR_NOT_FOUND]          = "no such path or solution",
    [GRAPH_ERR_NO_CONVERGENCE]     = "did not converge",
    [GRAPH_ERR_IO]                 = "I/O or format error",
    [GRAPH_ERR_NO_BACKEND_SUPPORT] = "operation not supported by this representation",
    [GRAPH_ERR_FROZEN]             = "graph is frozen (read-only)",
};

const char *graph_status_string(GraphStatus status)
{
    if ((unsigned)status >= (unsigned)GRAPH_ERR_COUNT) {
        return "unknown error";
    }
    return names[status];
}

GraphStatus graph_last_error(void)
{
    return last_error;
}

void graph_clear_error(void)
{
    last_error = GRAPH_OK;
}

GraphStatus graph_set_error(GraphStatus status)
{
    last_error = status;
    return status;
}
