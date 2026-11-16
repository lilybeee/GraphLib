/*
 * memory_budget.h — the library-wide memory budget.
 *
 * Every graph remembers, when it is created, how many bytes it may use.  A
 * mutation that would take it past that (adding an edge or vertex, a conversion
 * to the other representation, loading or generating a big graph) is refused
 * with GRAPH_ERR_TOO_LARGE, the graph stays exactly as it was, and the program
 * keeps running.  Without this, a request that needs more memory than the
 * machine has ends in the operating system killing the process.
 *
 * The default is 4 GiB, which a typical machine can give a single graph.
 * Change it once at start-up, before creating graphs; a graph keeps the budget
 * that was in force when it was created.  0 means unlimited.
 */

#ifndef MEMORY_BUDGET_H
#define MEMORY_BUDGET_H

#include <stddef.h>

#define GRAPH_DEFAULT_MEMORY_BUDGET ((size_t)4 << 30)   /* 4 GiB */

/** The budget new graphs get (bytes; 0 = unlimited). */
size_t graph_memory_budget(void);

/** Set the budget for graphs created from now on.  Not thread-safe: call it
 *  during start-up, before other threads create graphs. */
void graph_set_memory_budget(size_t bytes);

#endif /* MEMORY_BUDGET_H */
