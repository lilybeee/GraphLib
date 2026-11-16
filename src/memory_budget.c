/* memory_budget.c — see memory_budget.h. */

#include "memory_budget.h"

static size_t g_budget = GRAPH_DEFAULT_MEMORY_BUDGET;

size_t graph_memory_budget(void)
{
    return g_budget;
}

void graph_set_memory_budget(size_t bytes)
{
    g_budget = bytes;
}
