/*
 * graph_results.c — destructors for the shared result types (graph_results.h).
 */

#include "graph_results.h"

#include <stdlib.h>
#include <string.h>

void graph_groups_free(GraphGroups *groups)
{
    if (groups == NULL) {
        return;
    }
    free(groups->offsets);
    free(groups->members);
    memset(groups, 0, sizeof(*groups));
}

void graph_path_free(GraphPath *path)
{
    if (path == NULL) {
        return;
    }
    free(path->vertices);
    memset(path, 0, sizeof(*path));
}
