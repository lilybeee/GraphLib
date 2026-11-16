/*
 * list_backend.h — the adjacency list exposed as a GraphBackend.
 */

#ifndef LIST_BACKEND_H
#define LIST_BACKEND_H

#include "backend.h"

/** Returns the (static, immutable) vtable for the adjacency-list backend. */
const GraphBackend *list_backend(void);

#endif /* LIST_BACKEND_H */
