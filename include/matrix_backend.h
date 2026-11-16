/*
 * matrix_backend.h — the adjacency matrix exposed as a GraphBackend.
 */

#ifndef MATRIX_BACKEND_H
#define MATRIX_BACKEND_H

#include "backend.h"

/** Returns the (static, immutable) vtable for the adjacency-matrix backend. */
const GraphBackend *matrix_backend(void);

#endif /* MATRIX_BACKEND_H */
