/*
 * stub_backend.h — a deliberately simple dense-matrix GraphBackend.
 *
 * Test support only.  It stands in for the real adjacency-matrix backend
 * (which lives on the Adjacency-Matrix branch) so the observer, migration and
 * cache can be exercised on this branch.  Correctness over speed: every
 * operation is the obvious one.  It implements only the required vtable
 * slots (no algorithm kernels).
 *
 * Its cost model has two modes:
 *   STUB_COST_MATRIX  roughly what a matrix costs: O(1) point queries,
 *                     O(V) neighbour scans, O(V^2) whole-graph algorithms,
 *                     O(V^2) memory.
 *   STUB_COST_SCALED  `scale` times the list's cost model, so a test can put
 *                     the stub at an exact ratio from the list (e.g. 1.3x
 *                     cheaper) and probe the hysteresis boundary.
 */

#ifndef STUB_BACKEND_H
#define STUB_BACKEND_H

#include "backend.h"

typedef enum { STUB_COST_MATRIX, STUB_COST_SCALED } StubCostMode;

const GraphBackend *stub_backend(void);

/** Select the cost model; `scale` is used only by STUB_COST_SCALED. */
void stub_set_cost_mode(StubCostMode mode, double scale);

/** Make create() fail, to test that a failed migration leaves the graph alone. */
void stub_set_fail_create(bool fail);

#endif /* STUB_BACKEND_H */
