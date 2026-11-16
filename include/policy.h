/*
 * policy.h — the observer: decides when to change a graph's representation.
 *
 * The observer is passive and runs on the caller's thread.  The facade calls
 * observer_note() for every operation it performs; once enough operations have
 * accumulated, observer_maybe_migrate() asks each candidate backend what the
 * recent workload would cost on it, and migrates if the answer is clearly
 * better.  There is no background thread, so decisions are deterministic and
 * need no locking.
 *
 * A migration happens only if all of these hold:
 *
 *   1. Memory     the target's predicted footprint fits the memory budget
 *                 (a hard constraint, checked first);
 *   2. Hysteresis current cost >= hysteresis * target cost, so a graph sitting
 *                 near the break-even point does not flip back and forth;
 *   3. Payback    the saving over the next `horizon` windows exceeds the
 *                 predicted cost of rebuilding the graph;
 *   4. Cooldown   enough operations have passed since the last migration.
 *
 * If the current representation itself exceeds the memory budget, the
 * observer migrates to the cheapest candidate that fits, ignoring 2-4.
 */

#ifndef POLICY_H
#define POLICY_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include "backend.h"
#include "stats.h"

typedef struct {
    size_t   check_interval;  /* graph operations between evaluations          */
    double   decay;           /* counters are multiplied by this after a check */
    double   hysteresis;      /* required cost ratio current/target, > 1       */
    double   horizon;         /* windows over which a saving must repay a build */
    uint64_t cooldown_ops;    /* minimum operations between two migrations     */
    size_t   memory_budget;   /* bytes; 0 means unlimited                      */
    FILE    *log;             /* decisions that lead to a migration; or NULL   */
} ObserverConfig;

/** check_interval 1000, decay 0.5, hysteresis 1.5, horizon 8, cooldown 4000. */
ObserverConfig observer_default_config(void);

typedef enum {
    DECISION_STAY,           /* nothing better, or too early to tell          */
    DECISION_COOLDOWN,       /* a better backend exists but we just migrated  */
    DECISION_NOT_ENOUGH_GAIN,/* better, but not by the hysteresis margin      */
    DECISION_NOT_WORTH_IT,   /* big enough gain, but the rebuild won't repay  */
    DECISION_MIGRATE,        /* migrate: clear gain                           */
    DECISION_FORCED_MEMORY   /* migrate: current backend is over budget       */
} DecisionKind;

typedef struct {
    DecisionKind        kind;
    const GraphBackend *target;        /* best candidate, or NULL if none fits */
    double              cost_current;  /* predicted cost on the current backend */
    double              cost_target;
    double              build_cost;    /* predicted cost of migrating           */
} Decision;

typedef struct Observer Observer;

/** cfg may be NULL for the defaults.  Returns NULL on allocation failure. */
Observer *observer_create(const ObserverConfig *cfg);
void      observer_destroy(Observer *obs);

/** The memory budget this observer enforces (bytes; 0 = unlimited; the library default if obs is NULL). */
size_t observer_memory_budget(const Observer *obs);

/** Record one operation.  Cache hits and misses are recorded but do not
 *  advance the evaluation clock (a hit does no work on the graph). */
void observer_note(Observer *obs, OpKind op);

/** True once check_interval graph operations have passed since the last check. */
bool observer_due(const Observer *obs);

/** The shape of a graph held by 'backend'. */
GraphShape graph_shape_of(const GraphBackend *backend, const void *impl);

/**
 * observer_evaluate - decide what to do now, without changing anything.
 * @candidates: every backend that may be used, including the current one.
 */
Decision observer_evaluate(const Observer *obs, const GraphBackend *current,
                           const void *impl,
                           const GraphBackend *const *candidates, size_t n);

/**
 * observer_maybe_migrate - the call the facade makes after each operation.
 * Does nothing (returns false) until observer_due().  When due, evaluates,
 * decays the counters, and if the decision is to migrate, rebuilds the graph
 * in the target representation, destroys the old one, and updates *current
 * and *impl.  If the rebuild fails the old graph is left untouched.
 * @out: optional; receives the decision taken.
 * Returns true iff the graph was migrated.
 */
bool observer_maybe_migrate(Observer *obs, const GraphBackend **current, void **impl,
                            const GraphBackend *const *candidates, size_t n,
                            Decision *out);

size_t         observer_migrations(const Observer *obs);
const Workload *observer_workload(const Observer *obs);

#endif /* POLICY_H */
