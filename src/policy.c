/*
 * policy.c — the observer.  See policy.h.
 *
 * The observer answers one question, "should this live graph move to another
 * representation?", from three inputs: how the graph has been used (workload
 * counters, decayed so old behaviour fades), what it looks like now (GraphShape), and
 * each candidate backend's cost model.  observer_evaluate() applies these rules in order:
 *
 *   1. Over the memory budget            -> FORCED_MEMORY: leave, whatever it costs.
 *   2. A migration happened very recently -> COOLDOWN: do not flip-flop on noise.
 *   3. Gain below the hysteresis ratio    -> NOT_ENOUGH_GAIN: a 10% win is not worth a rebuild.
 *   4. Saving over `horizon` windows does not repay the rebuild -> NOT_WORTH_IT.
 *   5. Otherwise                          -> MIGRATE.
 *
 * Only candidates that fit the budget (and, for the matrix, are within its size limit; the
 * facade filters those) are considered.  Nothing here changes the graph except
 * observer_maybe_migrate(), which rebuilds the graph in the target representation first
 * and only then destroys the old one, so a failed rebuild leaves the graph untouched.
 */

#include "policy.h"
#include "memory_budget.h"

#include <stdlib.h>

#include "migrate.h"

struct Observer {
    ObserverConfig cfg;
    Workload       workload;
    uint64_t       ticks;           /* graph operations since the last check */
    uint64_t       total_ops;
    uint64_t       last_migration;  /* total_ops at the last migration       */
    size_t         migrations;
};

ObserverConfig observer_default_config(void)
{
    ObserverConfig cfg = {
        .check_interval = 1000,
        .decay          = 0.5,
        .hysteresis     = 1.5,
        .horizon        = 8.0,
        .cooldown_ops   = 4000,
        .memory_budget  = 0,
        .log            = NULL,
    };
    return cfg;
}

Observer *observer_create(const ObserverConfig *cfg)
{
    Observer *obs = calloc(1, sizeof(*obs));
    if (obs == NULL) {
        return NULL;
    }
    obs->cfg = cfg != NULL ? *cfg : observer_default_config();
    if (obs->cfg.check_interval == 0) {
        obs->cfg.check_interval = 1;
    }
    return obs;
}

void observer_destroy(Observer *obs)
{
    free(obs);
}

void observer_note(Observer *obs, OpKind op)
{
    if (obs == NULL) {
        return;
    }
    workload_note(&obs->workload, op);
    if (op != OP_CACHE_HIT && op != OP_CACHE_MISS) {
        ++obs->ticks;
        ++obs->total_ops;
    }
}

bool observer_due(const Observer *obs)
{
    return obs != NULL && obs->ticks >= obs->cfg.check_interval;
}

size_t observer_memory_budget(const Observer *obs)
{
    return obs == NULL ? graph_memory_budget() : obs->cfg.memory_budget;
}

size_t observer_migrations(const Observer *obs)
{
    return obs == NULL ? 0 : obs->migrations;
}

const Workload *observer_workload(const Observer *obs)
{
    return obs == NULL ? NULL : &obs->workload;
}

GraphShape graph_shape_of(const GraphBackend *backend, const void *impl)
{
    GraphShape s;
    s.live     = backend->vertex_count(impl);
    s.id_bound = backend->id_bound(impl);
    s.edges    = backend->edge_count(impl);
    s.directed = backend->is_directed(impl);
    s.unweighted = backend->is_unweighted != NULL && backend->is_unweighted(impl);
    return s;
}

/* ------------------------------------------------------------------ */
/*  Decision                                                            */
/* ------------------------------------------------------------------ */

static bool has_cost_model(const GraphBackend *b)
{
    return b->estimate_cost != NULL && b->estimate_memory != NULL &&
           b->estimate_build_cost != NULL;
}

Decision observer_evaluate(const Observer *obs, const GraphBackend *current,
                           const void *impl,
                           const GraphBackend *const *candidates, size_t n)
{
    Decision d = { DECISION_STAY, NULL, 0.0, 0.0, 0.0 };
    if (obs == NULL || current == NULL || impl == NULL || !has_cost_model(current)) {
        return d;
    }

    const GraphShape shape = graph_shape_of(current, impl);
    const size_t     budget = obs->cfg.memory_budget;
    d.cost_current = current->estimate_cost(&shape, &obs->workload);

    /* Cheapest candidate that fits in memory (never the current one). */
    for (size_t i = 0; i < n; ++i) {
        const GraphBackend *c = candidates[i];
        if (c == NULL || c == current || !has_cost_model(c)) {
            continue;
        }
        if (budget != 0 && c->estimate_memory(&shape) > budget) {
            continue;
        }
        const double cost = c->estimate_cost(&shape, &obs->workload);
        if (d.target == NULL || cost < d.cost_target) {
            d.target      = c;
            d.cost_target = cost;
        }
    }
    if (d.target == NULL) {
        return d;
    }
    d.build_cost = d.target->estimate_build_cost(&shape);

    /* Over budget: leave regardless of what it costs. */
    if (budget != 0 && current->memory_bytes(impl) > budget) {
        d.kind = DECISION_FORCED_MEMORY;
        return d;
    }

    if (obs->migrations > 0 &&
        obs->total_ops - obs->last_migration < obs->cfg.cooldown_ops) {
        d.kind = d.cost_current > d.cost_target ? DECISION_COOLDOWN : DECISION_STAY;
        return d;
    }
    if (d.cost_current < obs->cfg.hysteresis * d.cost_target) {
        d.kind = d.cost_current > d.cost_target ? DECISION_NOT_ENOUGH_GAIN
                                                : DECISION_STAY;
        return d;
    }
    if ((d.cost_current - d.cost_target) * obs->cfg.horizon <= d.build_cost) {
        d.kind = DECISION_NOT_WORTH_IT;
        return d;
    }
    d.kind = DECISION_MIGRATE;
    return d;
}

bool observer_maybe_migrate(Observer *obs, const GraphBackend **current, void **impl,
                            const GraphBackend *const *candidates, size_t n,
                            Decision *out)
{
    Decision d = { DECISION_STAY, NULL, 0.0, 0.0, 0.0 };
    if (out != NULL) {
        *out = d;
    }
    if (obs == NULL || current == NULL || impl == NULL || !observer_due(obs)) {
        return false;
    }

    d = observer_evaluate(obs, *current, *impl, candidates, n);
    obs->ticks = 0;
    workload_decay(&obs->workload, obs->cfg.decay);
    if (out != NULL) {
        *out = d;
    }
    if (d.kind != DECISION_MIGRATE && d.kind != DECISION_FORCED_MEMORY) {
        return false;
    }

    void *fresh = graph_clone(*current, *impl, d.target);
    if (fresh == NULL) {
        return false; /* old graph untouched; try again at the next check */
    }

    if (obs->cfg.log != NULL) {
        const GraphShape s = graph_shape_of(d.target, fresh);
        fprintf(obs->cfg.log,
                "[observer] op %llu: %s -> %s (%s; predicted cost %.3g -> %.3g, "
                "rebuild %.3g; V=%zu E=%zu density=%.3f)\n",
                (unsigned long long)obs->total_ops, (*current)->name, d.target->name,
                d.kind == DECISION_FORCED_MEMORY ? "over memory budget" : "workload",
                d.cost_current, d.cost_target, d.build_cost, s.live, s.edges,
                shape_density(&s));
    }

    (*current)->destroy(*impl);
    *current = d.target;
    *impl    = fresh;
    ++obs->migrations;
    obs->last_migration = obs->total_ops;
    return true;
}
