/*
 * stats.c — workload counters.  See stats.h.
 *
 * The counters are doubles, not integers, so that workload_decay() can multiply them by a
 * factor below 1: recent behaviour then outweighs old behaviour smoothly, with no hard
 * "window" edge where the picture suddenly changes.
 */

#include "stats.h"

void workload_note(Workload *w, OpKind op)
{
    if (w != NULL && (int)op >= 0 && op < OP_COUNT) {
        w->n[op] += 1.0;
    }
}

void workload_decay(Workload *w, double factor)
{
    if (w == NULL) {
        return;
    }
    for (int i = 0; i < OP_COUNT; ++i) {
        w->n[i] *= factor;
    }
}

void workload_reset(Workload *w)
{
    if (w == NULL) {
        return;
    }
    for (int i = 0; i < OP_COUNT; ++i) {
        w->n[i] = 0.0; /* assigned, not scaled: 0 * NaN would stay NaN */
    }
}

double workload_graph_ops(const Workload *w)
{
    double total = 0.0;
    for (int i = 0; i < OP_COUNT; ++i) {
        if (i != OP_CACHE_HIT && i != OP_CACHE_MISS) {
            total += w->n[i];
        }
    }
    return total;
}

double workload_hit_rate(const Workload *w)
{
    const double lookups = w->n[OP_CACHE_HIT] + w->n[OP_CACHE_MISS];
    return lookups == 0.0 ? 0.0 : w->n[OP_CACHE_HIT] / lookups;
}
