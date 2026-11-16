/*
 * graph_crossover.c — where does the matrix start beating the list?
 *
 * graph_crossover_report() builds random graphs over a grid of sizes and
 * densities, times every algorithm on both representations with
 * graph_compare(), prints the raw timings and, for each algorithm and size,
 * the density at which the faster representation changes.
 */

#include <math.h>
#include <stdlib.h>

#include "graph_algorithms.h"
#include "graph_internal.h"

typedef struct {
    double list_ns;
    double matrix_ns;
    size_t edges;
    bool   valid;       /* both representations produced a result */
    bool   agreed;
} Cell;

static const char *winner_name(const Cell *c)
{
    return c->matrix_ns < c->list_ns ? "matrix" : "list";
}

bool graph_crossover_report(FILE *out, const size_t *sizes, size_t n_sizes,
                            const double *densities, size_t n_dens,
                            unsigned long long seed, int repetitions, bool csv)
{
    if (out == NULL || sizes == NULL || densities == NULL || n_sizes == 0 || n_dens == 0) {
        graph_set_error(GRAPH_ERR_INVALID_ARGUMENT);
        return false;
    }
    Cell *cells = calloc(n_sizes * n_dens * GRAPH_ALGO_COUNT, sizeof(*cells));
    if (cells == NULL) {
        graph_set_error(GRAPH_ERR_OUT_OF_MEMORY);
        return false;
    }
#define CELL(si, di, a) cells[((si) * n_dens + (di)) * GRAPH_ALGO_COUNT + (size_t)(a)]

    for (size_t si = 0; si < n_sizes; ++si) {
        for (size_t di = 0; di < n_dens; ++di) {
            Graph *undirected = graph_gen_erdos_renyi(GRAPH_LIST, sizes[si], densities[di],
                                                      seed, false, 1.0, 100.0);
            Graph *dag = graph_gen_random_dag(GRAPH_LIST, sizes[si], densities[di], seed, 1.0, 100.0);
            if (undirected == NULL || dag == NULL) {
                graph_destroy(undirected);
                graph_destroy(dag);
                free(cells);
                return false;       /* graph_last_error() already set */
            }
            for (int a = 0; a < GRAPH_ALGO_COUNT; ++a) {
                GraphCompareResult r;
                const Graph *g = a == GRAPH_ALGO_TOPO ? dag : undirected;
                Cell *c = &CELL(si, di, a);
                if (graph_compare(g, (GraphAlgo)a, 0, repetitions, &r)) {
                    c->list_ns   = r.list_ns;
                    c->matrix_ns = r.matrix_ns;
                    c->edges     = graph_edge_count(g);
                    c->valid     = r.list_ok && r.matrix_ok;
                    c->agreed    = r.agreed;
                }
            }
            graph_destroy(undirected);
            graph_destroy(dag);
        }
    }

    if (csv) {
        fputs("algorithm,n,density,edges,list_ns,matrix_ns,faster\n", out);
    } else {
        fprintf(out, "%-22s %7s %9s %9s %13s %13s  %s\n", "algorithm", "n", "density",
                "edges", "list ns", "matrix ns", "faster");
    }
    for (int a = 0; a < GRAPH_ALGO_COUNT; ++a) {
        for (size_t si = 0; si < n_sizes; ++si) {
            for (size_t di = 0; di < n_dens; ++di) {
                const Cell *c = &CELL(si, di, a);
                if (!c->valid) {
                    continue;
                }
                if (csv) {
                    fprintf(out, "%s,%zu,%g,%zu,%.0f,%.0f,%s\n", graph_algo_name((GraphAlgo)a),
                            sizes[si], densities[di], c->edges, c->list_ns, c->matrix_ns,
                            winner_name(c));
                } else {
                    fprintf(out, "%-22s %7zu %9.4f %9zu %13.0f %13.0f  %s (%.2fx)%s\n",
                            graph_algo_name((GraphAlgo)a), sizes[si], densities[di], c->edges,
                            c->list_ns, c->matrix_ns, winner_name(c),
                            c->list_ns > c->matrix_ns ? c->list_ns / fmax(c->matrix_ns, 1.0)
                                                      : c->matrix_ns / fmax(c->list_ns, 1.0),
                            c->agreed ? "" : "  RESULTS DIFFER");
                }
            }
        }
    }

    if (!csv) {
        fputs("\nCrossover summary (density at which the faster representation changes):\n", out);
        for (int a = 0; a < GRAPH_ALGO_COUNT; ++a) {
            for (size_t si = 0; si < n_sizes; ++si) {
                const Cell *first = NULL, *last = NULL;
                bool found = false;
                double flip_lo = 0.0, flip_hi = 0.0;
                size_t valid = 0;
                const Cell *prev = NULL;
                double prev_d = 0.0;
                for (size_t di = 0; di < n_dens; ++di) {
                    const Cell *c = &CELL(si, di, a);
                    if (!c->valid) {
                        continue;
                    }
                    ++valid;
                    if (first == NULL) {
                        first = c;
                    }
                    if (prev != NULL && !found && winner_name(prev)[0] != winner_name(c)[0]) {
                        found   = true;
                        flip_lo = prev_d;
                        flip_hi = densities[di];
                    }
                    prev = c;
                    prev_d = densities[di];
                    last = c;
                }
                if (valid == 0) {
                    continue;
                }
                fprintf(out, "  %-22s n=%-6zu ", graph_algo_name((GraphAlgo)a), sizes[si]);
                if (found) {
                    fprintf(out, "%s -> %s between density %g and %g\n", winner_name(first),
                            winner_name(last), flip_lo, flip_hi);
                } else {
                    fprintf(out, "%s at every density tried\n", winner_name(first));
                }
            }
        }
    }
    free(cells);
    return true;
#undef CELL
}
