/*
 * crossover_report — prints, for each algorithm, where the matrix starts to beat
 * the list.  Usage: crossover_report [--csv] [n1 n2 ...]   (default sizes 200 1000 2000)
 * Build with -O2 for real numbers:  make crossover CFLAGS="-std=c11 -O2"
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "graph.h"

int main(int argc, char **argv)
{
    bool csv = false;
    size_t sizes[16];
    size_t n_sizes = 0;
    for (int i = 1; i < argc && n_sizes < 16; ++i) {
        if (strcmp(argv[i], "--csv") == 0) {
            csv = true;
        } else {
            sizes[n_sizes++] = (size_t)strtoul(argv[i], NULL, 10);
        }
    }
    if (n_sizes == 0) {
        sizes[0] = 200; sizes[1] = 1000; sizes[2] = 2000;
        n_sizes = 3;
    }
    static const double densities[] = { 0.001, 0.005, 0.02, 0.05, 0.1, 0.3 };
    if (!graph_crossover_report(stdout, sizes, n_sizes, densities,
                                sizeof(densities) / sizeof(densities[0]), 12345, 3, csv)) {
        fprintf(stderr, "crossover report failed: %s\n", graph_status_string(graph_last_error()));
        return 1;
    }
    return 0;
}
