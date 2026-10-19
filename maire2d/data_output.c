#include "data_output.h"

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

/*
 * Write mesh_t_<t>.csv in the current working directory.
 *
 * Node order: n = j * NUM_NODES_X + i (x index changes fastest).
 * coor_node[n][0] is x; coor_node[n][1] is y at time t.
 *
 * The caller must supply NUM_NODES > 0 and an array of that many nodes.
 * Use the default C numeric locale (decimal point '.').
 * Reusing a filename overwrites that file.
 * Errors are reported to stderr.
 */
void putOutNodeCoor(double t, size_t NUM_NODES_X, size_t NUM_NODES_Y,
                   size_t NUM_NODES, double coor_node[NUM_NODES][2])
{
    if (coor_node == NULL || !isfinite(t) ||
        NUM_NODES_X < 2 || NUM_NODES_Y < 2 ||
        NUM_NODES_X > SIZE_MAX / NUM_NODES_Y ||
        NUM_NODES != NUM_NODES_X * NUM_NODES_Y) {
        fprintf(stderr,
                "putOutNodeCoor: invalid time, dimensions, or array\n");
        return;
    }

    /* Invalid coordinates would break gnuplot's scan layout. */
    for (size_t n = 0; n < NUM_NODES; ++n) {
        if (!isfinite(coor_node[n][0]) ||
            !isfinite(coor_node[n][1])) {
            fprintf(stderr,
                    "putOutNodeCoor: nonfinite coordinate at node %zu\n",
                    n);
            return;
        }
    }

    char filename[128];
    int len = snprintf(filename, sizeof filename,
                       "build/mesh_t_%.17g.csv", t);

    if (len < 0 || (size_t)len >= sizeof filename) {
        fprintf(stderr,
                "putOutNodeCoor: cannot format output filename\n");
        return;
    }

    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
        perror(filename);
        return;
    }

    int ok = fprintf(fp, "t,node,i,j,x,y\n") >= 0;

    for (size_t j = 0; j < NUM_NODES_Y && ok; ++j) {
        for (size_t i = 0; i < NUM_NODES_X; ++i) {
            size_t n = j * NUM_NODES_X + i;

            if (fprintf(fp, "%.17g,%zu,%zu,%zu,%.17g,%.17g\n",
                        t, n, i, j,
                        coor_node[n][0], coor_node[n][1]) < 0) {
                ok = 0;
                break;
            }
        }

        /* Exactly one blank line separates successive mesh rows. */
        if (ok && j + 1 < NUM_NODES_Y && fputc('\n', fp) == EOF)
            ok = 0;
    }

    /* fclose also detects errors while flushing buffered output. */
    if (fclose(fp) != 0)
        ok = 0;

    if (!ok)
        fprintf(stderr,
                "putOutNodeCoor: could not finish writing %s\n",
                filename);
}
