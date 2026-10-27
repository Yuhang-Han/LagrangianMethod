#include "data_output.h"

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

/*
 * Write mesh_t_<t>.csv in the current working directory.
 *
 * x indices change fastest. CSV i/j are actual padded-array indices;
 * n is only a sequential output record number. coor_node[j][i][0] is x;
 * coor_node[j][i][1] is y at time t.
 *
 * Supply the padded array dimensions and inclusive physical-node bounds.
 * Use the default C numeric locale (decimal point '.').
 * Reusing a filename overwrites that file.
 * Errors are reported to stderr.
 */
void putOutNodeCoor(double t, size_t NUM_GHOST_NODES_X, size_t NUM_GHOST_NODES_Y,
                   size_t START_NODE_X, size_t END_NODE_X,
                   size_t START_NODE_Y, size_t END_NODE_Y,
                   double coor_node[NUM_GHOST_NODES_Y][NUM_GHOST_NODES_X][2])
{
    if (coor_node == NULL || !isfinite(t) ||
        NUM_GHOST_NODES_X < 2 || NUM_GHOST_NODES_Y < 2 ||
        NUM_GHOST_NODES_X > SIZE_MAX / NUM_GHOST_NODES_Y ||
        START_NODE_X >= END_NODE_X || START_NODE_Y >= END_NODE_Y ||
        END_NODE_X >= NUM_GHOST_NODES_X || END_NODE_Y >= NUM_GHOST_NODES_Y) {
        fprintf(stderr,
                "putOutNodeCoor: invalid time, dimensions, or array\n");
        return;
    }

    /* Invalid coordinates would break gnuplot's scan layout. */
    size_t n = 0;
    for (size_t j = START_NODE_Y; j <= END_NODE_Y; ++j) {
        for (size_t i = START_NODE_X; i <= END_NODE_X; ++i) {
            if (!isfinite(coor_node[j][i][0]) ||
                !isfinite(coor_node[j][i][1])) {
                fprintf(stderr,
                        "putOutNodeCoor: nonfinite coordinate at node %zu\n",
                        n);
                return;
            }
            ++n;
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

    n = 0;
    for (size_t j = START_NODE_Y; j <= END_NODE_Y && ok; ++j) {
        for (size_t i = START_NODE_X; i <= END_NODE_X; ++i) {
            if (fprintf(fp, "%.17g,%zu,%zu,%zu,%.17g,%.17g\n",
                        t, n, i, j,
                        coor_node[j][i][0], coor_node[j][i][1]) < 0) {
                ok = 0;
                break;
            }
            ++n;
        }

        /* Exactly one blank line separates successive mesh rows. */
        if (ok && j < END_NODE_Y && fputc('\n', fp) == EOF)
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
