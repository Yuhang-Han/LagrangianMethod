#include "data_output.h"

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

/*
 * x varies fastest. Only the inclusive physical-node window is written.
 * CSV i/j preserve the actual node_x/node_y array indices, including their
 * ghost-layer offset. The node column is a sequential output record number,
 * never a solver index. Blank lines preserve the gnuplot mesh-row scans.
 *
 * The build directory must exist. Use the default C numeric locale ('.').
 * Reusing a time overwrites that output file. Errors are reported to stderr.
 */
int putOutNodeCoor(
    double t, size_t NUM_TOTAL_NODES_X, size_t NUM_TOTAL_NODES_Y,
    size_t START_NODE_X, size_t END_NODE_X,
    size_t START_NODE_Y, size_t END_NODE_Y,
    double coor_node[NUM_TOTAL_NODES_Y][NUM_TOTAL_NODES_X][2])
{
    if (coor_node == NULL || !isfinite(t) || NUM_TOTAL_NODES_X < 2 ||
        NUM_TOTAL_NODES_Y < 2 ||
        NUM_TOTAL_NODES_X > SIZE_MAX / NUM_TOTAL_NODES_Y ||
        START_NODE_X >= END_NODE_X || START_NODE_Y >= END_NODE_Y ||
        END_NODE_X >= NUM_TOTAL_NODES_X || END_NODE_Y >= NUM_TOTAL_NODES_Y) {
        fprintf(stderr,
                "putOutNodeCoor: invalid time, bounds, or array\n");
        return 0;
    }

    /* Ghost coordinates are outside the output window and need not be set. */
    for (size_t node_y = START_NODE_Y; node_y <= END_NODE_Y; ++node_y) {
        for (size_t node_x = START_NODE_X; node_x <= END_NODE_X; ++node_x) {
            if (!isfinite(coor_node[node_y][node_x][0]) ||
                !isfinite(coor_node[node_y][node_x][1])) {
                fprintf(stderr,
                        "putOutNodeCoor: nonfinite coordinate at node (%zu,%zu)\n",
                        node_x, node_y);
                return 0;
            }
        }
    }

    char filename[128];
    int len = snprintf(filename, sizeof filename,
                       "build/mesh_t_%.17g.csv", t);

    if (len < 0 || (size_t)len >= sizeof filename) {
        fprintf(stderr,
                "putOutNodeCoor: cannot format output filename\n");
        return 0;
    }

    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
        perror(filename);
        return 0;
    }

    int ok = fprintf(fp, "t,node,i,j,x,y\n") >= 0;

    size_t node_record = 0;
    for (size_t node_y = START_NODE_Y; node_y <= END_NODE_Y && ok; ++node_y) {
        for (size_t node_x = START_NODE_X; node_x <= END_NODE_X; ++node_x) {
            if (fprintf(fp, "%.17g,%zu,%zu,%zu,%.17g,%.17g\n",
                        t, node_record, node_x, node_y,
                        coor_node[node_y][node_x][0],
                        coor_node[node_y][node_x][1]) < 0) {
                ok = 0;
                break;
            }
            ++node_record;
        }

        /* Exactly one blank line separates successive mesh rows. */
        if (ok && node_y < END_NODE_Y && fputc('\n', fp) == EOF)
            ok = 0;
    }

    /* fclose also detects errors while flushing buffered output. */
    if (fclose(fp) != 0)
        ok = 0;

    if (!ok)
        fprintf(stderr,
                "putOutNodeCoor: could not finish writing %s\n",
                filename);
    return ok;
}
