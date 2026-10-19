#ifndef DATA_OUTPUT_H
#define DATA_OUTPUT_H

#include <stddef.h>

/*
 * Write physical nodes within the inclusive bounds to build/mesh_t_<t>.csv.
 * CSV i and j are actual array indices; node is only a sequential record
 * number. The physical window has at least two nodes in each direction.
 * The caller supplies positive dimensions and matching storage.
 * Returns 1 on success, 0 on invalid input or an I/O error.
 */
int putOutNodeCoor(
    double t, size_t NUM_TOTAL_NODES_X, size_t NUM_TOTAL_NODES_Y,
    size_t START_NODE_X, size_t END_NODE_X,
    size_t START_NODE_Y, size_t END_NODE_Y,
    double coor_node[NUM_TOTAL_NODES_Y][NUM_TOTAL_NODES_X][2]);

#endif
