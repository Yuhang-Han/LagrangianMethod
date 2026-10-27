#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

void putOutNodeCoor(double t, size_t NUM_GHOST_NODES_X, size_t NUM_GHOST_NODES_Y,
                   size_t START_NODE_X, size_t END_NODE_X,
                   size_t START_NODE_Y, size_t END_NODE_Y,
                   double coor_node[NUM_GHOST_NODES_Y][NUM_GHOST_NODES_X][2]);
