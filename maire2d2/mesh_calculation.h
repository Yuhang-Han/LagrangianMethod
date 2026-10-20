// FOR QUAD MESH ONLY

#ifndef MESH_CALCULATION_H
#define MESH_CALCULATION_H

#include <stdio.h>
#include "linear_algebra_2d.h"


/*
 * Return actual array indices of the next (+1) or previous (-1) CCW vertex.
 * The node must belong to cell (cell_x, cell_y), whose southwest node has
 * the same x/y indices. Output pointers must be valid; there are no bounds
 * checks, as in the original helper.
 */
void calNeighVertexIdx(size_t node_x, size_t node_y,
                      size_t cell_x, size_t cell_y, int direction,
                      size_t *neigh_node_x, size_t *neigh_node_y);

double calCellVolume(const double coor_sw[2], const double coor_se[2],
                     const double coor_ne[2], const double coor_nw[2]);

double calCornerLength(double coor_p[2], double coor_pl[2], double coor_pr[2]);
void   calCornerNormal(double coor_p[2], double coor_pl[2], double coor_pr[2],
					   double normal_corner[2]);

#endif
