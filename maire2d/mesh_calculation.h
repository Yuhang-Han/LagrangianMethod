// FOR QUAD MESH ONLY

#ifndef MESH_CALCULATION_H
#define MESH_CALCULATION_H

#include <stdio.h>
#include "linear_algebra_2d.h"

// WARNING: These two function does not contain coundary safty check.
size_t localToGlobalCellIdx(size_t idxP, size_t idxCP, 
							size_t NUM_NODES_X, size_t NUM_CELLS_X);
size_t localToGlobalNodeIdx(size_t idxC, size_t idxPC, 
							size_t NUM_NODES_X, size_t NUM_CELLS_X);
size_t globalToLocalNodeIdx(size_t idxC, size_t idxP, 
							size_t NUM_NODES_X, size_t NUM_CELLS_X);

size_t calNeighVertexIdx(size_t idxP, size_t idxC, int direction,
						 size_t NUM_NODES_X, size_t NUM_CELLS_X);

double calCornerLength(double coor_p[2], double coor_pl[2], double coor_pr[2]);
void   calCornerNormal(double coor_p[2], double coor_pl[2], double coor_pr[2],
					   double normal_corner[2]);

#endif
