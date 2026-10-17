// FOR QUAD MESH ONLY

#ifndef MESH_CALCULATION_H
#define MESH_CALCULATION_H

#include <stdio.h>
#include "linear_algebra_2d.h"


void categorizeNodes(const size_t NUM_GHOST_NODES_X, const size_t NUM_GHOST_NODES_Y, 
					 const size_t NUM_NODES, size_t NODE_DIC[NUM_NODES]);
					 
void categorizeCells(const size_t NUM_GHOST_CELLS_X, const size_t NUM_GHOST_CELLS_Y, 
					 const size_t NUM_CELLS, const size_t NUM_BOUND_CELLS,
					 size_t CELL_DIC[NUM_CELLS], size_t BOUND_CELL_DIC[NUM_BOUND_CELLS]);


// WARNING: These two function does not contain coundary safty check.
size_t localToGlobalCellIdx(size_t iP, size_t iCP, 
							size_t NUM_NODES_X, size_t NUM_CELLS_X);
size_t localToGlobalNodeIdx(size_t iC, size_t iPC, 
							size_t NUM_NODES_X, size_t NUM_CELLS_X);
size_t globalToLocalNodeIdx(size_t iC, size_t iP, 
							size_t NUM_NODES_X, size_t NUM_CELLS_X);

size_t calNeighVertexIdx(size_t iP, size_t iC, int direction,
						 size_t NUM_NODES_X, size_t NUM_CELLS_X);

double calCornerLength(double coor_p[2], double coor_pl[2], double coor_pr[2]);
void   calCornerNormal(double coor_p[2], double coor_pl[2], double coor_pr[2],
					   double normal_corner[2]);

#endif
