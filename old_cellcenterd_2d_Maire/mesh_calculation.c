#include <stdio.h>

#include "linear_algebra_2d.h"
#include "mesh_calculation.h"

void categorizeNodes(const size_t NUM_GHOST_NODES_X, const size_t NUM_GHOST_NODES_Y, 
					 const size_t NUM_NODES, const size_t NUM_BOUND_NODES,
					 size_t NODE_DIC[NUM_NODES], size_t BOUND_NODE_DIC[NUM_BOUND_NODES])
{
	size_t num_ghost_nodes = NUM_GHOST_NODES_X * NUM_GHOST_NODES_Y;
	size_t iNode = 0;
	size_t iBound = 0;
	
	for(size_t iP = 0; iP < num_ghost_nodes; ++iP) {
		size_t row_node = iP / NUM_GHOST_NODES_X;
		size_t col_node	= iP % NUM_GHOST_NODES_Y;

		if ( (row_node == 0) || (col_node == 0)
			|| (row_node == NUM_GHOST_NODES_Y - 1)
			|| (col_node == NUM_GHOST_NODES_X - 1) )
		{
			//printf("iC: %zu (%zu, %zu) is %zu th bound node\n", 
			//		iC, row_node, col_node, iBound);
			BOUND_NODE_DIC[iBound] = iP;
			++iBound;
		}
		else
		{
			NODE_DIC[iNode] = iP;
			++iNode;
		}
	}
}



void categorizeCells(const size_t NUM_GHOST_CELLS_X, const size_t NUM_GHOST_CELLS_Y, 
					 const size_t NUM_CELLS, const size_t NUM_BOUND_CELLS,
					 size_t CELL_DIC[NUM_CELLS], size_t BOUND_CELL_DIC[NUM_BOUND_CELLS])
{
	size_t num_ghost_cells = NUM_GHOST_CELLS_X * NUM_GHOST_CELLS_Y;
	size_t iCell = 0;
	size_t iBound = 0;

	for(size_t iC = 0; iC < num_ghost_cells; ++iC) {
		size_t row_cell = iC / NUM_GHOST_CELLS_X;
		size_t col_cell = iC % NUM_GHOST_CELLS_Y;

		if ( (row_cell == 0) || (col_cell == 0)
			|| (row_cell == NUM_GHOST_CELLS_Y - 1)
			|| (col_cell == NUM_GHOST_CELLS_X - 1) )
		{
			//printf("iC: %zu (%zu, %zu) is %zu th bound cell\n", 
			//		iC, row_cell, col_cell, iBound);
			BOUND_CELL_DIC[iBound] = iC;
			++iBound;
		}
		else
		{
			CELL_DIC[iCell] = iC;
			++iCell;
		}
	}

}


// WARNING: This function does not contain coundary safty check.
size_t localToGlobalCellIdx(size_t iP, size_t iCP, 
							size_t NUM_NODES_X, size_t NUM_CELLS_X) 
{
	//printf("iP, iCP: %zu, %zu\n", iP, iCP);
	
	// locate the node P
	size_t row_node = iP / NUM_NODES_X;
	size_t col_node = iP % NUM_NODES_X;

	//printf("row_node, col_node: %zu, %zu\n", row_node, col_node);

	// left-down neigh cell
	size_t row_base = row_node - 1;
	size_t col_base = col_node - 1;

	//printf("row_base, col_base: %zu, %zu\n", row_base, col_base);

	size_t row_cell = row_base + iCP / 2;
	size_t col_cell = col_base + iCP % 2;

	//printf("row_cell, col_cell: %zu, %zu\n", row_cell, col_cell);

	size_t iC = row_cell * NUM_CELLS_X + col_cell;

	return iC;
}


size_t localToGlobalNodeIdx(size_t iC, size_t iPC, 
							size_t NUM_NODES_X, size_t NUM_CELLS_X) 
{
	//printf("iC, iPC: %zu, %zu\n", iC, iPC);
	
	// locate the cell C
	size_t row_cell = iC / NUM_CELLS_X;
	size_t col_cell = iC % NUM_CELLS_X;

	//printf("row_cell, col_cell: %zu, %zu\n", row_cell, col_cell);

	// left-down neigh vertex
	size_t row_base = row_cell;
	size_t col_base = col_cell;

	//printf("row_base, col_base: %zu, %zu\n", row_base, col_base);

	size_t row_node = row_base + iPC / 2;
	size_t col_node = col_base + iPC % 2;

	//printf("row_node, col_node: %zu, %zu\n", row_node, col_node);

	size_t iP = row_node * NUM_NODES_X + col_node;

	return iP;
}


size_t globalToLocalNodeIdx(size_t iC, size_t iP, 
							size_t NUM_NODES_X, size_t NUM_CELLS_X) 
{
	//printf("iC, iPC: %zu, %zu\n", iC, iPC);
	
	// locate the cell C
	size_t row_cell = iC / NUM_CELLS_X;
	size_t col_cell = iC % NUM_CELLS_X;

	// locate the node P
	size_t row_node = iP / NUM_NODES_X;
	size_t col_node = iP % NUM_NODES_X;

	//printf("row_cell, col_cell: %zu, %zu\n", row_cell, col_cell);


	//printf("row_base, col_base: %zu, %zu\n", row_base, col_base);

	size_t row_vertex = row_node - row_cell;
	size_t col_vertex = col_node - col_cell;

	//printf("row_node, col_node: %zu, %zu\n", row_node, col_node);

	size_t iPC = row_vertex * 2 + col_vertex;

	return iPC;
}


size_t calNeighVertexIdx(size_t iP, size_t iC, int direction,
						 size_t NUM_NODES_X, size_t NUM_CELLS_X)
{
	size_t iPC = globalToLocalNodeIdx(iC, iP, NUM_NODES_X, NUM_CELLS_X);

	size_t vertex_book[4] = {0, 1, 3, 2};

	int order_node = 0;
	for (int i = 0; i < 4; ++i) {
		if (vertex_book[i] == iPC)
			order_node = i;
	}

	size_t iPC_neigh = vertex_book[(order_node + direction + 4) % 4];

	return localToGlobalNodeIdx(iC, iPC_neigh, NUM_NODES_X, NUM_CELLS_X);
}


double calCornerLength(double coor_p[2], double coor_pl[2], double coor_pr[2])
{
	double vec_plp[2], vec_ppr[2];
	vecSubVec2(coor_p, coor_pl, vec_plp);
	vecSubVec2(coor_pr, coor_p, vec_ppr);

	double length_plp = calNorm2(vec_plp);
	double length_ppr = calNorm2(vec_ppr);

	double normal_plp[2], normal_ppr[2];
	calRHNormal2(vec_plp, normal_plp);
	calRHNormal2(vec_ppr, normal_ppr);

	double half_directed_normal_plp[2], half_directed_normal_ppr[2];
	numMulVec2(length_plp/2.0, normal_plp, half_directed_normal_plp);
	numMulVec2(length_ppr/2.0, normal_ppr, half_directed_normal_ppr);

	double directed_length_corner[2];
	vecAddVec2(half_directed_normal_plp, half_directed_normal_ppr,
				directed_length_corner);

	double length_corner = calNorm2(directed_length_corner);

	return length_corner;
}



void calCornerNormal(double coor_p[2], double coor_pl[2], double coor_pr[2],
						double normal_corner[2])
{
	double vec_plp[2], vec_ppr[2];
	vecSubVec2(coor_p, coor_pl, vec_plp);
	vecSubVec2(coor_pr, coor_p, vec_ppr);

	double length_plp = calNorm2(vec_plp);
	double length_ppr = calNorm2(vec_ppr);

	double normal_plp[2], normal_ppr[2];
	calRHNormal2(vec_plp, normal_plp);
	calRHNormal2(vec_ppr, normal_ppr);

	double half_directed_normal_plp[2], half_directed_normal_ppr[2];
	numMulVec2(length_plp/2.0, normal_plp, half_directed_normal_plp);
	numMulVec2(length_ppr/2.0, normal_ppr, half_directed_normal_ppr);

	double directed_length_corner[2];
	vecAddVec2(half_directed_normal_plp, half_directed_normal_ppr,
				directed_length_corner);

	double length_corner = calNorm2(directed_length_corner);

	numMulVec2(1.0/length_corner, directed_length_corner,
				normal_corner);
}
