#include <stdio.h>

#include "linear_algebra_2d.h"
#include "mesh_calculation.h"



// WARNING: This function does not contain coundary safty check.
size_t localToGlobalCellIdx(size_t idxP, size_t idxCP, 
							size_t NUM_NODES_X, size_t NUM_CELLS_X) 
{
	//printf("idxP, idxCP: %zu, %zu\n", idxP, idxCP);
	
	// locate the node P
	size_t row_node = idxP / NUM_NODES_X;
	size_t col_node = idxP % NUM_NODES_X;

	//printf("row_node, col_node: %zu, %zu\n", row_node, col_node);

	// left-down neigh cell
	size_t row_base = row_node - 1;
	size_t col_base = col_node - 1;

	//printf("row_base, col_base: %zu, %zu\n", row_base, col_base);

	size_t row_cell = row_base + idxCP / 2;
	size_t col_cell = col_base + idxCP % 2;

	//printf("row_cell, col_cell: %zu, %zu\n", row_cell, col_cell);

	size_t idxC = row_cell * NUM_CELLS_X + col_cell;

	return idxC;
}


size_t localToGlobalNodeIdx(size_t idxC, size_t idxPC, 
							size_t NUM_NODES_X, size_t NUM_CELLS_X) 
{
	printf("idxC, idxPC: %zu, %zu\n", idxC, idxPC);
	
	// locate the cell C
	size_t row_cell = idxC / NUM_CELLS_X;
	size_t col_cell = idxC % NUM_CELLS_X;

	//printf("row_cell, col_cell: %zu, %zu\n", row_cell, col_cell);

	// left-down neigh vertex
	size_t row_base = row_cell;
	size_t col_base = col_cell;

	//printf("row_base, col_base: %zu, %zu\n", row_base, col_base);

	size_t row_node = row_base + idxPC / 2;
	size_t col_node = col_base + idxPC % 2;

	//printf("row_node, col_node: %zu, %zu\n", row_node, col_node);

	size_t idxP = row_node * NUM_NODES_X + col_node;

	return idxP;
}


size_t globalToLocalNodeIdx(size_t idxC, size_t idxP, 
							size_t NUM_NODES_X, size_t NUM_CELLS_X) 
{
	//printf("idxC, idxPC: %zu, %zu\n", idxC, idxPC);
	
	// locate the cell C
	size_t row_cell = idxC / NUM_CELLS_X;
	size_t col_cell = idxC % NUM_CELLS_X;

	// locate the node P
	size_t row_node = idxP / NUM_NODES_X;
	size_t col_node = idxP % NUM_NODES_X;

	//printf("row_cell, col_cell: %zu, %zu\n", row_cell, col_cell);


	//printf("row_base, col_base: %zu, %zu\n", row_base, col_base);

	size_t row_vertex = row_node - row_cell;
	size_t col_vertex = col_node - col_cell;

	//printf("row_node, col_node: %zu, %zu\n", row_node, col_node);

	size_t idxPC = row_vertex * 2 + col_vertex;

	return idxPC;
}


size_t calNeighVertexIdx(size_t idxP, size_t idxC, int direction,
						 size_t NUM_NODES_X, size_t NUM_CELLS_X)
{
	size_t idxPC = globalToLocalNodeIdx(idxC, idxP, NUM_NODES_X, NUM_CELLS_X);

	size_t vertex_book[4] = {0, 1, 3, 2};

	int order_node = 0;
	for (int i = 0; i < 4; ++i) {
		if (vertex_book[i] == idxPC)
			order_node = i;
	}

	size_t idxPC_neigh = vertex_book[(order_node + direction + 4) % 4];

	return localToGlobalNodeIdx(idxC, idxPC_neigh, NUM_NODES_X, NUM_CELLS_X);
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
