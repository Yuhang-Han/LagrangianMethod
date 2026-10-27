#include <stdio.h>
#include "../mesh_calculation.h"

int main()
{
	const size_t NUM_NODES_X = 4;
	const size_t NUM_CELLS_X = 3;
	
/*
	for(size_t iP = 0; iP < 16; ++iP) {
		for(size_t iCP = 0; iCP < 4; ++iCP) {
			size_t idxC = localToGlobalCellIdx(iP, iCP, NUM_NODES_X, NUM_CELLS_X);
			printf("idxC: %zu\n", idxC);
		}
	}
*/

/*
	for(size_t iC = 0; iC < 9; ++iC) {
		for(int iPC = 0; iPC < 4; ++iPC) {
			size_t idxP = localToGlobalNodeIdx(iC, iPC, NUM_NODES_X, NUM_CELLS_X);
			printf("idxP: %zu\n", idxP);
			int idxPC = globalToLocalNodeIdx(iC, idxP, NUM_NODES_X, NUM_CELLS_X);
			printf("idxPC: %d\n", idxPC);

			if (idxPC != iPC) {
				printf("WRONG in globalToLocalNodeIdx()! \n\n");
			}
		}
	}
*/

	for(size_t iC = 0; iC < 9; ++iC) {
		for(size_t iPC = 0; iPC < 4; ++iPC) {
			size_t idxP = localToGlobalNodeIdx(iC, iPC, NUM_NODES_X, NUM_CELLS_X);
			printf("idxP: %zu\n", idxP);
			size_t idxPl = calNeighVertexIdx(idxP, iC, 1, 
											NUM_NODES_X, NUM_CELLS_X);
			size_t idxPr = calNeighVertexIdx(idxP, iC, -1, 
											NUM_NODES_X, NUM_CELLS_X);

			printf("Pl, Pr: %zu, %zu\n", idxPl, idxPr);

		}
	}


}
