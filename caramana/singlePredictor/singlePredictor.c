#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>

typedef struct {
    size_t (**vertices)[4][2];
    double  **mass;
	double  **pressure;
	double  **interior;
	
    size_t (*vertices_data)[4][2];
    double  *mass_data;
    double  *pressure_data;
    double  *interior_data;
} Grid;

typedef struct {
    size_t  (**center)[2];
    double  **mass;
    double (**velocity)[2];		// row first element's pointer
    double (**halftime_velocity)[2];
    
    size_t  (*center_data)[2];
    double  *mass_data;
    double (*velocity_data)[2]; // contiguous storage backing velocity
    double (*halftime_velocity_data)[2];
} DualGrid;

typedef struct {
	double (**vector)[2];
    double  **mass;
    double (**force)[2];
    
    double (*vector_data)[2];
    double  *mass_data;
    double (*force_data)[2];
} Corners;



double CalLength (double x0, double y0, double x1, double y1) {
	double length = sqrt( pow((x0-x1),2) + pow((y0-y1),2) );
	return length;
}


double  CalCellArea(Grid grid, double (**coor_node)[2], size_t idxCell0, size_t idxCell1) {
	double CrossProduct(double ax, double ay, double bx, double by, double cx, double cy) {
	    return (bx - ax) * (cy - ay)
	         - (by - ay) * (cx - ax);
	}

	int IsConvex(double ax, double ay, double bx, double by, double cx, double cy, double dx, double dy) {
	    double c1 = CrossProduct(ax, ay, bx, by, cx, cy);
	    double c2 = CrossProduct(bx, by, cx, cy, dx, dy);
	    double c3 = CrossProduct(cx, cy, dx, dy, ax, ay);
	    double c4 = CrossProduct(dx, dy, ax, ay, bx, by);

	    return (c1 > 0 && c2 > 0 && c3 > 0 && c4 > 0) ||
	           (c1 < 0 && c2 < 0 && c3 < 0 && c4 < 0);
	}
	
	size_t v0[2] = {grid.vertices[idxCell0][idxCell1][0][0], grid.vertices[idxCell0][idxCell1][0][1]};
	size_t v1[2] = {grid.vertices[idxCell0][idxCell1][1][0], grid.vertices[idxCell0][idxCell1][1][1]};
	size_t v2[2] = {grid.vertices[idxCell0][idxCell1][2][0], grid.vertices[idxCell0][idxCell1][2][1]};
	size_t v3[2] = {grid.vertices[idxCell0][idxCell1][3][0], grid.vertices[idxCell0][idxCell1][3][1]};

	double ax = coor_node[v0[0]][v0[1]][0];
	double ay = coor_node[v0[0]][v0[1]][1];
	double bx = coor_node[v1[0]][v1[1]][0];
	double by = coor_node[v1[0]][v1[1]][1];
	double cx = coor_node[v2[0]][v2[1]][0];
	double cy = coor_node[v2[0]][v2[1]][1];
	double dx = coor_node[v3[0]][v3[1]][0];
	double dy = coor_node[v3[0]][v3[1]][1];

	if ( IsConvex(ax, ay, bx, by, cx, cy, dx, dy) ) {
		double twiceArea = 
		          ax * by + bx * cy
		        + cx * dy + dx * ay
		        - ay * bx - by * cx
		        - cy * dx - dy * ax; 
		return 0.5 * twiceArea;
	}
	else {
		printf("ERROR: Get non-convex quad in CallCellArea().\n");
		return 0.0;
	}
}

double CalPressure(double density, double interior, double gamma) {
	return (gamma - 1) * density * interior;
}



double InitializeCellAveDensity(size_t idxCell0, size_t idxCell1) {
	
}

double InitializeDualAveMomentum(size_t idxCell0, size_t idxCell1) {
	
}

double InitializeCellPressure(size_t idxCell0, size_t idxCell1) {

}

double InitializeCellInterior(size_t idxCell0, size_t idxCell1) {
	
}


void MapNeighborDual(size_t iDual0, size_t iDual1, size_t iNeigh, size_t idxNeigh[2]) {
// map {0,1,2,3} to {(-1,0), (0,-1), (1,0), (0,1)}
	idxNeigh[0] = iDual0 + (iNeigh == 2) - (iNeigh == 0);
	idxNeigh[1] = iDual1 + (iNeigh == 3) - (iNeigh == 1);
}

double CalDistance(size_t node0[2], size_t node1[2], double (**coor_node)[2]) {
// general distance calculator of two 2d arrays
	double diffX = coor_node[node0[0]] - coor_node[node1[0]];
	double diffY = coor_node[node0[1]] - coor_node[node1[1]];
	return sqrt( pow(diffX,2) + pow(diffY,2) );
}

double GetNorm2(double vector[2]) {
// 2d vector 2-norm
	return sqrt( pow(vector[0],2) + pow(vector[1],2) );
}

double GetMinDouble(double a, double b) {
	if(a < b)
		return a;
	else 
		return b;
}

double GetMaxDouble(double a, double b) {
	if(a > b)
		return a;
	else 
		return b;
}

double CalNeighAcousticSpeed(size_t iDual0, size_t iDual1, Grid grid, double gamma) {
// use max of neighbor cell acou speed as dual's acou speed
	// find dual's neighbor cells
	double dualAcoustic = 0.0;
	for(size_t iNeigh = 0; iNeigh < 4; ++iNeigh) {
		size_t idxNeighCell0 = iDual0 - iNeigh/2;
		size_t idxNeighCell1 = iDual1 - iNeigh%2;

		double cellInterior = grid.interior[idxNeighCell0][idxNeighCell1];

		double cellAcoustic = sqrt(gamma * (gamma-1) * cellInterior);

		dualAcoustic = GetMaxDouble(cellAcoustic, dualAcoustic);
	}

	return dualAcoustic;
}

void ApplyBoundatyCondition(Grid grid, char *bc, size_t numGhCell0, size_t numGhCell1) {

	if(strcmp(bc, "free") == 0) {
		// (0,1) -- (0,N-1) row
		for(size_t iBoundCell = 1; iBoundCell < numGhCell1-1; ++iBoundCell) {

			grid.mass[iBoundCell][0] 		= grid.mass[iBoundCell][1];
			grid.pressure[iBoundCell][0] 	= grid.pressure[iBoundCell][1];
			grid.interior[iBoundCell][0] 	= grid.interior[iBoundCell][1];
		}
		// (N,1) -- (N,N-1) row
		for(size_t iBoundCell = 1; iBoundCell < numGhCell1-1; ++iBoundCell) {

			grid.mass[iBoundCell][numGhCell0-1] 		= grid.mass[iBoundCell][numGhCell0-2];
			grid.pressure[iBoundCell][numGhCell0-1] 	= grid.pressure[iBoundCell][numGhCell0-2];
			grid.interior[iBoundCell][numGhCell0-1] 	= grid.interior[iBoundCell][numGhCell0-2];
		}
		// (1,0) -- (N-1,0) column
		for(size_t iBoundCell = 1; iBoundCell < numGhCell0-1; ++iBoundCell) {

			grid.mass[0][iBoundCell] 		= grid.mass[1][iBoundCell];
			grid.pressure[0][iBoundCell]	= grid.pressure[1][iBoundCell];
			grid.interior[0][iBoundCell] 	= grid.interior[1][iBoundCell];
		}
		// (N,1) -- (N,N-1) column
		for(size_t iBoundCell = 1; iBoundCell < numGhCell0-1; ++iBoundCell) {

			grid.mass[iBoundCell][numGhCell0-1] 		= grid.mass[iBoundCell][numGhCell0-2];
			grid.pressure[iBoundCell][numGhCell0-1] 	= grid.pressure[iBoundCell][numGhCell0-2];
			grid.interior[iBoundCell][numGhCell0-1] 	= grid.interior[iBoundCell][numGhCell0-2];
		}

		// four cells at corners
		size_t cornerCells[4][2] = { {0,0}, {0,numGhCell0-1}, {numGhCell1-1,0}, {numGhCell1-1, numGhCell0-1} };
		size_t cornerNeigh[4][4] = { {1,0, 0,1}, 
									 {1,numGhCell0-1, 0,numGhCell0-2},
									 {numGhCell1-1,1, numGhCell1-2,0},
									 {numGhCell1-2, numGhCell0-1, numGhCell1-1, numGhCell0-2} };
		for(size_t iCornerCell = 0; iCornerCell < 4; ++iCornerCell) {

			grid.mass[ cornerCells[iCornerCell][0] ][ cornerCells[iCornerCell][1] ] = \
				0.5 * ( grid.mass[ cornerNeigh[iCornerCell][0] ][ cornerNeigh[iCornerCell][1] ] 
					   +grid.mass[ cornerNeigh[iCornerCell][2] ][ cornerNeigh[iCornerCell][3] ] );
					   
			grid.pressure[ cornerCells[iCornerCell][0] ][ cornerCells[iCornerCell][1] ] = \
				0.5 * ( grid.pressure[ cornerNeigh[iCornerCell][0] ][ cornerNeigh[iCornerCell][1] ] 
					   +grid.pressure[ cornerNeigh[iCornerCell][2] ][ cornerNeigh[iCornerCell][3] ] );
					   
			grid.interior[ cornerCells[iCornerCell][0] ][ cornerCells[iCornerCell][1] ] = \
				0.5 * ( grid.interior[ cornerNeigh[iCornerCell][0] ][ cornerNeigh[iCornerCell][1] ] 
					   +grid.interior[ cornerNeigh[iCornerCell][2] ][ cornerNeigh[iCornerCell][3] ] );
		}

		
		
	}
	else if (strcmp(bc, "periodic") == 0) {
		
	}
	else {
		printf("Wrong boundary condition in ApplyBoundatyCondition()\n");
	}
}



int main()
{

	const double GAMMA 	= 1.4;
	const double CFL	= 0.2;

	const double X0 = 1.0;
	const double X1 = 1.0;

	const double T  = 1.0;

	const size_t numCell0 = 10;
	const size_t numCell1 = 10;

	
	const size_t numGhCell0 = numCell0 + 2;
	const size_t numGhCell1 = numCell1 + 2;

	const size_t numDual0 = numCell0 + 1;
	const size_t numDual1 = numCell1 + 1;
	const size_t numGhDual0 = numDual0 + 2;
	const size_t numGhDual1 = numDual1 + 2;

	const size_t numCorner0 = 2 * numCell0;
	const size_t numCorner1 = 2 * numCell1;
	const size_t numGhCorner0 = 2 * numGhCell0;
	const size_t numGhCorner1 = 2 * numGhCell1;

	const size_t numNode0 = numCell0 + 1;
	const size_t numNode1 = numCell1 + 1;
	const size_t numGhNode0 = numGhCell0 + 1;
	const size_t numGhNode1 = numGhCell1 + 1;

	const size_t numCenter0 = numCell0;
	const size_t numCenter1 = numCell1;
	const size_t numGhCenter0 = numGhCell0;
	const size_t numGhCenter1 = numGhCell1;


/***********************
* Memory Allocate
************************/

	// Grid
	// dimension: [numGhCell0][numGhCell1]
	Grid grid;
	grid.vertices = malloc(numGhCell0 * sizeof(*grid.vertices));
	grid.mass     = malloc(numGhCell0 * sizeof(*grid.mass));
	grid.pressure = malloc(numGhCell0 * sizeof(*grid.pressure));
	grid.interior = malloc(numGhCell0 * sizeof(*grid.interior));

	grid.vertices_data 	= calloc(numGhCell0 * numGhCell1, sizeof(*grid.vertices_data));
	grid.mass_data 		= calloc(numGhCell0 * numGhCell1, sizeof(*grid.mass_data));
	grid.pressure_data 	= calloc(numGhCell0 * numGhCell1, sizeof(*grid.pressure_data));
	grid.interior_data 	= calloc(numGhCell0 * numGhCell1, sizeof(*grid.interior_data));

	for (size_t i = 0; i < numGhCell0; i++) {
	    grid.vertices[i] = grid.vertices_data + i * numGhCell1;
	    grid.mass[i]     = grid.mass_data     + i * numGhCell1;
	    grid.pressure[i] = grid.pressure_data + i * numGhCell1;
	    grid.interior[i] = grid.interior_data + i * numGhCell1;
	}


	// DualGrid
	// dimension: [numGhDual0][numGhDual1]
	DualGrid dualgrid;
	dualgrid.center            = malloc(numGhDual0 * sizeof(*dualgrid.center));
	dualgrid.mass              = malloc(numGhDual0 * sizeof(*dualgrid.mass));
	dualgrid.velocity          = malloc(numGhDual0 * sizeof(*dualgrid.velocity));
	dualgrid.halftime_velocity = malloc(numGhDual0 * sizeof(*dualgrid.halftime_velocity));

	dualgrid.center_data 			= calloc(numGhDual0 * numGhDual1, sizeof(*dualgrid.center_data));
	dualgrid.mass_data 				= calloc(numGhDual0 * numGhDual1, sizeof(*dualgrid.mass_data));
	dualgrid.velocity_data 			= calloc(numGhDual0 * numGhDual1, sizeof(*dualgrid.velocity_data));
	dualgrid.halftime_velocity_data = calloc(numGhDual0 * numGhDual1, sizeof(*dualgrid.halftime_velocity_data));

	for (size_t i = 0; i < numGhDual0; i++) {
	    dualgrid.center[i] 				= dualgrid.center_data 				+ i * numGhDual1;
	    dualgrid.mass[i] 				= dualgrid.mass_data 				+ i * numGhDual1;
	    dualgrid.velocity[i] 			= dualgrid.velocity_data 			+ i * numGhDual1;
	    dualgrid.halftime_velocity[i] 	= dualgrid.halftime_velocity_data 	+ i * numGhDual1;
	}


	// Corners
	// dimension: [numGhCorner0][numGhCorner1]
	Corners corners;
	corners.vector = malloc(numGhCorner0 * sizeof(*corners.vector));
	corners.mass   = malloc(numGhCorner0 * sizeof(*corners.mass));
	corners.force  = malloc(numGhCorner0 * sizeof(*corners.force));

	corners.vector_data = calloc(numGhCorner0 * numGhCorner1, sizeof(*corners.vector_data));
	corners.mass_data 	= calloc(numGhCorner0 * numGhCorner1, sizeof(*corners.mass_data));
	corners.force_data 	= calloc(numGhCorner0 * numGhCorner1, sizeof(*corners.force_data));

	for (size_t i = 0; i < numGhCorner0; i++) {
	    corners.vector[i] = corners.vector_data + i * numGhCorner1;
	    corners.mass[i]   = corners.mass_data   + i * numGhCorner1;
	    corners.force[i]  = corners.force_data  + i * numGhCorner1;
	}


	// node coor
	// dimension: [numGhNode0][numGhNode1]
	double (**coor_node)  [2] = malloc(numGhNode0 	* sizeof(*coor_node));
	double (**coor_center)[2] = malloc(numGhCenter0 * sizeof(*coor_center));
	
	double (*coor_node_data)  [2] 	= calloc(numGhNode0   *	numGhNode1, 	sizeof(*coor_node_data));
	double (*coor_center_data)[2] 	= calloc(numGhCenter0 * numGhCenter1, 	sizeof(*coor_center_data));
	
	for (size_t i = 0; i < numGhNode0; i++) {
		coor_node[i] 	= coor_node_data 	+ i * numGhNode1;
		coor_center[i] 	= coor_center_data 	+ i * numGhCenter1;
	}


/*****************************
* Initialize
******************************/
//TODO: 检查所有索引应该用Gh 还是非Gh



// (uniform) Grid geometry

	// node and center coordinates
	double cellH0 = X0 / numGhCell0;
	double cellH1 = X1 / numGhCell1;
	
	for (size_t iGhNode0 = 0; iGhNode0 < numGhNode0; ++iGhNode0) {
		for (size_t iGhNode1 = 0; iGhNode1 < numGhNode1; ++iGhNode1) {
			coor_node[iGhNode0][iGhNode1][0] = iGhNode0 * cellH0;
			coor_node[iGhNode0][iGhNode1][1] = iGhNode1 * cellH1;
		}
	}

	for (size_t iGhCenter0 = 0; iGhCenter0 < numGhCenter0; ++iGhCenter0) {
		for (size_t iGhCenter1 = 0; iGhCenter1 < numGhCenter1; ++iGhCenter1) {
			coor_center[iGhCenter0][iGhCenter1][0] = (iGhCenter0 + 0.5) * cellH0;
			coor_center[iGhCenter0][iGhCenter1][1] = (iGhCenter1 + 0.5) * cellH1;
		}
	}

	// cell vertices
	for (size_t iGhCell0 = 0; iGhCell0 < numGhCell0; ++iGhCell0) {
		for (size_t iGhCell1 = 0; iGhCell1 < numGhCell1; ++iGhCell1) {
		
			for (size_t iCellVertex = 0; iCellVertex < 4; ++iCellVertex) {
				grid.vertices[iGhCell0][iGhCell1][iCellVertex][0] = iGhCell0 + iCellVertex / 2;
				grid.vertices[iGhCell0][iGhCell1][iCellVertex][1] = iGhCell1 + iCellVertex % 2;
				/*  2   	3
					*-------*
					| 	| 	|
					--------- 
					|   | 	|
					*-------*
					0    	1	
				*/
			}
		}
	}

	// dual centers
	for (size_t iGhDual0 = 0; iGhDual0 < numGhDual0; ++iGhDual0) {
		for (size_t iGhDual1 = 0; iGhDual1 < numGhDual1; ++iGhDual1) {
		
			dualgrid.center[iGhDual0][iGhDual1][0] = iGhDual0;
			dualgrid.center[iGhDual0][iGhDual1][1] = iGhDual1;
		}
	}


// Initial condition

	// corner mass
	for (size_t iGhCorner0 = 0; iGhCorner0 < numGhCorner0; ++iGhCorner0) {
		for (size_t iGhCorner1 = 0; iGhCorner1 < numGhCorner1; ++iGhCorner1) {
			
			size_t idxFatherCell0 = iGhCorner0 / 2;
			size_t idxFatherCell1 = iGhCorner1 / 2;

			double cornerArea = 1.0/4 * cellH0 * cellH1;
			
			corners.mass[iGhCorner0][iGhCorner0] = cornerArea * InitializeCellAveDensity(idxFatherCell0, idxFatherCell1);
		}
	}

	// cell
	for (size_t iGhCell0 = 0; iGhCell0 < numGhCell0; ++iGhCell0) {
		for (size_t iGhCell1 = 0; iGhCell1 < numGhCell1; ++iGhCell1) {

			// cell mass
			grid.mass[iGhCell0][iGhCell1] = 0.0;
			for (size_t iSubCorner = 0; iSubCorner < 4; ++iSubCorner) {
				size_t idxCorner0 = 2*iGhCell0 + iSubCorner / 2;
				size_t idxCorner1 = 2*iGhCell1 + iSubCorner % 2;
				grid.mass[iGhCell0][iGhCell1] += corners.mass[idxCorner0][idxCorner1];
			}

			// cell pressure
			grid.pressure[iGhCell0][iGhCell1] = InitializeCellPressure(iGhCell0, iGhCell1);

			// cell interior energy
			grid.interior[iGhCell0][iGhCell1] = InitializeCellInterior(iGhCell0, iGhCell1);
		}
	}

	// dual
	for (size_t iGhDual0 = 0; iGhDual0 < numGhDual0; ++iGhDual0) {
		for (size_t iGhDual1 = 0; iGhDual1 < numGhDual1; ++iGhDual1){

			// dual mass
			dualgrid.mass[iGhDual0][iGhDual1] = 0.0;
			for (size_t iSubCorner = 0; iSubCorner < 4; ++iSubCorner) {
				size_t idxCorner0 = 2*iGhDual0 - 1 + iSubCorner / 2;
				size_t idxCorner1 = 2*iGhDual1 - 1 + iSubCorner % 2;
				dualgrid.mass[iGhDual0][iGhDual1] += corners.mass[idxCorner0][idxCorner1];
			}

		}
	}




/****************************
* Time Advance
*****************************/

	double t = 0.0;

	while (t < T) {
		// Apply boundary condition
		ApplyBoundatyCondition(grid, "free", numGhCell0, numGhCell1);	// "free" or "periodic"

	
// CFL 条件这里有待考虑, 是否能更细致

		// Determine time step length
		double tau = T;
		
		for (size_t iDual0 = 0; iDual0 < numDual0; ++iDual0) {
			for (size_t iDual1 = 0; iDual1 < numDual1; ++iDual1){

				double localMinDistance = X0;
				for(size_t iNeigh = 0; iNeigh < 4; ++iNeigh) {
				// compare every physical dual with its neighbors
					size_t idxNeigh[2];
					MapNeighborDual(iDual0, iDual1, iNeigh, idxNeigh);
					// map local idx to global idx
// 这里距离用了对偶单元中心距离, 不知是否合理
					double localDistance = CalDistance(dualgrid.center[iDual0][iDual1], 
													   dualgrid.center[idxNeigh[0]][idxNeigh[1]], coor_node);
					localMinDistance = GetMinDouble(localDistance, localMinDistance);
				}
			
				double localAbsSpeed = GetNorm2(dualgrid.velocity[iDual0][iDual1]);
				double localAcousticSpeed = CalNeighAcousticSpeed(iDual0, iDual1, grid, GAMMA);
// 这里也有一些小问题, 声速和速度的定义不在同一个地方

				double localMaxSpeed = localAbsSpeed + localAcousticSpeed;
				double localTau = localMinDistance / localMaxSpeed / CFL;

				tau = GetMinDouble(localTau, tau);
			}
		}

		t += tau;
		
		

	
		// Calculate corner force
		for (size_t iGhCorner0 = 0; iGhCorner0 < numGhCorner0; ++iGhCorner0) {
			for (size_t iGhCorner1 = 0; iGhCorner1 < numGhCorner1; ++iGhCorner1){
			
				// calculate corner vector
				size_t idxFatherCell0 = iGhCorner0 / 2;
				size_t idxFatherCell1 = iGhCorner1 / 2;

				double lenCellEdge[4];
				size_t refNode0 = grid.vertices[idxFatherCell0][idxFatherCell1][0][0];
				size_t refNode1 = grid.vertices[idxFatherCell0][idxFatherCell1][0][1];
				for(size_t iCellEdge = 0; iCellEdge < 4; ++iCellEdge) {
					size_t tmp0 = iCellEdge ^ (iCellEdge >> 1);
					size_t tmp1 = (iCellEdge+1) ^ ((iCellEdge+1) >> 1);
					lenCellEdge[iCellEdge] = CalLength(refNode0 + ((tmp0 >> 1) & 1), refNode1 + (tmp0 & 1),
													   refNode0 + ((tmp1 >> 1) & 1), refNode1 + (tmp1 & 1));
				}
					/*     2
						*-------*
						| 2	| 3	|
					3	--------- 1
						| 0 | 1	|
						*-------*
						    0		
					*/
				double cornerVector[2];
				size_t subCornerId = iGhCorner0%2 + 2 * iGhCorner1%2;
				double norm;
				switch(subCornerId){
					case 0: 
						norm = sqrt(pow(lenCellEdge[3],2) + pow(lenCellEdge[0],2));
						cornerVector[0] = -lenCellEdge[3] / norm;
						cornerVector[1] = -lenCellEdge[0] / norm;
						break;
					case 1:
						norm = sqrt(pow(lenCellEdge[0],2) + pow(lenCellEdge[1],2));
						cornerVector[0] =  lenCellEdge[1] / norm;
						cornerVector[1] = -lenCellEdge[0] / norm;	
						break;
					case 2:
						norm = sqrt(pow(lenCellEdge[2],2) + pow(lenCellEdge[3],2));
						cornerVector[0] = -lenCellEdge[3] / norm;
						cornerVector[1] =  lenCellEdge[2] / norm;	
						break;
					case 3:
						norm = sqrt(pow(lenCellEdge[1],2) + pow(lenCellEdge[2],2));
						cornerVector[0] =  lenCellEdge[1] / norm;
						cornerVector[1] =  lenCellEdge[2] / norm;	
						break;
				}

				// corner force
				double cellPressure = grid.pressure[idxFatherCell0][idxFatherCell1];
				corners.force[iGhCorner0][iGhCorner1][0] = cellPressure * cornerVector[0];
				corners.force[iGhCorner0][iGhCorner1][1] = cellPressure * cornerVector[1];
				
			}
		}
		

		// Update physical dual momentum 
		for (size_t iDual0 = 0; iDual0 < numGhDual0; ++iDual0) {
			for (size_t iDual1 = 0; iDual1 < numGhDual1; ++iDual1) {

				double acceleration[2] = {0.0, 0.0};
				for (size_t iSubCorner = 0; iSubCorner < 4; ++iSubCorner) {
					size_t idxCorner0 = 2*iDual0 - 1 + iSubCorner / 2;
					size_t idxCorner1 = 2*iDual1 - 1 + iSubCorner % 2;
					
					acceleration[0] += corners.force[idxCorner0][idxCorner1][0] / dualgrid.mass[iDual0][iDual1];
					acceleration[1] += corners.force[idxCorner0][idxCorner1][1] / dualgrid.mass[iDual0][iDual1];
				}

				dualgrid.halftime_velocity[iDual0][iDual1][0] = \
						 dualgrid.velocity[iDual0][iDual1][0] + tau/2.0 * acceleration[0];
				dualgrid.halftime_velocity[iDual0][iDual1][1] = \
						 dualgrid.velocity[iDual0][iDual1][1] + tau/2.0 * acceleration[1];

				dualgrid.velocity[iDual0][iDual1][0] += tau * acceleration[0];
				dualgrid.velocity[iDual0][iDual1][1] += tau * acceleration[1];

			}
		}

		// Update physical cell interior energy
		for (size_t iCell0 = 0; iCell0 < numCell0; ++iCell0) {
			for (size_t iCell1 = 0; iCell1 < numCell1; ++iCell1) {

				for (size_t iSubCorner = 0; iSubCorner < 4; ++iSubCorner) {

					double corner_velocity[2];
					size_t idxFatherDual0 = iCell0 + iSubCorner / 2;
					size_t idxFatherDual1 = iCell1 + iSubCorner % 2;
					corner_velocity[0] = dualgrid.halftime_velocity[idxFatherDual0][idxFatherDual1][0];
					corner_velocity[1] = dualgrid.halftime_velocity[idxFatherDual0][idxFatherDual1][1];


					size_t idxCorner0 = 2*iCell0 + iSubCorner / 2;
					size_t idxCorner1 = 2*iCell1 + iSubCorner % 2;
					
					double cornerWork = \
							corners.force[idxCorner0][idxCorner1][0] * corner_velocity[0] \
						  + corners.force[idxCorner0][idxCorner1][1] * corner_velocity[1];

					grid.interior[iCell0][iCell1] += - tau * cornerWork / grid.mass[iCell0][iCell1];
			
				}
			}
		}

		// Update dual center (cell vertex) position
		for (size_t iDual0 = 0; iDual0 < numDual0; ++iDual0) {
			for (size_t iDual1 = 0; iDual1 < numDual1; ++iDual1) {
				size_t idxNode0 = dualgrid.center[iDual0][iDual1][0];
				size_t idxNode1 = dualgrid.center[iDual0][iDual1][1];
				coor_node[idxNode0][idxNode1][0] += tau * dualgrid.halftime_velocity[iDual0][iDual1][0];
				coor_node[idxNode0][idxNode1][1] += tau * dualgrid.halftime_velocity[iDual0][iDual1][1];
			}
		}

		// Update cell center position
// 似乎不会直接用到单元中心?

		

		// Update cell pressure
		for (size_t iCell0 = 0; iCell0 < numCell0; ++iCell0) {
			for (size_t iCell1 = 0; iCell1 < numCell1; ++iCell1) {
			
				double cellArea = CalCellArea(grid, coor_node, iCell0, iCell1);
				double cellDensity  = grid.mass[iCell0][iCell1] / cellArea;

				grid.pressure[iCell0][iCell1] = CalPressure(cellDensity, grid.interior[iCell0][iCell1], GAMMA);
			}
		}

	}
	







    



/**************************
* Free Memory
**************************/

	free(grid.vertices_data);
	free(grid.mass_data);
	free(grid.pressure_data);
	free(grid.interior_data);

	free(grid.vertices);
	free(grid.mass);
	free(grid.pressure);
	free(grid.interior);


	free(dualgrid.center_data);
	free(dualgrid.mass_data);
	free(dualgrid.velocity_data);
	free(dualgrid.halftime_velocity_data);

	free(dualgrid.center);
	free(dualgrid.mass);
	free(dualgrid.velocity);
	free(dualgrid.halftime_velocity);


	free(corners.vector_data);
	free(corners.mass_data);
	free(corners.force_data);

	free(corners.vector);
	free(corners.mass);
	free(corners.force);
	
}
