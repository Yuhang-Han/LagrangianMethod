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

static const int cornerOffset[4][2] = {
    {0, 0},   // 0: bottom-left
    {1, 0},   // 1: bottom-right
    {1, 1},   // 2: top-right
    {0, 1}    // 3: top-left
};



double CalLength (size_t x0, size_t y0, size_t x1, size_t y1, double (**coor_node)[2]) {
	double coor_x0 = coor_node[x0][y0][0];
	double coor_y0 = coor_node[x0][y0][1];
	double coor_x1 = coor_node[x1][y1][0];
	double coor_y1 = coor_node[x1][y1][1];
	double length = sqrt( pow((coor_x0-coor_x1),2) + pow((coor_y0-coor_y1),2) );
	return length;
}

void CalCellCornerVectors(
    Grid grid,
    double (**coor_node)[2],
    size_t ic,
    size_t jc,
    double C[4][2])
{
    double X[4][2];
    double N[4][2];

    /* Obtain the four physical vertex coordinates.
       Vertex order: BL, BR, TR, TL = CCW. */
    for (size_t k = 0; k < 4; ++k) {
        size_t ni = grid.vertices[ic][jc][k][0];
        size_t nj = grid.vertices[ic][jc][k][1];

        X[k][0] = coor_node[ni][nj][0];
        X[k][1] = coor_node[ni][nj][1];
    }

    /* Full outward edge normals.
       For CCW polygon:
           edge = (dx,dy)
           outward normal = (dy,-dx)
    */
    for (size_t k = 0; k < 4; ++k) {
        size_t kp = (k + 1) % 4;

        double dx = X[kp][0] - X[k][0];
        double dy = X[kp][1] - X[k][1];

        N[k][0] =  dy;
        N[k][1] = -dx;
    }

    /* Corner vector = sum of the two half-edge normals */
    for (size_t k = 0; k < 4; ++k) {
        size_t km = (k + 3) % 4;

        C[k][0] = 0.5 * (N[km][0] + N[k][0]);
        C[k][1] = 0.5 * (N[km][1] + N[k][1]);
    }
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
	return 1.0;
}

double InitializeDualAveMomentum(size_t idxCell0, size_t idxCell1) {
	return 0.0;
}

double InitializeCellPressure(size_t idxCell0, size_t idxCell1) {
	return 1.0;
}

double InitializeCellInterior(size_t idxCell0, size_t idxCell1) {
	return 2.5;
}


void MapNeighborDual(size_t iDual0, size_t iDual1, size_t iNeigh, size_t idxNeigh[2]) {
// map {0,1,2,3} to {(-1,0), (0,-1), (1,0), (0,1)}
	idxNeigh[0] = iDual0 + (iNeigh == 2) - (iNeigh == 0);
	idxNeigh[1] = iDual1 + (iNeigh == 3) - (iNeigh == 1);
}

double CalDistance(size_t node0[2], size_t node1[2], double (**coor_node)[2]) {
// general distance calculator of two 2d arrays
	double diffX = coor_node[node0[0]][node0[1]][0]
	      		 - coor_node[node1[0]][node1[1]][0];
	double diffY = coor_node[node0[0]][node0[1]][1]
	      		 - coor_node[node1[0]][node1[1]][1];
	
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

double CalTotalEnergy(
    Grid grid,
    DualGrid dualgrid,
    size_t numCell0,
    size_t numCell1,
    size_t numDual0,
    size_t numDual1)
{
    double Eint = 0.0;
    double Ekin = 0.0;

    for (size_t i = 1; i <= numCell0; ++i) {
        for (size_t j = 1; j <= numCell1; ++j) {

            Eint += grid.mass[i][j]
                  * grid.interior[i][j];
        }
    }

    for (size_t i = 1; i <= numDual0; ++i) {
        for (size_t j = 1; j <= numDual1; ++j) {

            double vx = dualgrid.velocity[i][j][0];
            double vy = dualgrid.velocity[i][j][1];

            Ekin += 0.5 * dualgrid.mass[i][j]
                         * (vx*vx + vy*vy);
        }
    }

    return Eint + Ekin;
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
//	double (**coor_center)[2] = malloc(numGhCenter0 * sizeof(*coor_center));
	
	double (*coor_node_data)  [2] 	= calloc(numGhNode0   *	numGhNode1, 	sizeof(*coor_node_data));
//	double (*coor_center_data)[2] 	= calloc(numGhCenter0 * numGhCenter1, 	sizeof(*coor_center_data));
	
	for (size_t i = 0; i < numGhNode0; i++) {
		coor_node[i] 	= coor_node_data 	+ i * numGhNode1;
		
	}
/*
	for (size_t i = 0; i < numGhCenter0; i++) {
		coor_center[i] 	= coor_center_data 	+ i * numGhCenter1;
	}
*/
	


/*****************************
* Initialize
******************************/
//TODO: 检查所有索引应该用Gh 还是非Gh



// (uniform) Grid geometry

	// node and center coordinates
	double cellH0 = X0 / numCell0;
	double cellH1 = X1 / numCell1;
	
	for (size_t iGhNode0 = 0; iGhNode0 < numGhNode0; ++iGhNode0) {
		for (size_t iGhNode1 = 0; iGhNode1 < numGhNode1; ++iGhNode1) {
			coor_node[iGhNode0][iGhNode1][0] = ((double)iGhNode0 -1.0) * cellH0;
			coor_node[iGhNode0][iGhNode1][1] = ((double)iGhNode1 -1.0) * cellH1;
		}
	}
/*
	for (size_t iGhCenter0 = 0; iGhCenter0 < numGhCenter0; ++iGhCenter0) {
		for (size_t iGhCenter1 = 0; iGhCenter1 < numGhCenter1; ++iGhCenter1) {
			coor_center[iGhCenter0][iGhCenter1][0] = (iGhCenter0 + 0.5) * cellH0;
			coor_center[iGhCenter0][iGhCenter1][1] = (iGhCenter1 + 0.5) * cellH1;
		}
	}
*/
	// cell vertices
	for (size_t iGhCell0 = 0; iGhCell0 < numGhCell0; ++iGhCell0) {
		for (size_t iGhCell1 = 0; iGhCell1 < numGhCell1; ++iGhCell1) {
		
			for (size_t iCellVertex = 0; iCellVertex < 4; ++iCellVertex) {
				grid.vertices[iGhCell0][iGhCell1][iCellVertex][0] = iGhCell0 + cornerOffset[iCellVertex][0];
				grid.vertices[iGhCell0][iGhCell1][iCellVertex][1] = iGhCell1 + cornerOffset[iCellVertex][1];
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
			
			corners.mass[iGhCorner0][iGhCorner1] = cornerArea * InitializeCellAveDensity(idxFatherCell0, idxFatherCell1);
		}
	}

	// cell
	for (size_t iGhCell0 = 0; iGhCell0 < numGhCell0; ++iGhCell0) {
		for (size_t iGhCell1 = 0; iGhCell1 < numGhCell1; ++iGhCell1) {

			// cell mass
			grid.mass[iGhCell0][iGhCell1] = 0.0;
			for (size_t iSubCorner = 0; iSubCorner < 4; ++iSubCorner) {
				size_t idxCorner0 = 2*iGhCell0 + cornerOffset[iSubCorner][0];
				size_t idxCorner1 = 2*iGhCell1 + cornerOffset[iSubCorner][1];
				grid.mass[iGhCell0][iGhCell1] += corners.mass[idxCorner0][idxCorner1];
			}

			// cell pressure
			grid.pressure[iGhCell0][iGhCell1] = InitializeCellPressure(iGhCell0, iGhCell1);

			// cell interior energy
			grid.interior[iGhCell0][iGhCell1] = InitializeCellInterior(iGhCell0, iGhCell1);
		}
	}

	// dual
	for (size_t iDual0 = 1; iDual0 <= numDual0; ++iDual0) {
		for (size_t iDual1 = 1; iDual1 <= numDual1; ++iDual1){

			// dual mass
			dualgrid.mass[iDual0][iDual1] = 0.0;
			for (size_t iSubCorner = 0; iSubCorner < 4; ++iSubCorner) {
				size_t idxCorner0 = 2*iDual0 - 1 + cornerOffset[iSubCorner][0];
				size_t idxCorner1 = 2*iDual1 - 1 + cornerOffset[iSubCorner][1];
				dualgrid.mass[iDual0][iDual1] += corners.mass[idxCorner0][idxCorner1];
			}

		}
	}




/****************************
* Time Advance
*****************************/

	double t = 0.0;
	size_t tcounter = 0;
	
	double E0 = CalTotalEnergy(grid, dualgrid, numCell0, numCell1, numDual0, numDual1);
	printf(
	    "t=%e E=%20.16e \n",
	    t, E0);

	while (t < T) {
		// Apply boundary condition
		ApplyBoundatyCondition(grid, "free", numGhCell0, numGhCell1);	// "free" or "periodic"
// 边界处理需要细考虑


		// Determine time step length
/*
		double tau = T;
		
		for (size_t iDual0 = 1; iDual0 < numDual0; ++iDual0) {
			for (size_t iDual1 = 1; iDual1 < numDual1; ++iDual1){

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
				double localTau = CFL * localMinDistance / localMaxSpeed;

				tau = GetMinDouble(localTau, tau);
			}
		}
*/
		double tau = 1E-4;
// 拉氏的CFL 条件与Euler 不同


		t += tau;
		++tcounter;
		
		

	
		// Calculate corner force
		for (size_t ic = 0; ic < numGhCell0; ++ic) {
		    for (size_t jc = 0; jc < numGhCell1; ++jc) {

		        double C[4][2];
		        CalCellCornerVectors(grid, coor_node, ic, jc, C);

		        for (size_t k = 0; k < 4; ++k) {

		            size_t I = 2*ic + cornerOffset[k][0];
		            size_t J = 2*jc + cornerOffset[k][1];

		            corners.vector[I][J][0] = C[k][0];
		            corners.vector[I][J][1] = C[k][1];

		            corners.force[I][J][0]
		                = grid.pressure[ic][jc] * C[k][0];

		            corners.force[I][J][1]
		                = grid.pressure[ic][jc] * C[k][1];
		        }
		    }
		}
		

		// Update physical dual momentum 
		for (size_t iDual0 = 1; iDual0 <= numDual0; ++iDual0) {
			for (size_t iDual1 = 1; iDual1 <= numDual1; ++iDual1) {

				double acceleration[2] = {0.0, 0.0};
				for (size_t iSubCorner = 0; iSubCorner < 4; ++iSubCorner) {
					size_t idxCorner0 = 2*iDual0 - 1 + cornerOffset[iSubCorner][0];
					size_t idxCorner1 = 2*iDual1 - 1 + cornerOffset[iSubCorner][1];
					
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
		for (size_t iCell0 = 1; iCell0 <= numCell0; ++iCell0) {
			for (size_t iCell1 = 1; iCell1 <= numCell1; ++iCell1) {

				for (size_t iSubCorner = 0; iSubCorner < 4; ++iSubCorner) {

					double corner_velocity[2];
					size_t idxFatherDual0 = iCell0 + cornerOffset[iSubCorner][0];
					size_t idxFatherDual1 = iCell1 + cornerOffset[iSubCorner][1];
					corner_velocity[0] = dualgrid.halftime_velocity[idxFatherDual0][idxFatherDual1][0];
					corner_velocity[1] = dualgrid.halftime_velocity[idxFatherDual0][idxFatherDual1][1];


					size_t idxCorner0 = 2*iCell0 + cornerOffset[iSubCorner][0];
					size_t idxCorner1 = 2*iCell1 + cornerOffset[iSubCorner][1];
					
					double cornerWork = \
							corners.force[idxCorner0][idxCorner1][0] * corner_velocity[0] \
						  + corners.force[idxCorner0][idxCorner1][1] * corner_velocity[1];

					grid.interior[iCell0][iCell1] += - tau * cornerWork / grid.mass[iCell0][iCell1];
			
				}
			}
		}

		// Update dual center (cell vertex) position
		for (size_t iDual0 = 1; iDual0 <= numDual0; ++iDual0) {
			for (size_t iDual1 = 1; iDual1 <= numDual1; ++iDual1) {
				size_t idxNode0 = dualgrid.center[iDual0][iDual1][0];
				size_t idxNode1 = dualgrid.center[iDual0][iDual1][1];
				coor_node[idxNode0][idxNode1][0] += tau * dualgrid.halftime_velocity[iDual0][iDual1][0];
				coor_node[idxNode0][idxNode1][1] += tau * dualgrid.halftime_velocity[iDual0][iDual1][1];
			}
		}

		// Update cell center position
// 似乎不会直接用到单元中心?

		

		// Update cell pressure
		for (size_t iCell0 = 1; iCell0 <= numCell0; ++iCell0) {
			for (size_t iCell1 = 1; iCell1 <= numCell1; ++iCell1) {
			
				double cellArea = CalCellArea(grid, coor_node, iCell0, iCell1);
				double cellDensity  = grid.mass[iCell0][iCell1] / cellArea;

				grid.pressure[iCell0][iCell1] = CalPressure(cellDensity, grid.interior[iCell0][iCell1], GAMMA);
			}
		}


		if (tcounter % (10*numCell0) == 0) {
			double E = CalTotalEnergy(grid, dualgrid, numCell0, numCell1, numDual0, numDual1);

			printf("t=%e E=%20.16e relerr=%20.16e\n",
			    	t,E,(E-E0)/E0);
		}

	}





// Test  diagnostics

double maxVelocity = 0.0;

for (size_t i = 1; i <= numDual0; ++i) {
    for (size_t j = 1; j <= numDual1; ++j) {

        double vx = dualgrid.velocity[i][j][0];
        double vy = dualgrid.velocity[i][j][1];

        double v = sqrt(vx*vx + vy*vy);

        if (v > maxVelocity)
            maxVelocity = v;
    }
}

printf("t = %.8e, max|v| = %.16e\n",
       t, maxVelocity);






    



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

	free(coor_node_data);
	free(coor_node);
	
}
