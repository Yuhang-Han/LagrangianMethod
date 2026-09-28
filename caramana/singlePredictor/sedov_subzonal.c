/* Modified from sedov_singlePredictor_checked.c.
   Adds Caramana--Shashkov subquadrilateral pressure forces:
   J. Comput. Phys. 142 (1998), 521--561, Eqs. (18)--(20),
   with the factor-of-two merit normalization stated on p. 537.
   Source: https://cnls.lanl.gov/~shashkov/papers/ed1.pdf

   The original single predictor, ideal-gas EOS, scalar artificial viscosity,
   fixed corner masses, free exterior boundary, and compatible work update
   are retained. This is not a reproduction of the published axisymmetric
   Sedov calculation. Changes beyond the force are documented in the diff.

   Build/run:
     mkdir -p build/output
     gcc -std=c11 -O2 sedov_subzonal.c -lm -o build/sedov_subzonal
     ./build/sedov_subzonal
*/

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


/* Signed polygon area: a concave quad is not automatically inverted.
   Never turn a negative area into a positive one with fabs().
   Subzone areas are checked independently in CalCellSubzonalState(). */
double CalCellArea(Grid grid, double (**coor_node)[2], size_t ic, size_t jc)
{
    double X[4][2];
    for (size_t k = 0; k < 4; ++k) {
        size_t ni = grid.vertices[ic][jc][k][0];
        size_t nj = grid.vertices[ic][jc][k][1];
        X[k][0] = coor_node[ni][nj][0];
        X[k][1] = coor_node[ni][nj][1];
    }
    double twice_area = 0.0;
    for (size_t k = 1; k < 3; ++k)
        twice_area += (X[k][0]-X[0][0])*(X[k+1][1]-X[0][1])
                    - (X[k][1]-X[0][1])*(X[k+1][0]-X[0][0]);
    return 0.5 * twice_area;
}

double CalPressure(double density, double interior, double gamma) {
	return (gamma - 1) * density * interior;
}




/* Current subzone geometry/pressures. Corner masses remain Lagrangian.
   Reference: Caramana & Shashkov, JCP 142 (1998), pp. 529, 536-537.
   Vertex numbering is BL, BR, TR, TL (counterclockwise).
*/
typedef struct {
    double area[4];
    double pressure[4];
    double S[4][2];           /* oriented median-segment normal */
    double min_length;       /* cell area / longest edge: a practical CFL scale */
} SubzonalState;

static void CalCellSubzonalState(
    Grid grid, Corners corners, double (**coor_node)[2],
    size_t ic, size_t jc, double gamma, SubzonalState *sz)
{
    double X[4][2], B[4][2], O[2] = {0.0, 0.0};
    double longest_edge = 0.0;
    const double e = grid.interior[ic][jc];

    if (!(e >= 0.0) || !isfinite(e) || !(gamma > 1.0)) {
        fprintf(stderr, "Invalid internal energy/EOS in cell (%zu,%zu): e=%g\n",
                ic, jc, e);
        exit(EXIT_FAILURE);
    }
    for (size_t k = 0; k < 4; ++k) {
        size_t ni = grid.vertices[ic][jc][k][0];
        size_t nj = grid.vertices[ic][jc][k][1];
        for (size_t d = 0; d < 2; ++d) {
            X[k][d] = coor_node[ni][nj][d];
            O[d] += 0.25 * X[k][d];
        }
    }
    for (size_t k = 0; k < 4; ++k) {
        size_t kp = (k + 1) % 4;
        for (size_t d = 0; d < 2; ++d)
            B[k][d] = 0.5 * (X[k][d] + X[kp][d]);
        longest_edge = fmax(longest_edge,
                            hypot(X[kp][0] - X[k][0], X[kp][1] - X[k][1]));
    }

    double total_area = 0.0;
    for (size_t k = 0; k < 4; ++k) {
        size_t km = (k + 3) % 4;
        size_t I = 2*ic + cornerOffset[k][0];
        size_t J = 2*jc + cornerOffset[k][1];
        /* Signed area of [X_k, B_k, O, B_{k-1}].
           Two triangle cross products, with X_k as local origin. */
        double ax = B[k][0] - X[k][0], ay = B[k][1] - X[k][1];
        double bx = O[0] - X[k][0],    by = O[1] - X[k][1];
        double cx = B[km][0] - X[k][0], cy = B[km][1] - X[k][1];
        sz->area[k] = 0.5 * (ax*by - ay*bx + bx*cy - by*cx);
        if (!(sz->area[k] > 0.0) || !isfinite(sz->area[k]) ||
            !(corners.mass[I][J] > 0.0)) {
            fprintf(stderr, "Invalid subzone (%zu,%zu,%zu): A=%g, m=%g\n",
                    ic, jc, k, sz->area[k], corners.mass[I][J]);
            exit(EXIT_FAILURE);
        }
        total_area += sz->area[k];
        sz->pressure[k] = CalPressure(corners.mass[I][J] / sz->area[k], e, gamma);
        if (!isfinite(sz->pressure[k])) {
            fprintf(stderr, "Nonfinite subzonal pressure in cell (%zu,%zu)\n", ic, jc);
            exit(EXIT_FAILURE);
        }
        /* S_k = rotate_clockwise(B_k - O).
           With this orientation, C_k = S_k - S_{k-1}. */
        sz->S[k][0] = B[k][1] - O[1];
        sz->S[k][1] = O[0] - B[k][0];
    }
    sz->min_length = total_area / longest_edge;
    if (!(sz->min_length > 0.0) || !isfinite(sz->min_length)) {
        fprintf(stderr, "Invalid cell length in cell (%zu,%zu)\n", ic, jc);
        exit(EXIT_FAILURE);
    }
}

static void AddCellSubzonalForces(
    const SubzonalState *sz, double pzone, double merit, double F[4][2])
{
    /* Preferred subquadrilateral differencing, Eq. (20), multiplied by
       2*merit as specified on p. 537. Equal/opposite edge contributions
       conserve zone momentum; do NOT use p_corner*C_corner instead. */
    for (size_t k = 0; k < 4; ++k) {
        size_t kp = (k + 1) % 4;
        double dp_k  = sz->pressure[k]  - pzone;
        double dp_kp = sz->pressure[kp] - pzone;
        double weight = merit * (dp_k + dp_kp);
        for (size_t d = 0; d < 2; ++d) {
            double g = weight * sz->S[k][d];
            F[k][d]  += g;
            F[kp][d] -= g;
        }
    }
}


double InitializeCellAveDensity(size_t idxCell0, size_t idxCell1) {
	return 1.0;
}

double InitializeDualAveMomentum(size_t idxCell0, size_t idxCell1) {
	return 0.0;
}

double InitializeCellInterior(size_t idxCell0, size_t idxCell1, size_t numCell0, size_t numCell1, double mass) {
	const double e_floor = 1.0e-12;
	const double E_blast = 1.0;

	// hot cell
    if (idxCell0 == (numCell0 + 1)/2 && idxCell1 == (numCell1 + 1)/2)
        return e_floor + E_blast / mass;

    return e_floor;
}

double InitializeCellPressure(size_t idxCell0, size_t idxCell1, double interior, double GAMMA) {
	double initial_density = InitializeCellAveDensity(idxCell0, idxCell1);
    return CalPressure(initial_density, interior, GAMMA);
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

void ApplyVacuumBoundary(
    Grid grid,
    size_t numGhCell0,
    size_t numGhCell1)
{
    for (size_t i = 0; i < numGhCell0; ++i) {
        for (size_t j = 0; j < numGhCell1; ++j) {

            int physical =
                (i >= 1 && i < numGhCell0-1 &&
                 j >= 1 && j < numGhCell1-1);

            if (!physical) {
                grid.mass[i][j] = 0.0;
                grid.pressure[i][j] = 0.0;
                grid.interior[i][j] = 0.0;
            }
        }
    }
}


double CalCellDivergence(
    Grid grid,
    DualGrid dualgrid,
    double (**coor_node)[2],
    size_t ic,
    size_t jc)
{
    double C[4][2];

    CalCellCornerVectors(grid, coor_node, ic, jc, C);

    double area = CalCellArea(grid, coor_node, ic, jc);

    double sum = 0.0;

    for (size_t k = 0; k < 4; ++k) {

        size_t iDual0 = ic + cornerOffset[k][0];
        size_t iDual1 = jc + cornerOffset[k][1];

        double vx = dualgrid.velocity[iDual0][iDual1][0];
        double vy = dualgrid.velocity[iDual0][iDual1][1];

        sum += C[k][0] * vx
             + C[k][1] * vy;
    }

    return sum / area;
}


double CalArtificialPressure(
    Grid grid,
    DualGrid dualgrid,
    double (**coor_node)[2],
    size_t ic,
    size_t jc,
    double Cq)
{
    double area = CalCellArea(grid, coor_node, ic, jc);

    double rho = grid.mass[ic][jc] / area;

    double divv =
        CalCellDivergence(grid, dualgrid,
                          coor_node, ic, jc);

    if (divv >= 0.0)
        return 0.0;

    double l = sqrt(area);

    return Cq * rho * l*l * divv*divv;
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


static void WriteCellCSV(
    const char *filename,
    double t,
    Grid grid,
    DualGrid dualgrid,
    double (**coor_node)[2],
    size_t numCell0,
    size_t numCell1,
    double Cq)
{
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        perror(filename);
        exit(EXIT_FAILURE);
    }

    fprintf(fp, "time,i,j,x,y,area,rho,pressure,q,ptotal,internal,divv\n");

    for (size_t i = 1; i <= numCell0; ++i) {
        for (size_t j = 1; j <= numCell1; ++j) {
            double xc = 0.0, yc = 0.0;
            for (size_t k = 0; k < 4; ++k) {
                size_t ni = grid.vertices[i][j][k][0];
                size_t nj = grid.vertices[i][j][k][1];
                xc += 0.25 * coor_node[ni][nj][0];
                yc += 0.25 * coor_node[ni][nj][1];
            }

            double area = CalCellArea(grid, coor_node, i, j);
            if (!(area >= 0.0) || !isfinite(area)) {
                fprintf(stderr, "Invalid cell area while writing CSV: (%zu,%zu), area=%g\n", i, j, area);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            double rho = grid.mass[i][j] / area;
            double divv = CalCellDivergence(grid, dualgrid, coor_node, i, j);
            double q = CalArtificialPressure(grid, dualgrid, coor_node, i, j, Cq);

            fprintf(fp,
                    "%.17e,%zu,%zu,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e\n",
                    t, i, j, xc, yc, area, rho,
                    grid.pressure[i][j], q, grid.pressure[i][j] + q,
                    grid.interior[i][j], divv);
        }
    }

    fclose(fp);
}

static void WriteNodeCSV(
    const char *filename,
    double t,
    DualGrid dualgrid,
    double (**coor_node)[2],
    size_t numDual0,
    size_t numDual1)
{
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        perror(filename);
        exit(EXIT_FAILURE);
    }

    fprintf(fp, "time,i,j,x,y,vx,vy,mass\n");

    for (size_t i = 1; i <= numDual0; ++i) {
        for (size_t j = 1; j <= numDual1; ++j) {
            size_t ni = dualgrid.center[i][j][0];
            size_t nj = dualgrid.center[i][j][1];

            fprintf(fp,
                    "%.17e,%zu,%zu,%.17e,%.17e,%.17e,%.17e,%.17e\n",
                    t, i, j,
                    coor_node[ni][nj][0], coor_node[ni][nj][1],
                    dualgrid.velocity[i][j][0], dualgrid.velocity[i][j][1],
                    dualgrid.mass[i][j]);
        }
    }

    fclose(fp);
}

static void InitHistoryCSV(const char *filename)
{
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        perror(filename);
        exit(EXIT_FAILURE);
    }
    fprintf(fp, "step,time,dt,total_energy,relative_energy_error,max_velocity,min_area\n");
    fclose(fp);
}

static void AppendHistoryCSV(
    const char *filename,
    size_t step,
    double t,
    double dt,
    double E0,
    Grid grid,
    DualGrid dualgrid,
    double (**coor_node)[2],
    size_t numCell0,
    size_t numCell1,
    size_t numDual0,
    size_t numDual1)
{
    double E = CalTotalEnergy(grid, dualgrid,
                              numCell0, numCell1,
                              numDual0, numDual1);

    double maxv = 0.0;
    for (size_t i = 1; i <= numDual0; ++i) {
        for (size_t j = 1; j <= numDual1; ++j) {
            double vx = dualgrid.velocity[i][j][0];
            double vy = dualgrid.velocity[i][j][1];
            double v = sqrt(vx*vx + vy*vy);
            if (v > maxv) maxv = v;
        }
    }

    double minArea = 1.0e300;
    for (size_t i = 1; i <= numCell0; ++i) {
        for (size_t j = 1; j <= numCell1; ++j) {
            double area = CalCellArea(grid, coor_node, i, j);
            if (area < minArea) minArea = area;
        }
    }

    FILE *fp = fopen(filename, "a");
    if (!fp) {
        perror(filename);
        exit(EXIT_FAILURE);
    }
    fprintf(fp, "%zu,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e\n",
            step, t, dt, E, (E-E0)/E0, maxv, minArea);
    fclose(fp);
}


int main()
{

	const double GAMMA = 5.0/3.0;

	const double T = 1.0;

	const double CFL = 0.2;
	const double CQ = 2.0;
    const double MF = 0.5;   /* tested on this 11x11 setup; 0 disables only the force */
	const size_t OUTPUT_EVERY = 1;

	const double L = 0.2;
	const double X0 = 2.0*L;
	const double X1 = 2.0*L;

	const size_t numCell0 = 101;
	const size_t numCell1 = 101;

	
	const size_t numGhCell0 = numCell0 + 2;
	const size_t numGhCell1 = numCell1 + 2;

	const size_t numDual0 = numCell0 + 1;
	const size_t numDual1 = numCell1 + 1;
	const size_t numGhDual0 = numDual0 + 2;
	const size_t numGhDual1 = numDual1 + 2;

	const size_t numGhCorner0 = 2 * numGhCell0;
	const size_t numGhCorner1 = 2 * numGhCell1;

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
	double (*coor_node_data)  [2] 	= calloc(numGhNode0   *	numGhNode1, 	sizeof(*coor_node_data));
	
	for (size_t i = 0; i < numGhNode0; i++) {
		coor_node[i] 	= coor_node_data 	+ i * numGhNode1;
	}


/*****************************
* Initialize
******************************/

// (uniform) Grid geometry

	// node and center coordinates
	double cellH0 = X0 / numCell0;
	double cellH1 = X1 / numCell1;
	
	for (size_t iGhNode0 = 0; iGhNode0 < numGhNode0; ++iGhNode0) {
		for (size_t iGhNode1 = 0; iGhNode1 < numGhNode1; ++iGhNode1) {
			coor_node[iGhNode0][iGhNode1][0] = -L + ((double)iGhNode0 -1.0) * cellH0;
			coor_node[iGhNode0][iGhNode1][1] = -L + ((double)iGhNode1 -1.0) * cellH1;
		}
	}

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

			if (idxFatherCell0 >= 1 && idxFatherCell0 <= numCell0 &&
			    idxFatherCell1 >= 1 && idxFatherCell1 <= numCell1) {

			    corners.mass[iGhCorner0][iGhCorner1]
			        = cornerArea * InitializeCellAveDensity(idxFatherCell0, idxFatherCell1);
			}
			else {
			    corners.mass[iGhCorner0][iGhCorner1] = 0.0;
			}
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

			// cell interior energy
			grid.interior[iGhCell0][iGhCell1] = InitializeCellInterior(iGhCell0, iGhCell1, 
												numCell0, numCell1, grid.mass[iGhCell0][iGhCell1]);
			// cell pressure
			grid.pressure[iGhCell0][iGhCell1] = InitializeCellPressure(iGhCell0, iGhCell1, 
												grid.interior[iGhCell0][iGhCell1], GAMMA);

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

	InitHistoryCSV("build/output/history.csv");
	AppendHistoryCSV("build/output/history.csv", 0, t, 0.0, E0,
	                 grid, dualgrid, coor_node,
	                 numCell0, numCell1, numDual0, numDual1);
	WriteCellCSV("build/output/cells_000000.csv", t, grid, dualgrid, coor_node,
	             numCell0, numCell1, CQ);
	WriteNodeCSV("build/output/nodes_000000.csv", t, dualgrid, coor_node,
	             numDual0, numDual1);

	while (t < T) {
		// Apply boundary condition
		ApplyVacuumBoundary(grid, numGhCell0, numGhCell1);
		// Determine time step length

		double tau = T - t;

		/* Practical CFL estimate including the subzonal pressure scale. */
		for (size_t ic = 1; ic <= numCell0; ++ic) {
			for (size_t jc = 1; jc <= numCell1; ++jc) {
				double area = CalCellArea(grid, coor_node, ic, jc);
				if (!(area > 0.0) || !isfinite(area)) {
					fprintf(stderr, "Mesh tangling/non-positive area at cell (%zu,%zu), t=%g\n", ic, jc, t);
					exit(EXIT_FAILURE);
				}

				double rho = grid.mass[ic][jc] / area;
				double q = CalArtificialPressure(grid, dualgrid, coor_node, ic, jc, CQ);
                SubzonalState sz;
                CalCellSubzonalState(grid, corners, coor_node, ic, jc, GAMMA, &sz);
                double pmax = grid.pressure[ic][jc];
                for (size_t k = 0; k < 4; ++k)
                    pmax = fmax(pmax, sz.pressure[k]);
                /* Practical estimate, not a proved stability bound.
                   Use corner pressure maxima and include the strength multiplier.
                   For MF=0 retain the original acoustic pressure estimate. */
                double pstar = grid.pressure[ic][jc] + q;
                if (MF > 0.0)
                    pstar = fmax(grid.pressure[ic][jc], fmax(1.0, 2.0*MF)*pmax) + q;
                double cstar = (pstar > 0.0) ? sqrt(GAMMA * pstar / rho) : 0.0;
                double l = sz.min_length;

				if (cstar > 0.0)
					tau = GetMinDouble(tau, CFL * l / cstar);

				double divv = CalCellDivergence(grid, dualgrid, coor_node, ic, jc);
				if (fabs(divv) > 0.0)
					tau = GetMinDouble(tau, 0.8 * CFL / fabs(divv));
			}
		}

		if (t+tau > T) {
			tau = T - t;
		}

		t += tau;
		++tcounter;
		
        // Form the COMPLETE corner force once, before changing velocity.
        memset(corners.force_data, 0,
               numGhCorner0 * numGhCorner1 * sizeof(*corners.force_data));
        for (size_t ic = 1; ic <= numCell0; ++ic) {
            for (size_t jc = 1; jc <= numCell1; ++jc) {
                double C[4][2], F[4][2];
                SubzonalState sz;
                CalCellCornerVectors(grid, coor_node, ic, jc, C);
                CalCellSubzonalState(grid, corners, coor_node, ic, jc, GAMMA, &sz);
                double q = CalArtificialPressure(grid, dualgrid, coor_node, ic, jc, CQ);
                double Ptotal = grid.pressure[ic][jc] + q;
                for (size_t k = 0; k < 4; ++k)
                    for (size_t d = 0; d < 2; ++d)
                        F[k][d] = Ptotal * C[k][d];

                AddCellSubzonalForces(&sz, grid.pressure[ic][jc], MF, F);
                for (size_t k = 0; k < 4; ++k) {
                    size_t I = 2*ic + cornerOffset[k][0];
                    size_t J = 2*jc + cornerOffset[k][1];
                    for (size_t d = 0; d < 2; ++d) {
                        corners.vector[I][J][d] = C[k][d];
                        corners.force[I][J][d] = F[k][d];
                    }
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
		

		// Update cell pressure
		for (size_t iCell0 = 1; iCell0 <= numCell0; ++iCell0) {
			for (size_t iCell1 = 1; iCell1 <= numCell1; ++iCell1) {
			
                /* Validate the advanced state, including the final timestep.
                   Stop rather than silently taking abs(area) or flooring e. */
                SubzonalState advanced;
                CalCellSubzonalState(grid, corners, coor_node,
                                    iCell0, iCell1, GAMMA, &advanced);
				double cellArea = CalCellArea(grid, coor_node, iCell0, iCell1);
                if (!(cellArea > 0.0) || !isfinite(cellArea)) {
                    fprintf(stderr, "Invalid advanced cell area (%zu,%zu), t=%g\n",
                            iCell0, iCell1, t);
                    exit(EXIT_FAILURE);
                }
				double cellDensity  = grid.mass[iCell0][iCell1] / cellArea;

				grid.pressure[iCell0][iCell1] = CalPressure(cellDensity, grid.interior[iCell0][iCell1], GAMMA);
			}
		}


		AppendHistoryCSV("build/output/history.csv", tcounter, t, tau, E0,
		                 grid, dualgrid, coor_node,
		                 numCell0, numCell1, numDual0, numDual1);

		if (tcounter % OUTPUT_EVERY == 0 || t >= T) {
			char cellName[128];
			char nodeName[128];
			snprintf(cellName, sizeof(cellName), "build/output/cells_%06zu.csv", tcounter);
			snprintf(nodeName, sizeof(nodeName), "build/output/nodes_%06zu.csv", tcounter);

			WriteCellCSV(cellName, t, grid, dualgrid, coor_node,
			             numCell0, numCell1, CQ);
			WriteNodeCSV(nodeName, t, dualgrid, coor_node,
			             numDual0, numDual1);
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

double E = CalTotalEnergy(grid, dualgrid, numCell0, numCell1, numDual0, numDual1);

printf("t = %.8e, max|v| = %.16e, E=%20.16e\n",
       t, maxVelocity, E);


printf("v_LD: (%f, %f) \n", dualgrid.velocity[5][5][0], dualgrid.velocity[5][5][1]);
printf("v_LU: (%f, %f) \n", dualgrid.velocity[6][5][0], dualgrid.velocity[6][5][1]);
printf("v_RD: (%f, %f) \n", dualgrid.velocity[5][6][0], dualgrid.velocity[5][6][1]);
printf("v_RU: (%f, %f) \n", dualgrid.velocity[6][6][0], dualgrid.velocity[6][6][1]);



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
