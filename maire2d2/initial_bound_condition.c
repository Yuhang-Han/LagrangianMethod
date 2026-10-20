/*****************************
* 2D Cartesian Sedov problem
*
* domain: [0,1]^2
* explosion center: (0.5, 0.5)
*
* boundary: constant 
******************************/

#include <math.h>

#include "linear_algebra_2d.h"
#include "initial_bound_condition.h"

const double EXPLOSION_CENTER[2] = {0.5, 0.5};

const double INITIAL_EXPLOSION_RADIUS = 0.1;
const double AMBIENT_PRESSURE = 1e-5;
const double EXPLOSION_ENERGY = 1.;


void initializeNodeCoor(size_t node_x, size_t node_y,
						double cell_h_0, double cell_h_1,
						size_t START_NODE_X, size_t START_NODE_Y,
						double coor_node[2])
{
	coor_node[0] = 0. + (node_x - START_NODE_X) * cell_h_0;
	coor_node[1] = 0. + (node_y - START_NODE_Y) * cell_h_1;
}


double initializeCellVolume(double cell_h_0, double cell_h_1)
{
	return cell_h_0 * cell_h_1;
}


double initializeCellDensity()
{
	return 1.;
}


double initializeCellMass(double volume_cell, double density_cell)
{
	return volume_cell * density_cell;
}


void initializeCellSpMomentum(double sp_momentum_cell[2])
{
	sp_momentum_cell[0] = 0.;
	sp_momentum_cell[1] = 0.;
}

double initializeCellPresure(size_t cell_x, size_t cell_y,
							double cell_h_0, double cell_h_1,
						    size_t NUM_GHOST_NODES_X, size_t NUM_GHOST_NODES_Y,
						    double coor_node[NUM_GHOST_NODES_Y][NUM_GHOST_NODES_X][2],
						    double GAMMA)
{
	double explosion_pressure = AMBIENT_PRESSURE \
								+ (GAMMA - 1) * EXPLOSION_ENERGY \
									/ (M_PI * pow(INITIAL_EXPLOSION_RADIUS, 2));

	double min_vertex[2] = {cell_h_0, cell_h_1};
	if( calNorm2(min_vertex) > INITIAL_EXPLOSION_RADIUS ) {
		printf("WARNING: Explotion radius is smaller than cell width.\n");
		printf("Using pressure/area to initialize.\n");

		explosion_pressure = explosion_pressure / (cell_h_0*cell_h_1);
	}
	
    /* cell_x/cell_y identify a physical cell with initialized vertices. */
    double coor_cell[2];

    /* Center of the initial Cartesian cell, using its four vertices. */
    for (size_t component = 0; component < 2; ++component) {
        coor_cell[component] = 0.25 * (
            coor_node[cell_y][cell_x][component] +
            coor_node[cell_y][cell_x + 1][component] +
            coor_node[cell_y + 1][cell_x + 1][component] +
            coor_node[cell_y + 1][cell_x][component]);
    }

    double radius_vector[2];
    vecSubVec2(coor_cell, EXPLOSION_CENTER, radius_vector);
    double distance_to_center = calNorm2(radius_vector);

    if (distance_to_center < INITIAL_EXPLOSION_RADIUS) {
        return explosion_pressure;
    }

    return AMBIENT_PRESSURE;
}



double initializeCellSpEnergy(double density_cell, double sp_momentum_cell[2], 
							  double pressure_cell, double GAMMA)
{
	double interior = pressure_cell / (GAMMA - 1) / density_cell;
	double sp_total_energy = interior \
						+ 0.5 * pow(calNorm2(sp_momentum_cell), 2);

	return sp_total_energy;
}



double boundarilizeMass(double cell_h_0, double cell_h_1)
{
	return cell_h_0 * cell_h_1 * 1.;
}

double boundarilizeDensity()
{
	return 1.;
}
	
	
double boundarilizePressure()
{
	return AMBIENT_PRESSURE;
}
	
	
double boundarilizeVolume(double cell_h_0, double cell_h_1)
{
	return cell_h_0 * cell_h_1;
}

double boundarilizeSpEnergy(double GAMMA) 
{
	double interior = AMBIENT_PRESSURE / (GAMMA - 1) / 1.;
	double sp_total_energy = interior + 0.;

	return sp_total_energy;
}

void boundarilizeSpMomentum(double sp_momentum_cell[2]) 
{
	sp_momentum_cell[0] = 0.;
	sp_momentum_cell[1] = 0.;
}


void boundarilizeCoor(size_t node_x, size_t node_y,
						double cell_h_0, double cell_h_1,
						size_t START_NODE_X, size_t START_NODE_Y,
						double coor_node[2])
{
	initializeNodeCoor(node_x, node_y, cell_h_0, cell_h_1, 
						START_NODE_X, START_NODE_Y, 
						coor_node);
}
