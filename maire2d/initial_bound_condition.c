/*****************************
* 2D Cartesian Sedov problem
*
* domain: [0,1]^2
* explosion center: (0.5, 0.5)
*
*
******************************/

#include <math.h>

#include "linear_algebra_2d.h"
#include "initial_bound_condition.h"

const double EXPLOSION_CENTER[2] = {0.5, 0.5};

const double INITIAL_EXPLOSION_RADIUS = 0.02;
const double AMBIENT_PRESSURE = 1e-5;
const double EXPLOSION_ENERGY = 1.;


void initializeNodeCoor(size_t idxP, size_t num_nodes_x, 
						double cell_h_0, double cell_h_1,
						double coor_node[2])
{
	size_t row_node = idxP / num_nodes_x;
	size_t col_node = idxP % num_nodes_x;

	coor_node[0] = 0. + col_node * cell_h_0;
	coor_node[1] = 0. + row_node * cell_h_1;
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


double initializeCellPresure(double coor_node[2], double GAMMA)
{

	
	double explosion_pressure = AMBIENT_PRESSURE \
								+ (GAMMA - 1) * EXPLOSION_ENERGY \
									/ (M_PI * pow(INITIAL_EXPLOSION_RADIUS, 2));

	double radius_vecter[2];
	vecSubVec2(coor_node, EXPLOSION_CENTER, radius_vecter);
	double distance_to_center = calNorm2(radius_vecter);

	if(distance_to_center < INITIAL_EXPLOSION_RADIUS) {
		return explosion_pressure;
	}
	else {
		return AMBIENT_PRESSURE;
	}
}


double initializeCellSpEnergy(double density_cell, double sp_momentum_cell[2], 
							  double pressure_cell, double GAMMA)
{
	double interior = pressure_cell / (GAMMA - 1) / density_cell;
	double sp_total_energy = interior \
						+ 0.5 * calNorm2(sp_momentum_cell);

	return sp_total_energy;
}
