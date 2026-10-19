/*****************************
* 2D Cartesian Sedov problem
*
* domain: [0,1]^2
* explosion center: (0.5, 0.5)
*
* boundary: constant ambient state, extrapolated ghost coordinates
******************************/

#include <math.h>

#include "linear_algebra_2d.h"
#include "initial_bound_condition.h"

const double EXPLOSION_CENTER[2] = {0.5, 0.5};

const double INITIAL_EXPLOSION_RADIUS = 0.02;
const double AMBIENT_PRESSURE = 1e-5;
const double EXPLOSION_ENERGY = 1.;


void initializeNodeCoor(size_t node_x, size_t node_y,
                        size_t START_NODE_X, size_t START_NODE_Y,
                        double cell_h_x, double cell_h_y,
                        double coor_node[2])
{
    /* Convert before subtraction: ghost coordinates can be negative. */
    coor_node[0] = ((double)node_x - (double)START_NODE_X) * cell_h_x;
    coor_node[1] = ((double)node_y - (double)START_NODE_Y) * cell_h_y;
}

int fillGhostNodeCoordinates(
    size_t NUM_TOTAL_NODES_X, size_t NUM_TOTAL_NODES_Y,
    size_t START_NODE_X, size_t END_NODE_X,
    size_t START_NODE_Y, size_t END_NODE_Y,
    double coor_node[NUM_TOTAL_NODES_Y][NUM_TOTAL_NODES_X][2])
{
    if (coor_node == NULL || NUM_TOTAL_NODES_X < 2 ||
        NUM_TOTAL_NODES_Y < 2 || START_NODE_X >= END_NODE_X ||
        START_NODE_Y >= END_NODE_Y || END_NODE_X >= NUM_TOTAL_NODES_X ||
        END_NODE_Y >= NUM_TOTAL_NODES_Y)
        return 0;

    /* Ghost values may be uninitialized; validate only physical inputs. */
    for (size_t node_y = START_NODE_Y; node_y <= END_NODE_Y; ++node_y) {
        for (size_t node_x = START_NODE_X; node_x <= END_NODE_X; ++node_x) {
            if (!isfinite(coor_node[node_y][node_x][0]) ||
                !isfinite(coor_node[node_y][node_x][1]))
                return 0;
        }
    }

    for (size_t node_y = START_NODE_Y; node_y <= END_NODE_Y; ++node_y) {
        for (size_t node_x = 0; node_x < START_NODE_X; ++node_x) {
            double layers = (double)(START_NODE_X - node_x);
            for (size_t component = 0; component < 2; ++component) {
                double boundary = coor_node[node_y][START_NODE_X][component];
                double interior = coor_node[node_y][START_NODE_X + 1][component];
                coor_node[node_y][node_x][component] =
                    boundary + layers * (boundary - interior);
            }
        }
        for (size_t node_x = END_NODE_X + 1;
             node_x < NUM_TOTAL_NODES_X; ++node_x) {
            double layers = (double)(node_x - END_NODE_X);
            for (size_t component = 0; component < 2; ++component) {
                double boundary = coor_node[node_y][END_NODE_X][component];
                double interior = coor_node[node_y][END_NODE_X - 1][component];
                coor_node[node_y][node_x][component] =
                    boundary + layers * (boundary - interior);
            }
        }
    }

    /* Include the x halos here so all four exterior corners are filled. */
    for (size_t node_y = 0; node_y < START_NODE_Y; ++node_y) {
        double layers = (double)(START_NODE_Y - node_y);
        for (size_t node_x = 0; node_x < NUM_TOTAL_NODES_X; ++node_x) {
            for (size_t component = 0; component < 2; ++component) {
                double boundary = coor_node[START_NODE_Y][node_x][component];
                double interior = coor_node[START_NODE_Y + 1][node_x][component];
                coor_node[node_y][node_x][component] =
                    boundary + layers * (boundary - interior);
            }
        }
    }
    for (size_t node_y = END_NODE_Y + 1;
         node_y < NUM_TOTAL_NODES_Y; ++node_y) {
        double layers = (double)(node_y - END_NODE_Y);
        for (size_t node_x = 0; node_x < NUM_TOTAL_NODES_X; ++node_x) {
            for (size_t component = 0; component < 2; ++component) {
                double boundary = coor_node[END_NODE_Y][node_x][component];
                double interior = coor_node[END_NODE_Y - 1][node_x][component];
                coor_node[node_y][node_x][component] =
                    boundary + layers * (boundary - interior);
            }
        }
    }

    for (size_t node_y = 0; node_y < NUM_TOTAL_NODES_Y; ++node_y) {
        for (size_t node_x = 0; node_x < NUM_TOTAL_NODES_X; ++node_x) {
            if (!isfinite(coor_node[node_y][node_x][0]) ||
                !isfinite(coor_node[node_y][node_x][1]))
                return 0;
        }
    }
    return 1;
}


double initializeCellVolume(double cell_h_x, double cell_h_y)
{
	return cell_h_x * cell_h_y;
}


double initializeCellDensity(void)
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


double initializeCellPresure(double coor_cell[2], double GAMMA)
{

	
	double explosion_pressure = AMBIENT_PRESSURE \
								+ (GAMMA - 1) * EXPLOSION_ENERGY \
									/ (acos(-1.0) * pow(INITIAL_EXPLOSION_RADIUS, 2));

	double radius_vecter[2];
	vecSubVec2(coor_cell, EXPLOSION_CENTER, radius_vecter);
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



double boundarilizeMass(double cell_h_x, double cell_h_y)
{
	return cell_h_x * cell_h_y * 1.;
}

double boundarilizeDensity(void)
{
	return 1.;
}
	
	
double boundarilizePressure(void)
{
	return AMBIENT_PRESSURE;
}
	
	
double boundarilizeVolume(double cell_h_x, double cell_h_y)
{
	return cell_h_x * cell_h_y;
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
