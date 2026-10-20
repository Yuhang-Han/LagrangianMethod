#include <stddef.h>


void initializeNodeCoor(size_t node_x, size_t node_y,
						double cell_h_0, double cell_h_1,
						size_t START_NODE_X, size_t START_NODE_Y,
						double coor_node[2]);

double initializeCellVolume(double cell_h_0, double cell_h_1);
double initializeCellDensity();
double initializeCellMass(double volume_cell, double density_cell);

void initializeCellSpMomentum(double sp_momentum_cell[2]);

double initializeCellPresure(size_t cell_x, size_t cell_y,
							double cell_h_0, double cell_h_1,
						    size_t NUM_GHOST_NODES_X, size_t NUM_GHOST_NODES_Y,
						    double coor_node[NUM_GHOST_NODES_Y][NUM_GHOST_NODES_X][2],
						    double GAMMA);

double initializeCellSpEnergy(double density_cell, double sp_momentum_cell[2], 
							  double pressure_cell, double GAMMA);

double boundarilizeMass(double cell_h_0, double cell_h_1);
double boundarilizeDensity();
double boundarilizePressure();
double boundarilizeVolume(double cell_h_0, double cell_h_1);
double boundarilizeSpEnergy(double GAMMA);

void boundarilizeSpMomentum(double sp_momentum_cell[2]);

void boundarilizeCoor(size_t node_x, size_t node_y,
						double cell_h_0, double cell_h_1,
						size_t START_NODE_X, size_t START_NODE_Y,
						double coor_node[2]);
