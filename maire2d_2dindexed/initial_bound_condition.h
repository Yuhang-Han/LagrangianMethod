#ifndef INITIAL_BOUND_CONDITION_H
#define INITIAL_BOUND_CONDITION_H

#include <stddef.h>

/* node_x/node_y are actual array indices; START_NODE_* locates (0, 0). */
void initializeNodeCoor(size_t node_x, size_t node_y,
                        size_t START_NODE_X, size_t START_NODE_Y,
                        double cell_h_x, double cell_h_y,
                        double coor_node[2]);

/*
 * All bounds are inclusive. There must be at least two physical nodes in
 * each direction. Extend the current physical mesh linearly into every
 * exterior ghost layer, first in x and then in y, including all corners.
 * Physical nodes are never modified. Returns 1 on success, 0 for invalid
 * bounds or nonfinite geometry. On failure, ghost values may be partial.
 * The caller supplies positive array dimensions and matching storage.
 */
int fillGhostNodeCoordinates(
    size_t NUM_TOTAL_NODES_X, size_t NUM_TOTAL_NODES_Y,
    size_t START_NODE_X, size_t END_NODE_X,
    size_t START_NODE_Y, size_t END_NODE_Y,
    double coor_node[NUM_TOTAL_NODES_Y][NUM_TOTAL_NODES_X][2]);

double initializeCellVolume(double cell_h_x, double cell_h_y);
double initializeCellDensity(void);
double initializeCellMass(double volume_cell, double density_cell);
void initializeCellSpMomentum(double sp_momentum_cell[2]);

/* coor_cell is the position of the cell being initialized. */
double initializeCellPresure(double coor_cell[2], double GAMMA);
double initializeCellSpEnergy(double density_cell, double sp_momentum_cell[2],
                              double pressure_cell, double GAMMA);

double boundarilizeMass(double cell_h_x, double cell_h_y);
double boundarilizeDensity(void);
double boundarilizePressure(void);
double boundarilizeVolume(double cell_h_x, double cell_h_y);
double boundarilizeSpEnergy(double GAMMA);
void boundarilizeSpMomentum(double sp_momentum_cell[2]);

#endif
