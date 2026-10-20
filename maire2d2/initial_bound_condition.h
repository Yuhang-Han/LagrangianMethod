#ifndef INITIAL_BOUND_CONDITION_H
#define INITIAL_BOUND_CONDITION_H

#include <stddef.h>

/* Maire's HDR, printed pp. 109 and 112: Cartesian quadrant Sedov test. */
enum { SEDOV_NUM_CELLS = 30 };
extern const double SEDOV_GAMMA;
extern const double SEDOV_FINAL_TIME;
extern const double SEDOV_DOMAIN_LENGTH;
extern const double SEDOV_AMBIENT_DENSITY;
extern const double SEDOV_AMBIENT_PRESSURE;
extern const double SEDOV_EXPLOSION_ENERGY;

/* START_NODE_X/Y is the physical node at the origin; lower ghosts are negative. */
void initializeNodeCoor(size_t node_x, size_t node_y,
                        double cell_h_0, double cell_h_1,
                        size_t START_NODE_X, size_t START_NODE_Y,
                        double coor_node[2]);

double initializeCellVolume(double cell_h_0, double cell_h_1);
double initializeCellDensity(void);
double initializeCellMass(double volume_cell, double density_cell);
void initializeCellSpMomentum(double sp_momentum_cell[2]);

/*
 * Initial Cartesian physical cells only: all four vertices must be initialized.
 * The cell with southwest vertex (0, 0) receives the entire source energy.
 * The existing spelling of Presure is retained for compatibility.
 */
double initializeCellPresure(size_t cell_x, size_t cell_y,
                             double cell_h_0, double cell_h_1,
                             size_t NUM_TOTAL_NODES_X, size_t NUM_TOTAL_NODES_Y,
                             double coor_node[NUM_TOTAL_NODES_Y][NUM_TOTAL_NODES_X][2],
                             double GAMMA);

/* Specific momentum is velocity; this returns internal plus kinetic energy. */
double initializeCellSpEnergy(double density_cell,
                              const double sp_momentum_cell[2],
                              double pressure_cell, double GAMMA);

/* Constant ambient state for the right/top exterior, not symmetry boundaries. */
double boundarilizeMass(double cell_h_0, double cell_h_1);
double boundarilizeDensity(void);
double boundarilizePressure(void);
double boundarilizeVolume(double cell_h_0, double cell_h_1);
double boundarilizeSpEnergy(double GAMMA);
void boundarilizeSpMomentum(double sp_momentum_cell[2]);

/* Initial Cartesian ghost coordinates; use applySedovNodeBoundaries in time. */
void boundarilizeCoor(size_t node_x, size_t node_y,
                      double cell_h_0, double cell_h_1,
                      size_t START_NODE_X, size_t START_NODE_Y,
                      double coor_node[2]);

/*
 * One ghost layer: physical cells [1..num_cells_y][1..num_cells_x],
 * physical nodes [1..num_cells_y+1][1..num_cells_x+1].
 * Counts must be positive, spacings positive, and GAMMA > 1.
 *
 * The Sedov paragraph does not prescribe boundaries explicitly. We use quadrant
 * symmetry at x=0 and y=0 and undisturbed ambient ghost states at right/top.
 * Reflect CURRENT cell states and node geometry before every nodal solve.
 * Physical entries are left unchanged. At mixed corners, ambient cell data
 * take precedence; node reflection and outer extrapolation are composed.
 */
void applySedovCellBoundaries(size_t num_cells_x, size_t num_cells_y,
                              double cell_h_0, double cell_h_1, double GAMMA,
                              double mass_cell[num_cells_y + 2][num_cells_x + 2],
                              double density_cell[num_cells_y + 2][num_cells_x + 2],
                              double pressure_cell[num_cells_y + 2][num_cells_x + 2],
                              double volume_cell[num_cells_y + 2][num_cells_x + 2],
                              double sp_energy_cell[num_cells_y + 2][num_cells_x + 2],
                              double sp_momentum_cell[num_cells_y + 2][num_cells_x + 2][2]);

void applySedovNodeBoundaries(size_t num_cells_x, size_t num_cells_y,
                              double coor_node[num_cells_y + 3][num_cells_x + 3][2]);

/* Apply after the nodal solve: constrain only the normal velocity on each axis. */
void applySedovVelocityBoundary(size_t node_x, size_t node_y,
                                size_t START_NODE_X, size_t START_NODE_Y,
                                double velocity_node[2]);

#endif
