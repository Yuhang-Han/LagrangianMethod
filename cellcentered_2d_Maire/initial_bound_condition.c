/*
 * Cartesian Sedov test: Maire's HDR, printed pp. 109 and 112.
 * Domain [0, 1.2]^2, 30 x 30 cells, rho=1, p=1e-6, u=0, gamma=7/5.
 * Deposit E0=0.244816 in the single physical cell at the origin.
 *
 * Boundary interpretation: reflect across the coordinate axes and retain
 * ambient states beyond the right/top boundaries. The exact shock reaches
 * R=1 at t=1, inside this domain. First-order acoustic GLACE can still fold
 * the mesh on this test, as reported on printed p. 112.
 */

#include "initial_bound_condition.h"

const double SEDOV_GAMMA = 7.0 / 5.0;
const double SEDOV_FINAL_TIME = 1.0;
const double SEDOV_DOMAIN_LENGTH = 1.2;
const double SEDOV_AMBIENT_DENSITY = 1.0;
const double SEDOV_AMBIENT_PRESSURE = 1e-6;
const double SEDOV_EXPLOSION_ENERGY = 0.244816;

void initializeNodeCoor(size_t node_x, size_t node_y,
                        double cell_h_0, double cell_h_1,
                        size_t START_NODE_X, size_t START_NODE_Y,
                        double coor_node[2])
{
    /* Convert before subtracting: size_t would wrap for lower ghost indices. */
    coor_node[0] = ((double)node_x - (double)START_NODE_X) * cell_h_0;
    coor_node[1] = ((double)node_y - (double)START_NODE_Y) * cell_h_1;
}


double initializeCellVolume(double cell_h_0, double cell_h_1)
{
    return cell_h_0 * cell_h_1;
}


double initializeCellDensity(void)
{
    return SEDOV_AMBIENT_DENSITY;
}


double initializeCellMass(double volume_cell, double density_cell)
{
    return volume_cell * density_cell;
}


void initializeCellSpMomentum(double sp_momentum_cell[2])
{
    sp_momentum_cell[0] = 0.0;
    sp_momentum_cell[1] = 0.0;
}


double initializeCellPresure(size_t cell_x, size_t cell_y,
                             double cell_h_0, double cell_h_1,
                             size_t NUM_TOTAL_NODES_X, size_t NUM_TOTAL_NODES_Y,
                             double coor_node[NUM_TOTAL_NODES_Y][NUM_TOTAL_NODES_X][2],
                             double GAMMA)
{
    /* initializeNodeCoor places the southwest physical vertex exactly at 0. */
    if (coor_node[cell_y][cell_x][0] == 0.0 &&
        coor_node[cell_y][cell_x][1] == 0.0) {
        const double volume = initializeCellVolume(cell_h_0, cell_h_1);

        /* Literal p_or=(gamma-1)*rho_or*E0/v_or from p. 109.
         * E0 is used directly in this quadrant: no extra factor of 1/4,
         * no finite-radius disk, and no added ambient pressure.
         */
        return (GAMMA - 1.0) * initializeCellDensity()
               * SEDOV_EXPLOSION_ENERGY / volume;
    }

    return SEDOV_AMBIENT_PRESSURE;
}


double initializeCellSpEnergy(double density_cell,
                              const double sp_momentum_cell[2],
                              double pressure_cell, double GAMMA)
{
    const double internal = pressure_cell / ((GAMMA - 1.0) * density_cell);
    const double kinetic = 0.5 * (sp_momentum_cell[0] * sp_momentum_cell[0]
                                + sp_momentum_cell[1] * sp_momentum_cell[1]);
    return internal + kinetic;
}


double boundarilizeMass(double cell_h_0, double cell_h_1)
{
    return initializeCellMass(boundarilizeVolume(cell_h_0, cell_h_1),
                              boundarilizeDensity());
}


double boundarilizeDensity(void)
{
    return SEDOV_AMBIENT_DENSITY;
}


double boundarilizePressure(void)
{
    return SEDOV_AMBIENT_PRESSURE;
}


double boundarilizeVolume(double cell_h_0, double cell_h_1)
{
    return initializeCellVolume(cell_h_0, cell_h_1);
}


double boundarilizeSpEnergy(double GAMMA)
{
    const double velocity[2] = {0.0, 0.0};
    return initializeCellSpEnergy(boundarilizeDensity(), velocity,
                                  boundarilizePressure(), GAMMA);
}


void boundarilizeSpMomentum(double sp_momentum_cell[2])
{
    initializeCellSpMomentum(sp_momentum_cell);
}


void boundarilizeCoor(size_t node_x, size_t node_y,
                      double cell_h_0, double cell_h_1,
                      size_t START_NODE_X, size_t START_NODE_Y,
                      double coor_node[2])
{
    initializeNodeCoor(node_x, node_y, cell_h_0, cell_h_1,
                        START_NODE_X, START_NODE_Y, coor_node);
}


void applySedovCellBoundaries(size_t num_cells_x, size_t num_cells_y,
                              double cell_h_0, double cell_h_1, double GAMMA,
                              double mass_cell[num_cells_y + 2][num_cells_x + 2],
                              double density_cell[num_cells_y + 2][num_cells_x + 2],
                              double pressure_cell[num_cells_y + 2][num_cells_x + 2],
                              double volume_cell[num_cells_y + 2][num_cells_x + 2],
                              double sp_energy_cell[num_cells_y + 2][num_cells_x + 2],
                              double sp_momentum_cell[num_cells_y + 2][num_cells_x + 2][2])
{
    for (size_t y = 0; y <= num_cells_y + 1; ++y) {
        for (size_t x = 0; x <= num_cells_x + 1; ++x) {
            if (x >= 1 && x <= num_cells_x &&
                y >= 1 && y <= num_cells_y)
                continue;

            if (x == num_cells_x + 1 || y == num_cells_y + 1) {
                /* Outer ambient state, including its two mixed corners. */
                mass_cell[y][x] = boundarilizeMass(cell_h_0, cell_h_1);
                density_cell[y][x] = boundarilizeDensity();
                pressure_cell[y][x] = boundarilizePressure();
                volume_cell[y][x] = boundarilizeVolume(cell_h_0, cell_h_1);
                sp_energy_cell[y][x] = boundarilizeSpEnergy(GAMMA);
                boundarilizeSpMomentum(sp_momentum_cell[y][x]);
            } else {
                /* Even thermodynamic state, odd normal velocity at each axis. */
                const size_t source_x = x == 0 ? 1 : x;
                const size_t source_y = y == 0 ? 1 : y;
                mass_cell[y][x] = mass_cell[source_y][source_x];
                density_cell[y][x] = density_cell[source_y][source_x];
                pressure_cell[y][x] = pressure_cell[source_y][source_x];
                volume_cell[y][x] = volume_cell[source_y][source_x];
                sp_energy_cell[y][x] = sp_energy_cell[source_y][source_x];
                sp_momentum_cell[y][x][0] = (x == 0 ? -1.0 : 1.0)
                    * sp_momentum_cell[source_y][source_x][0];
                sp_momentum_cell[y][x][1] = (y == 0 ? -1.0 : 1.0)
                    * sp_momentum_cell[source_y][source_x][1];
            }
        }
    }
}


void applySedovNodeBoundaries(size_t num_cells_x, size_t num_cells_y,
                              double coor_node[num_cells_y + 3][num_cells_x + 3][2])
{
    const size_t right = num_cells_x + 1;
    const size_t top = num_cells_y + 1;

    for (size_t y = 1; y <= top; ++y) {
        /* Reflect the current left interior node, preserving tangential motion. */
        coor_node[y][0][0] = -coor_node[y][2][0];
        coor_node[y][0][1] = coor_node[y][2][1];

        for (size_t component = 0; component < 2; ++component) {
            coor_node[y][right + 1][component] =
                2.0 * coor_node[y][right][component]
                - coor_node[y][right - 1][component];
        }
    }

    /* Include the new side ghosts so corner operations compose consistently. */
    for (size_t x = 0; x <= right + 1; ++x) {
        coor_node[0][x][0] = coor_node[2][x][0];
        coor_node[0][x][1] = -coor_node[2][x][1];

        for (size_t component = 0; component < 2; ++component) {
            coor_node[top + 1][x][component] =
                2.0 * coor_node[top][x][component]
                - coor_node[top - 1][x][component];
        }
    }
}


void applySedovVelocityBoundary(size_t node_x, size_t node_y,
                                size_t START_NODE_X, size_t START_NODE_Y,
                                double velocity_node[2])
{
    if (node_x == START_NODE_X)
        velocity_node[0] = 0.0;
    if (node_y == START_NODE_Y)
        velocity_node[1] = 0.0;
}
