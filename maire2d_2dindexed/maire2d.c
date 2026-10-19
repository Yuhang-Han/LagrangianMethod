#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "linear_algebra_2d.h"
#include "mesh_calculation.h"
#include "physical_calculation.h"
#include "initial_bound_condition.h"
#include "data_output.h"

/* Arithmetic helpers take actual array elements, never converted mesh IDs.
 * next/prev follow the cell counterclockwise. Keep the original inward normal. */
static int accumulateNodeCell(
    double coor_node[2], double coor_next[2], double coor_prev[2],
    double density_cell, double pressure_cell, double sp_momentum_cell[2],
    double GAMMA, double matrix_node[2][2],
    double pressure_gradient[2], double ave_velocity_node[2])
{
    double length_corner, normal_corner[2];
    if (!calCornerGeometry(coor_node, coor_next, coor_prev,
                           &length_corner, normal_corner))
        return 0;
    double impedance = calImpendance(density_cell, pressure_cell, GAMMA);
    if (!isfinite(impedance) || impedance <= 0.0)
        return 0;

    double matrix_corner[2][2];
    calMatCorner_glace(length_corner, normal_corner, impedance, matrix_corner);
    double contrib_pressure_gradient[2];
    numMulVec2(length_corner * pressure_cell, normal_corner,
               contrib_pressure_gradient);
    attVec2(pressure_gradient, contrib_pressure_gradient);
    double contrib_ave_velocity_node[2];
    matMulVec2(matrix_corner, sp_momentum_cell, contrib_ave_velocity_node);
    attVec2(ave_velocity_node, contrib_ave_velocity_node);
    attMat2(matrix_node, matrix_corner);
    return 1;
}

static int accumulateCellCorner(
    double coor_node[2], double coor_next[2], double coor_prev[2],
    double impedance, double pressure_cell, double sp_momentum_cell[2],
    double velocity_node[2], double *volume_flux,
    double resul_force[2], double *expan_power)
{
    double length_corner, normal_corner[2];
    if (!calCornerGeometry(coor_node, coor_next, coor_prev,
                           &length_corner, normal_corner))
        return 0;
    double corner_force[2];
    calCornerForce(corner_force, impedance, sp_momentum_cell, pressure_cell,
                   length_corner, normal_corner, velocity_node);
    *volume_flux += length_corner * vecMulVec2(normal_corner, velocity_node);
    attVec2(resul_force, corner_force);
    *expan_power += vecMulVec2(corner_force, velocity_node);
    return 1;
}

/* A function allows the same solver to be checked on rectangular meshes.
 * The editable example parameters remain together in main(). */
static int runSimulation(const size_t NUM_CELLS_X,
                         const size_t NUM_CELLS_Y, const double T)
{
    const size_t OUTPUT_FREQUENCY = 100;
    const double GAMMA = 1.4;
    const double LENGTH_X = 1.0;
    const double LENGTH_Y = 1.0;
    const size_t NUM_GHOST_LAYERS = 1;

    /* Space for a ghost cell on each side and the extra node. */
    if (NUM_CELLS_X == 0 || NUM_CELLS_Y == 0 ||
        NUM_CELLS_X > SIZE_MAX - 3 || NUM_CELLS_Y > SIZE_MAX - 3 ||
        !isfinite(T) || T < 0.0) {
        fprintf(stderr, "Invalid physical mesh dimensions or final time.\n");
        return EXIT_FAILURE;
    }

    const size_t NUM_NODES_X = NUM_CELLS_X + 1;
    const size_t NUM_NODES_Y = NUM_CELLS_Y + 1;
    const size_t NUM_TOTAL_CELLS_X = NUM_CELLS_X + 2 * NUM_GHOST_LAYERS;
    const size_t NUM_TOTAL_CELLS_Y = NUM_CELLS_Y + 2 * NUM_GHOST_LAYERS;
    const size_t NUM_TOTAL_NODES_X = NUM_TOTAL_CELLS_X + 1;
    const size_t NUM_TOTAL_NODES_Y = NUM_TOTAL_CELLS_Y + 1;

    /***********************************************************************
     * All START/END bounds are INCLUSIVE actual array subscripts.
     * Field access: field[y][x]. Physical boundary nodes are physical nodes.
     * Ghost nodes lie outside START_NODE_* ... END_NODE_*.
     ***********************************************************************/
    const size_t START_CELL_X = NUM_GHOST_LAYERS;
    const size_t END_CELL_X = START_CELL_X + NUM_CELLS_X - 1;
    const size_t START_CELL_Y = NUM_GHOST_LAYERS;
    const size_t END_CELL_Y = START_CELL_Y + NUM_CELLS_Y - 1;
    const size_t START_NODE_X = NUM_GHOST_LAYERS;
    const size_t END_NODE_X = START_NODE_X + NUM_NODES_X - 1;
    const size_t START_NODE_Y = NUM_GHOST_LAYERS;
    const size_t END_NODE_Y = START_NODE_Y + NUM_NODES_Y - 1;

    const size_t START_TOTAL_CELL_X = 0;
    const size_t END_TOTAL_CELL_X = NUM_TOTAL_CELLS_X - 1;
    const size_t START_TOTAL_CELL_Y = 0;
    const size_t END_TOTAL_CELL_Y = NUM_TOTAL_CELLS_Y - 1;
    const size_t START_TOTAL_NODE_X = 0;
    const size_t END_TOTAL_NODE_X = NUM_TOTAL_NODES_X - 1;
    const size_t START_TOTAL_NODE_Y = 0;
    const size_t END_TOTAL_NODE_Y = NUM_TOTAL_NODES_Y - 1;

    /* Check the largest field before evaluating variable-length array types. */
    if (NUM_TOTAL_NODES_X > SIZE_MAX / (2 * sizeof(double)) ||
        NUM_TOTAL_NODES_Y >
            SIZE_MAX / (NUM_TOTAL_NODES_X * 2 * sizeof(double))) {
        fprintf(stderr, "Mesh allocation size overflows size_t.\n");
        return EXIT_FAILURE;
    }

    /* Contiguous heap storage with direct multidimensional array syntax. */
    double (*coor_node)[NUM_TOTAL_NODES_X][2] =
        calloc(NUM_TOTAL_NODES_Y, sizeof *coor_node);
    double (*velocity_node)[NUM_TOTAL_NODES_X][2] =
        calloc(NUM_TOTAL_NODES_Y, sizeof *velocity_node);
    double (*mass_cell)[NUM_TOTAL_CELLS_X] =
        calloc(NUM_TOTAL_CELLS_Y, sizeof *mass_cell);
    double (*density_cell)[NUM_TOTAL_CELLS_X] =
        calloc(NUM_TOTAL_CELLS_Y, sizeof *density_cell);
    double (*pressure_cell)[NUM_TOTAL_CELLS_X] =
        calloc(NUM_TOTAL_CELLS_Y, sizeof *pressure_cell);
    double (*volume_cell)[NUM_TOTAL_CELLS_X] =
        calloc(NUM_TOTAL_CELLS_Y, sizeof *volume_cell);
    double (*sp_energy_cell)[NUM_TOTAL_CELLS_X] =
        calloc(NUM_TOTAL_CELLS_Y, sizeof *sp_energy_cell);
    double (*sp_momentum_cell)[NUM_TOTAL_CELLS_X][2] =
        calloc(NUM_TOTAL_CELLS_Y, sizeof *sp_momentum_cell);

    int status = EXIT_FAILURE;
    if (!coor_node || !velocity_node || !mass_cell || !density_cell ||
        !pressure_cell || !volume_cell || !sp_energy_cell || !sp_momentum_cell) {
        fprintf(stderr, "Could not allocate the mesh fields.\n");
        goto cleanup;
    }
    printf("Physical cells: x=%zu..%zu, y=%zu..%zu (%zu x %zu)\n",
           START_CELL_X, END_CELL_X, START_CELL_Y, END_CELL_Y,
           NUM_CELLS_X, NUM_CELLS_Y);
    printf("Physical nodes: x=%zu..%zu, y=%zu..%zu (%zu x %zu)\n",
           START_NODE_X, END_NODE_X, START_NODE_Y, END_NODE_Y,
           NUM_NODES_X, NUM_NODES_Y);
    printf("Allocated cells: %zu x %zu; nodes: %zu x %zu\n",
           NUM_TOTAL_CELLS_X, NUM_TOTAL_CELLS_Y,
           NUM_TOTAL_NODES_X, NUM_TOTAL_NODES_Y);

    /******************* Initialization *******************/
    const double cell_h_x = LENGTH_X / (double)NUM_CELLS_X;
    const double cell_h_y = LENGTH_Y / (double)NUM_CELLS_Y;
    for (size_t node_y = START_TOTAL_NODE_Y;
         node_y <= END_TOTAL_NODE_Y; ++node_y) {
        for (size_t node_x = START_TOTAL_NODE_X;
             node_x <= END_TOTAL_NODE_X; ++node_x) {
            initializeNodeCoor(node_x, node_y, START_NODE_X, START_NODE_Y,
                               cell_h_x, cell_h_y, coor_node[node_y][node_x]);
        }
    }
    for (size_t cell_y = START_CELL_Y; cell_y <= END_CELL_Y; ++cell_y) {
        for (size_t cell_x = START_CELL_X; cell_x <= END_CELL_X; ++cell_x) {
            double center_cell[2];
            calCellCenter(coor_node[cell_y][cell_x],
                          coor_node[cell_y][cell_x + 1],
                          coor_node[cell_y + 1][cell_x + 1],
                          coor_node[cell_y + 1][cell_x], center_cell);
            volume_cell[cell_y][cell_x] = initializeCellVolume(cell_h_x, cell_h_y);
            density_cell[cell_y][cell_x] = initializeCellDensity();
            mass_cell[cell_y][cell_x] = initializeCellMass(
                volume_cell[cell_y][cell_x], density_cell[cell_y][cell_x]);
            initializeCellSpMomentum(sp_momentum_cell[cell_y][cell_x]);
            pressure_cell[cell_y][cell_x] = initializeCellPresure(center_cell, GAMMA);
            sp_energy_cell[cell_y][cell_x] = initializeCellSpEnergy(
                density_cell[cell_y][cell_x], sp_momentum_cell[cell_y][cell_x],
                pressure_cell[cell_y][cell_x], GAMMA);
        }
    }
    if (!putOutNodeCoor(0.0, NUM_TOTAL_NODES_X, NUM_TOTAL_NODES_Y,
                       START_NODE_X, END_NODE_X, START_NODE_Y, END_NODE_Y,
                       coor_node))
        goto cleanup;

    /******************* Time advancing *******************/
    double t = 0.0;
    size_t time_step = 0;
    while (t < T) {
        if (!fillGhostNodeCoordinates(
                NUM_TOTAL_NODES_X, NUM_TOTAL_NODES_Y,
                START_NODE_X, END_NODE_X, START_NODE_Y, END_NODE_Y, coor_node)) {
            fprintf(stderr, "Invalid ghost-node geometry at t=%.17g.\n", t);
            goto cleanup;
        }

        /* Ambient states on every ghost side and all four ghost corners. */
        for (size_t cell_y = START_TOTAL_CELL_Y;
             cell_y <= END_TOTAL_CELL_Y; ++cell_y) {
            for (size_t cell_x = START_TOTAL_CELL_X;
                 cell_x <= END_TOTAL_CELL_X; ++cell_x) {
                double area_cell = calCellArea(
                    coor_node[cell_y][cell_x], coor_node[cell_y][cell_x + 1],
                    coor_node[cell_y + 1][cell_x + 1], coor_node[cell_y + 1][cell_x]);
                if (!isfinite(area_cell) || area_cell <= 0.0) {
                    fprintf(stderr, "Invalid cell geometry at (%zu,%zu), t=%.17g.\n",
                            cell_x, cell_y, t);
                    goto cleanup;
                }
                if (cell_x >= START_CELL_X && cell_x <= END_CELL_X &&
                    cell_y >= START_CELL_Y && cell_y <= END_CELL_Y)
                    continue;
                volume_cell[cell_y][cell_x] = area_cell;
                density_cell[cell_y][cell_x] = boundarilizeDensity();
                pressure_cell[cell_y][cell_x] = boundarilizePressure();
                mass_cell[cell_y][cell_x] = initializeCellMass(
                    area_cell, density_cell[cell_y][cell_x]);
                sp_energy_cell[cell_y][cell_x] = boundarilizeSpEnergy(GAMMA);
                boundarilizeSpMomentum(sp_momentum_cell[cell_y][cell_x]);
            }
        }

        /* Keep the original fixed step; CFL control is separate work. */
        const double tau = fmin(1e-3, T - t);

        /* Physical nodes include the physical boundary. Each neighboring
         * cell is addressed directly, in the original accumulation order. */
        for (size_t node_y = START_NODE_Y; node_y <= END_NODE_Y; ++node_y) {
            for (size_t node_x = START_NODE_X; node_x <= END_NODE_X; ++node_x) {
                double pressure_gradient[2] = {0.0, 0.0};
                double ave_velocity_node[2] = {0.0, 0.0};
                double matrix_node[2][2] = {{0.0, 0.0}, {0.0, 0.0}};

                /* Southwest cell: the node is its northeast vertex. */
                int valid = accumulateNodeCell(
                    coor_node[node_y][node_x], coor_node[node_y][node_x - 1],
                    coor_node[node_y - 1][node_x],
                    density_cell[node_y - 1][node_x - 1],
                    pressure_cell[node_y - 1][node_x - 1],
                    sp_momentum_cell[node_y - 1][node_x - 1],
                    GAMMA, matrix_node, pressure_gradient, ave_velocity_node);
                /* Southeast cell: the node is its northwest vertex. */
                valid = valid && accumulateNodeCell(
                    coor_node[node_y][node_x], coor_node[node_y - 1][node_x],
                    coor_node[node_y][node_x + 1],
                    density_cell[node_y - 1][node_x],
                    pressure_cell[node_y - 1][node_x],
                    sp_momentum_cell[node_y - 1][node_x],
                    GAMMA, matrix_node, pressure_gradient, ave_velocity_node);
                /* Northwest cell: the node is its southeast vertex. */
                valid = valid && accumulateNodeCell(
                    coor_node[node_y][node_x], coor_node[node_y + 1][node_x],
                    coor_node[node_y][node_x - 1],
                    density_cell[node_y][node_x - 1],
                    pressure_cell[node_y][node_x - 1],
                    sp_momentum_cell[node_y][node_x - 1],
                    GAMMA, matrix_node, pressure_gradient, ave_velocity_node);
                /* Northeast cell: the node is its southwest vertex. */
                valid = valid && accumulateNodeCell(
                    coor_node[node_y][node_x], coor_node[node_y][node_x + 1],
                    coor_node[node_y + 1][node_x],
                    density_cell[node_y][node_x], pressure_cell[node_y][node_x],
                    sp_momentum_cell[node_y][node_x],
                    GAMMA, matrix_node, pressure_gradient, ave_velocity_node);

                double determinant = matrix_node[0][0] * matrix_node[1][1]
                                   - matrix_node[0][1] * matrix_node[1][0];
                if (!valid || !isfinite(determinant) || determinant <= 0.0) {
                    fprintf(stderr, "Invalid nodal system at (%zu,%zu), t=%.17g.\n",
                            node_x, node_y, t);
                    goto cleanup;
                }
                double node_matrix_inverse[2][2], corrected_force[2];
                inverseMatrix2(matrix_node, node_matrix_inverse);
                vecAddVec2(pressure_gradient, ave_velocity_node, corrected_force);
                matMulVec2(node_matrix_inverse, corrected_force,
                           velocity_node[node_y][node_x]);
                if (!isfinite(velocity_node[node_y][node_x][0]) ||
                    !isfinite(velocity_node[node_y][node_x][1])) {
                    fprintf(stderr, "Nonfinite velocity at node (%zu,%zu).\n",
                            node_x, node_y);
                    goto cleanup;
                }
            }
        }

        for (size_t node_y = START_NODE_Y; node_y <= END_NODE_Y; ++node_y) {
            for (size_t node_x = START_NODE_X; node_x <= END_NODE_X; ++node_x) {
                double position_update[2];
                numMulVec2(tau, velocity_node[node_y][node_x], position_update);
                attVec2(coor_node[node_y][node_x], position_update);
            }
        }
        if (!fillGhostNodeCoordinates(
                NUM_TOTAL_NODES_X, NUM_TOTAL_NODES_Y,
                START_NODE_X, END_NODE_X, START_NODE_Y, END_NODE_Y, coor_node)) {
            fprintf(stderr, "Invalid ghost-node geometry after mesh motion.\n");
            goto cleanup;
        }

        /* Preserve the original update order/formulas. Unresolved volume,
         * normal-sign, and EOS consistency issues are recorded in INDEXING.md. */
        for (size_t cell_y = START_CELL_Y; cell_y <= END_CELL_Y; ++cell_y) {
            for (size_t cell_x = START_CELL_X; cell_x <= END_CELL_X; ++cell_x) {
                double area_cell = calCellArea(
                    coor_node[cell_y][cell_x], coor_node[cell_y][cell_x + 1],
                    coor_node[cell_y + 1][cell_x + 1], coor_node[cell_y + 1][cell_x]);
                if (!isfinite(area_cell) || area_cell <= 0.0) {
                    fprintf(stderr, "Mesh inverted at cell (%zu,%zu), t=%.17g.\n",
                            cell_x, cell_y, t + tau);
                    goto cleanup;
                }
                double volume_flux = 0.0, expan_power = 0.0;
                double resul_force[2] = {0.0, 0.0};
                double impedance = calImpendance(
                    density_cell[cell_y][cell_x], pressure_cell[cell_y][cell_x], GAMMA);

                /* Southwest vertex: next southeast, previous northwest. */
                int valid = accumulateCellCorner(
                    coor_node[cell_y][cell_x], coor_node[cell_y][cell_x + 1],
                    coor_node[cell_y + 1][cell_x], impedance,
                    pressure_cell[cell_y][cell_x], sp_momentum_cell[cell_y][cell_x],
                    velocity_node[cell_y][cell_x], &volume_flux, resul_force, &expan_power);
                /* Southeast vertex: next northeast, previous southwest. */
                valid = valid && accumulateCellCorner(
                    coor_node[cell_y][cell_x + 1], coor_node[cell_y + 1][cell_x + 1],
                    coor_node[cell_y][cell_x], impedance,
                    pressure_cell[cell_y][cell_x], sp_momentum_cell[cell_y][cell_x],
                    velocity_node[cell_y][cell_x + 1], &volume_flux, resul_force, &expan_power);
                /* Northwest vertex: next southwest, previous northeast. */
                valid = valid && accumulateCellCorner(
                    coor_node[cell_y + 1][cell_x], coor_node[cell_y][cell_x],
                    coor_node[cell_y + 1][cell_x + 1], impedance,
                    pressure_cell[cell_y][cell_x], sp_momentum_cell[cell_y][cell_x],
                    velocity_node[cell_y + 1][cell_x], &volume_flux, resul_force, &expan_power);
                /* Northeast vertex: next northwest, previous southeast. */
                valid = valid && accumulateCellCorner(
                    coor_node[cell_y + 1][cell_x + 1], coor_node[cell_y + 1][cell_x],
                    coor_node[cell_y][cell_x + 1], impedance,
                    pressure_cell[cell_y][cell_x], sp_momentum_cell[cell_y][cell_x],
                    velocity_node[cell_y + 1][cell_x + 1], &volume_flux, resul_force, &expan_power);
                if (!valid || !isfinite(volume_flux) || !isfinite(expan_power) ||
                    !isfinite(resul_force[0]) || !isfinite(resul_force[1])) {
                    fprintf(stderr, "Invalid corner update at cell (%zu,%zu).\n",
                            cell_x, cell_y);
                    goto cleanup;
                }
                volume_cell[cell_y][cell_x] +=
                    tau * volume_flux / mass_cell[cell_y][cell_x];
                double contrib_sp_momentum[2];
                numMulVec2(-tau / mass_cell[cell_y][cell_x], resul_force,
                           contrib_sp_momentum);
                attVec2(sp_momentum_cell[cell_y][cell_x], contrib_sp_momentum);
                sp_energy_cell[cell_y][cell_x] -=
                    tau * expan_power / mass_cell[cell_y][cell_x];
            }
        }
        t += tau;
        ++time_step;
        if (time_step % OUTPUT_FREQUENCY == 0 || t >= T) {
            if (!putOutNodeCoor(t, NUM_TOTAL_NODES_X, NUM_TOTAL_NODES_Y,
                               START_NODE_X, END_NODE_X, START_NODE_Y, END_NODE_Y,
                               coor_node))
                goto cleanup;
        }
    }
    printf("Completed %zu time steps to t=%.17g.\n", time_step, t);
    status = EXIT_SUCCESS;

cleanup:
    free(sp_momentum_cell);
    free(sp_energy_cell);
    free(volume_cell);
    free(pressure_cell);
    free(density_cell);
    free(mass_cell);
    free(velocity_node);
    free(coor_node);
    return status;
}

int main(void)
{
    /* Physical counts can be changed independently for a rectangular mesh. */
    const size_t NUM_CELLS_X = 5;
    const size_t NUM_CELLS_Y = 5;
    const double T = 1.0;
    return runSimulation(NUM_CELLS_X, NUM_CELLS_Y, T);
}
