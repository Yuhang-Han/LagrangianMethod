#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

/* Exercise the production solver with independent rectangular dimensions.
 * Compile this file with the normal modules, without a second maire2d.c. */
#define main maire2dExampleMain
#include "../maire2d.c"
#undef main

/* Independent Cartesian first-step reference. It uses the four adjacent
 * cell coordinates directly and does not call production mesh, nodal,
 * pressure, impedance, or matrix helpers. Only physical cell (3, 2) is hot
 * on the 5 x 3 mesh; the initial density is 1 and all velocities are zero.
 * C is the original inward length-times-normal vector. */
static void expectedFirstStepVelocity(
    size_t node_x, size_t node_y,
    size_t NUM_CELLS_X, size_t NUM_CELLS_Y,
    double velocity_expected[2])
{
    const double cell_width = 1.0 / (double)NUM_CELLS_X;
    const double cell_height = 1.0 / (double)NUM_CELLS_Y;
    const double corner_length = 0.5 * hypot(cell_width, cell_height);
    const double hot_pressure = 1e-5 + 0.4 / (acos(-1.0) * 0.02 * 0.02);
    double matrix_xx = 0.0, matrix_xy = 0.0, matrix_yy = 0.0;
    double rhs_x = 0.0, rhs_y = 0.0;

    for (size_t cell_y = node_y - 1; cell_y <= node_y; ++cell_y) {
        for (size_t cell_x = node_x - 1; cell_x <= node_x; ++cell_x) {
            const double pressure =
                (cell_x == 3 && cell_y == 2) ? hot_pressure : 1e-5;
            const double impedance = sqrt(1.4 * pressure);
            const double corner_x =
                (cell_x == node_x ? 0.5 : -0.5) * cell_height;
            const double corner_y =
                (cell_y == node_y ? 0.5 : -0.5) * cell_width;
            const double coefficient = impedance / corner_length;

            matrix_xx += coefficient * corner_x * corner_x;
            matrix_xy += coefficient * corner_x * corner_y;
            matrix_yy += coefficient * corner_y * corner_y;
            rhs_x += pressure * corner_x;
            rhs_y += pressure * corner_y;
        }
    }

    const double determinant = matrix_xx * matrix_yy - matrix_xy * matrix_xy;
    assert(isfinite(determinant) && determinant > 0.0);
    velocity_expected[0] =
        (matrix_yy * rhs_x - matrix_xy * rhs_y) / determinant;
    velocity_expected[1] =
        (matrix_xx * rhs_y - matrix_xy * rhs_x) / determinant;
}

static void checkMeshOutput(size_t NUM_CELLS_X, size_t NUM_CELLS_Y,
                            double time, int use_hot_reference)
{
    const size_t START_NODE_X = 1;
    const size_t END_NODE_X = NUM_CELLS_X + 1;
    const size_t START_NODE_Y = 1;
    const size_t END_NODE_Y = NUM_CELLS_Y + 1;
    char filename[128];
    int length = snprintf(filename, sizeof filename,
                          "build/mesh_t_%.17g.csv", time);
    assert(length > 0 && (size_t)length < sizeof filename);
    FILE *fp = fopen(filename, "r");
    assert(fp != NULL);
    char line[256];
    assert(fgets(line, sizeof line, fp) != NULL);
    assert(strcmp(line, "t,node,i,j,x,y\n") == 0);

    size_t expected_record = 0;
    size_t moving_nodes = 0;
    for (size_t node_y = START_NODE_Y; node_y <= END_NODE_Y; ++node_y) {
        for (size_t node_x = START_NODE_X; node_x <= END_NODE_X; ++node_x) {
            double output_time, x, y;
            size_t node_record, output_x, output_y;
            char newline;
            assert(fgets(line, sizeof line, fp) != NULL);
            assert(sscanf(line, "%lf,%zu,%zu,%zu,%lf,%lf%c",
                          &output_time, &node_record, &output_x, &output_y,
                          &x, &y, &newline) == 7);
            assert(newline == '\n');
            assert(isfinite(output_time) && isfinite(x) && isfinite(y));
            assert(output_time == time);
            assert(node_record == expected_record);
            assert(output_x == node_x && output_y == node_y);

            const double initial_x =
                (double)(node_x - START_NODE_X) / (double)NUM_CELLS_X;
            const double initial_y =
                (double)(node_y - START_NODE_Y) / (double)NUM_CELLS_Y;
            double velocity_expected[2] = {0.0, 0.0};
            if (use_hot_reference) {
                expectedFirstStepVelocity(node_x, node_y,
                                           NUM_CELLS_X, NUM_CELLS_Y,
                                           velocity_expected);
                if (hypot(velocity_expected[0], velocity_expected[1]) > 1.0) {
                    ++moving_nodes;
                    assert(hypot(x - initial_x, y - initial_y) > 1e-4);
                }
            }
            const double expected_x = initial_x + time * velocity_expected[0];
            const double expected_y = initial_y + time * velocity_expected[1];
            /* Includes rounding in the ill-conditioned hot/cold nodal solve;
             * small enough to detect a swapped neighbor or normal component. */
            assert(fabs(x - expected_x) <= 5e-13 * (1.0 + fabs(expected_x)));
            assert(fabs(y - expected_y) <= 5e-13 * (1.0 + fabs(expected_y)));
            ++expected_record;
        }
        if (node_y < END_NODE_Y) {
            assert(fgets(line, sizeof line, fp) != NULL);
            assert(strcmp(line, "\n") == 0);
        }
    }
    assert(expected_record == (NUM_CELLS_X + 1) * (NUM_CELLS_Y + 1));
    if (use_hot_reference)
        assert(moving_nodes == 4);
    assert(fgets(line, sizeof line, fp) == NULL);
    assert(feof(fp));
    assert(fclose(fp) == 0);
    assert(remove(filename) == 0);
}

static void testStationaryRectangularMesh(void)
{
    /* No 3 x 2 cell center falls in the explosion radius. Two steps also
     * exercise refreshed ghosts and cell updates in a uniform state. */
    assert(runSimulation(3, 2, 0.002) == EXIT_SUCCESS);
    checkMeshOutput(3, 2, 0.0, 0);
    checkMeshOutput(3, 2, 0.002, 0);
}

static void testActiveRectangularFirstStep(void)
{
    /* Unequal x/y spacing makes swapped dimensions observable. The hot
     * center cell produces nonzero motion in all four corner orientations.
     * This tests the existing discrete formula, not full Sedov correctness. */
    assert(runSimulation(5, 3, 0.001) == EXIT_SUCCESS);
    checkMeshOutput(5, 3, 0.0, 0);
    checkMeshOutput(5, 3, 0.001, 1);
}

static void testInvalidSimulationInputs(void)
{
    assert(runSimulation(0, 2, 0.0) == EXIT_FAILURE);
    assert(runSimulation(3, 0, 0.0) == EXIT_FAILURE);
    assert(runSimulation(SIZE_MAX, 2, 0.0) == EXIT_FAILURE);
    /* Valid additions can still overflow the contiguous allocation size. */
    assert(runSimulation(SIZE_MAX / (2 * sizeof(double)), 1, 0.0) == EXIT_FAILURE);
    assert(runSimulation(1, SIZE_MAX / (8 * sizeof(double)), 0.0) == EXIT_FAILURE);
    assert(runSimulation(3, 2, NAN) == EXIT_FAILURE);
    assert(runSimulation(3, 2, -0.001) == EXIT_FAILURE);
}

int main(void)
{
    testStationaryRectangularMesh();
    testActiveRectangularFirstStep();
    testInvalidSimulationInputs();
    puts("solver_test: PASS");
    return 0;
}
