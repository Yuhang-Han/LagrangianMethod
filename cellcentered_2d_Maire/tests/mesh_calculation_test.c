#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>

#include "../mesh_calculation.h"

static void checkNeighbor(size_t node_x, size_t node_y,
                          size_t cell_x, size_t cell_y, int direction,
                          size_t expected_node_x, size_t expected_node_y)
{
    size_t neigh_node_x, neigh_node_y;
    calNeighVertexIdx(node_x, node_y, cell_x, cell_y, direction,
                     &neigh_node_x, &neigh_node_y);
    assert(neigh_node_x == expected_node_x);
    assert(neigh_node_y == expected_node_y);
}

static void testEveryCornerAndDirection(void)
{
    /* Cell (2,1): SW=(2,1), SE=(3,1), NE=(3,2), NW=(2,2). */
    checkNeighbor(2, 1, 2, 1,  1, 3, 1);
    checkNeighbor(2, 1, 2, 1, -1, 2, 2);
    checkNeighbor(3, 1, 2, 1,  1, 3, 2);
    checkNeighbor(3, 1, 2, 1, -1, 2, 1);
    checkNeighbor(3, 2, 2, 1,  1, 2, 2);
    checkNeighbor(3, 2, 2, 1, -1, 3, 1);
    checkNeighbor(2, 2, 2, 1,  1, 2, 1);
    checkNeighbor(2, 2, 2, 1, -1, 3, 2);
}

static void testRectangularMeshAndGhostAdjacentNodes(void)
{
    const size_t NUM_CELLS_X = 3;
    const size_t NUM_CELLS_Y = 2;
    const size_t NUM_NODES_X = NUM_CELLS_X + 1;
    const size_t NUM_NODES_Y = NUM_CELLS_Y + 1;
    const size_t NUM_GHOST_LAYERS = 1;
    const size_t NUM_TOTAL_CELLS_X = NUM_CELLS_X + 2 * NUM_GHOST_LAYERS;
    const size_t NUM_TOTAL_CELLS_Y = NUM_CELLS_Y + 2 * NUM_GHOST_LAYERS;
    const size_t NUM_TOTAL_NODES_X = NUM_TOTAL_CELLS_X + 1;
    const size_t NUM_TOTAL_NODES_Y = NUM_TOTAL_CELLS_Y + 1;
    const size_t START_NODE_X = NUM_GHOST_LAYERS;
    const size_t END_NODE_X = START_NODE_X + NUM_NODES_X - 1;
    const size_t START_NODE_Y = NUM_GHOST_LAYERS;
    const size_t END_NODE_Y = START_NODE_Y + NUM_NODES_Y - 1;

    /* Include physical boundary nodes and their neighboring ghost cells. */
    for (size_t node_y = START_NODE_Y; node_y <= END_NODE_Y; ++node_y) {
        for (size_t node_x = START_NODE_X; node_x <= END_NODE_X; ++node_x) {
            assert(node_x > 0 && node_y > 0);
            assert(node_x < NUM_TOTAL_CELLS_X);
            assert(node_y < NUM_TOTAL_CELLS_Y);
            assert(node_x + 1 < NUM_TOTAL_NODES_X);
            assert(node_y + 1 < NUM_TOTAL_NODES_Y);

            /* The node is the NE vertex of the cell to its southwest. */
            checkNeighbor(node_x, node_y, node_x - 1, node_y - 1,
                           1, node_x - 1, node_y);
            checkNeighbor(node_x, node_y, node_x - 1, node_y - 1,
                          -1, node_x, node_y - 1);

            /* The node is the NW vertex of the cell to its southeast. */
            checkNeighbor(node_x, node_y, node_x, node_y - 1,
                           1, node_x, node_y - 1);
            checkNeighbor(node_x, node_y, node_x, node_y - 1,
                          -1, node_x + 1, node_y);

            /* The node is the SW vertex of the cell to its northeast. */
            checkNeighbor(node_x, node_y, node_x, node_y,
                           1, node_x + 1, node_y);
            checkNeighbor(node_x, node_y, node_x, node_y,
                          -1, node_x, node_y + 1);

            /* The node is the SE vertex of the cell to its northwest. */
            checkNeighbor(node_x, node_y, node_x - 1, node_y,
                           1, node_x, node_y + 1);
            checkNeighbor(node_x, node_y, node_x - 1, node_y,
                          -1, node_x - 1, node_y);
        }
    }
}

static void expectNear(double actual, double expected)
{
    assert(isfinite(actual));
    assert(fabs(actual - expected) <= 2e-12 * fmax(1.0, fabs(expected)));
}

static void testAnalyticCellAreas(void)
{
    const struct {
        double vertices[4][2];
        double area;
    } cases[] = {
        {{{0.0, 0.0}, {2.0, 0.0}, {2.0, 3.0}, {0.0, 3.0}}, 6.0},
        /* Parallelogram with spanning vectors (3,1) and (2,3). */
        {{{0.0, 0.0}, {3.0, 1.0}, {5.0, 4.0}, {2.0, 3.0}}, 7.0},
        /* A 4x3 rectangle minus triangles of areas 1.5 and 2.0. */
        {{{0.0, 0.0}, {4.0, 0.0}, {3.0, 2.0}, {0.0, 3.0}}, 8.5},
        {{{0.0, 0.0}, {0.5, 0.0}, {0.5, 0.25}, {0.0, 0.25}}, 0.125}
    };
    const double offsets[][2] = {
        {0.0, 0.0},
        {7.25, -13.5},
        /* Dyadic dimensions remain exactly representable at this offset. */
        {1099511627776.0, -549755813888.0}
    };

    for (size_t c = 0; c < sizeof cases / sizeof cases[0]; ++c) {
        for (size_t shift = 0; shift < sizeof offsets / sizeof offsets[0]; ++shift) {
            double q[4][2];
            for (size_t vertex = 0; vertex < 4; ++vertex)
                for (size_t component = 0; component < 2; ++component)
                    q[vertex][component] = cases[c].vertices[vertex][component]
                                           + offsets[shift][component];
            expectNear(calCellVolume(q[0], q[1], q[2], q[3]), cases[c].area);
            /* Clockwise orientation must stay negative, not become fabs(area). */
            expectNear(calCellVolume(q[0], q[3], q[2], q[1]), -cases[c].area);
        }
    }

    const double line[4][2] = {{0.0, 0.0}, {1.0, 1.0}, {3.0, 3.0}, {2.0, 2.0}};
    expectNear(calCellVolume(line[0], line[1], line[2], line[3]), 0.0);
    expectNear(calCellVolume(line[0], line[0], line[0], line[0]), 0.0);
}

static double cornerVolumeFlux(double coor[4][2], double velocity[4][2])
{
    double flux = 0.0;
    for (size_t vertex = 0; vertex < 4; ++vertex) {
        double normal[2];
        const size_t next = (vertex + 1) % 4;
        const size_t prev = (vertex + 3) % 4;
        const double length = calCornerLength(coor[vertex], coor[next], coor[prev]);
        calCornerNormal(coor[vertex], coor[next], coor[prev], normal);
        flux += length * (normal[0] * velocity[vertex][0]
                          + normal[1] * velocity[vertex][1]);
    }
    return flux;
}

static void moveVertices(double coor[4][2], double velocity[4][2], double dt)
{
    for (size_t vertex = 0; vertex < 4; ++vertex)
        for (size_t component = 0; component < 2; ++component)
            coor[vertex][component] += dt * velocity[vertex][component];
}

static void testExpandingCell(void)
{
    double q[4][2] = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};
    double velocity[4][2] = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};
    const double initial_volume = calCellVolume(q[0], q[1], q[2], q[3]);
    const double initial_flux = cornerVolumeFlux(q, velocity);
    moveVertices(q, velocity, 0.1);
    const double final_flux = cornerVolumeFlux(q, velocity);
    const double final_volume = calCellVolume(q[0], q[1], q[2], q[3]);

    /* x_new=1.1*x_old: both lengths grow by 1.1, so area grows by 1.21. */
    expectNear(initial_volume, 1.0);
    expectNear(initial_flux, 2.0);
    expectNear(final_flux, 2.2);
    expectNear(final_volume, 1.21);
    expectNear(initial_volume + 0.1 * final_flux, 1.22);
    expectNear(final_volume - initial_volume, 0.05 * (initial_flux + final_flux));
}

static void testRepeatedAffineMotion(void)
{
    const struct {
        double gradient[2][2];
        double determinant;
    } motions[] = {
        /* Pure shear preserves area. */
        {{{0.0, 0.7}, {0.0, 0.0}}, 1.0},
        /* det(I+0.2*A)=det([[1.06,-0.08],[0.04,0.98]])=1.042. */
        {{{0.3, -0.4}, {0.2, -0.1}}, 1.042}
    };

    for (size_t motion = 0; motion < sizeof motions / sizeof motions[0]; ++motion) {
        double q[4][2] = {{0.0, 0.0}, {2.0, 0.0}, {2.0, 3.0}, {0.0, 3.0}};
        for (size_t step = 0; step < 20; ++step) {
            double velocity[4][2];
            for (size_t vertex = 0; vertex < 4; ++vertex) {
                velocity[vertex][0] = motions[motion].gradient[0][0] * q[vertex][0]
                                      + motions[motion].gradient[0][1] * q[vertex][1] + 0.7;
                velocity[vertex][1] = motions[motion].gradient[1][0] * q[vertex][0]
                                      + motions[motion].gradient[1][1] * q[vertex][1] - 0.2;
            }
            moveVertices(q, velocity, 0.2);
            /* Translation has no effect; the determinant controls each step. */
            expectNear(calCellVolume(q[0], q[1], q[2], q[3]),
                       6.0 * pow(motions[motion].determinant, (double)(step + 1)));
        }
    }
}

static void testNonAffineGeometricConservation(void)
{
    double q[4][2] = {{0.0, 0.0}, {2.0, 0.0}, {2.5, 1.5}, {0.0, 1.0}};
    /* These velocities cannot all be represented by one affine velocity field. */
    double velocity[4][2] = {{0.1, -0.2}, {0.3, 0.4}, {-0.2, 0.5}, {0.4, -0.1}};
    const double initial_volume = calCellVolume(q[0], q[1], q[2], q[3]);
    const double initial_flux = cornerVolumeFlux(q, velocity);
    moveVertices(q, velocity, 0.2);
    const double final_volume = calCellVolume(q[0], q[1], q[2], q[3]);
    const double final_flux = cornerVolumeFlux(q, velocity);

    expectNear(initial_volume, 2.75);
    expectNear(final_volume, 2.7216);
    /* For linearly moving vertices, corner volume flux is linear in time. */
    expectNear(final_volume - initial_volume, 0.1 * (initial_flux + final_flux));
    assert(fabs(final_volume - (initial_volume + 0.2 * final_flux)) > 1e-4);
}

int main(void)
{
    testEveryCornerAndDirection();
    testRectangularMeshAndGhostAdjacentNodes();
    testAnalyticCellAreas();
    testExpandingCell();
    testRepeatedAffineMotion();
    testNonAffineGeometricConservation();
    puts("mesh_calculation_test: PASS (neighbors, signed areas, and geometric conservation)");
    return 0;
}
