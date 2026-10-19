#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "../initial_bound_condition.h"

static void assertNear(double actual, double expected)
{
    assert(isfinite(actual));
    assert(fabs(actual - expected) < 1e-12);
}

static void testCoordinateOrigin(void)
{
    double coor_node[2];
    initializeNodeCoor(2, 3, 2, 3, 0.2, 0.25, coor_node);
    assertNear(coor_node[0], 0.0);
    assertNear(coor_node[1], 0.0);

    initializeNodeCoor(5, 5, 2, 3, 0.2, 0.25, coor_node);
    assertNear(coor_node[0], 0.6);
    assertNear(coor_node[1], 0.5);

    initializeNodeCoor(0, 1, 2, 3, 0.2, 0.25, coor_node);
    assertNear(coor_node[0], -0.4);
    assertNear(coor_node[1], -0.5);
}

/* A deformed bilinear mesh has an independently known linear extension. */
static void deformedPosition(double x, double y, double position[2])
{
    position[0] = 1.25 + 0.5 * x + 0.125 * y + 0.0625 * x * y;
    position[1] = -0.75 + 0.25 * x + 0.75 * y - 0.03125 * x * y;
}

static void testGhostExtension(size_t NUM_CELLS_X, size_t NUM_CELLS_Y,
                               size_t NUM_GHOST_LAYERS)
{
    const size_t NUM_TOTAL_NODES_X = NUM_CELLS_X + 1 + 2 * NUM_GHOST_LAYERS;
    const size_t NUM_TOTAL_NODES_Y = NUM_CELLS_Y + 1 + 2 * NUM_GHOST_LAYERS;
    const size_t START_NODE_X = NUM_GHOST_LAYERS;
    const size_t START_NODE_Y = NUM_GHOST_LAYERS;
    const size_t END_NODE_X = START_NODE_X + NUM_CELLS_X;
    const size_t END_NODE_Y = START_NODE_Y + NUM_CELLS_Y;
    double coor_node[NUM_TOTAL_NODES_Y][NUM_TOTAL_NODES_X][2];
    double original_node[NUM_TOTAL_NODES_Y][NUM_TOTAL_NODES_X][2];

    for (size_t node_y = 0; node_y < NUM_TOTAL_NODES_Y; ++node_y) {
        for (size_t node_x = 0; node_x < NUM_TOTAL_NODES_X; ++node_x) {
            coor_node[node_y][node_x][0] = NAN;
            coor_node[node_y][node_x][1] = NAN;
        }
    }
    for (size_t node_y = START_NODE_Y; node_y <= END_NODE_Y; ++node_y) {
        for (size_t node_x = START_NODE_X; node_x <= END_NODE_X; ++node_x) {
            deformedPosition((double)node_x - (double)START_NODE_X,
                             (double)node_y - (double)START_NODE_Y,
                             coor_node[node_y][node_x]);
        }
    }
    memcpy(original_node, coor_node, sizeof coor_node);

    assert(fillGhostNodeCoordinates(NUM_TOTAL_NODES_X, NUM_TOTAL_NODES_Y,
                                    START_NODE_X, END_NODE_X,
                                    START_NODE_Y, END_NODE_Y, coor_node));

    /* Check every ghost, including all four corners and the outer layer. */
    for (size_t node_y = 0; node_y < NUM_TOTAL_NODES_Y; ++node_y) {
        for (size_t node_x = 0; node_x < NUM_TOTAL_NODES_X; ++node_x) {
            double expected[2];
            deformedPosition((double)node_x - (double)START_NODE_X,
                             (double)node_y - (double)START_NODE_Y, expected);
            assertNear(coor_node[node_y][node_x][0], expected[0]);
            assertNear(coor_node[node_y][node_x][1], expected[1]);
            if (node_x >= START_NODE_X && node_x <= END_NODE_X &&
                node_y >= START_NODE_Y && node_y <= END_NODE_Y) {
                assert(memcmp(coor_node[node_y][node_x],
                              original_node[node_y][node_x],
                              sizeof coor_node[node_y][node_x]) == 0);
            }
        }
    }

    /* A later boundary update must follow the moved physical mesh. */
    for (size_t node_y = START_NODE_Y; node_y <= END_NODE_Y; ++node_y) {
        for (size_t node_x = START_NODE_X; node_x <= END_NODE_X; ++node_x) {
            coor_node[node_y][node_x][0] += 0.375;
            coor_node[node_y][node_x][1] -= 0.25;
        }
    }
    assert(fillGhostNodeCoordinates(NUM_TOTAL_NODES_X, NUM_TOTAL_NODES_Y,
                                    START_NODE_X, END_NODE_X,
                                    START_NODE_Y, END_NODE_Y, coor_node));
    for (size_t node_y = 0; node_y < NUM_TOTAL_NODES_Y; ++node_y) {
        for (size_t node_x = 0; node_x < NUM_TOTAL_NODES_X; ++node_x) {
            double expected[2];
            deformedPosition((double)node_x - (double)START_NODE_X,
                             (double)node_y - (double)START_NODE_Y, expected);
            assertNear(coor_node[node_y][node_x][0], expected[0] + 0.375);
            assertNear(coor_node[node_y][node_x][1], expected[1] - 0.25);
        }
    }
}

static void testInvalidGhostInputs(void)
{
    double coor_node[4][4][2] = {{{0}}};
    double original_node[4][4][2];
    memcpy(original_node, coor_node, sizeof coor_node);
    assert(!fillGhostNodeCoordinates(4, 4, 2, 1, 1, 2, coor_node));
    assert(!fillGhostNodeCoordinates(4, 4, 1, 4, 1, 2, coor_node));
    assert(!fillGhostNodeCoordinates(4, 4, 1, 2, 1, 4, coor_node));
    assert(!fillGhostNodeCoordinates(4, 4, 1, 1, 1, 2, coor_node));
    assert(memcmp(original_node, coor_node, sizeof coor_node) == 0);

    coor_node[1][2][1] = NAN;
    memcpy(original_node, coor_node, sizeof coor_node);
    assert(!fillGhostNodeCoordinates(4, 4, 1, 2, 1, 2, coor_node));
    assert(memcmp(original_node, coor_node, sizeof coor_node) == 0);
}

int main(void)
{
    testCoordinateOrigin();
    testGhostExtension(3, 2, 1);
    testGhostExtension(3, 2, 2);
    testGhostExtension(1, 1, 1);
    testGhostExtension(1, 1, 2);
    testInvalidGhostInputs();
    puts("initial_boundary_test: PASS");
    return 0;
}
