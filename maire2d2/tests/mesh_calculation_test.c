#include <assert.h>
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

int main(void)
{
    testEveryCornerAndDirection();
    testRectangularMeshAndGhostAdjacentNodes();
    puts("mesh_calculation_test: neighbor-index assertions passed");
    return 0;
}
