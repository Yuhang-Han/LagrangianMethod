#include <assert.h>
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>

#include "../mesh_calculation.h"

static void assertClose(double actual, double expected)
{
    const double tolerance = 64.0 * DBL_EPSILON * fmax(1.0, fabs(expected));
    assert(isfinite(actual));
    assert(fabs(actual - expected) <= tolerance);
}

/* Independent reference: the original two normalized half-edge normals. */
static void originalCornerGeometry(const double node[2],
                                   const double next[2],
                                   const double prev[2],
                                   double *length, double normal[2])
{
    const double edge_next[2] = {node[0] - next[0], node[1] - next[1]};
    const double edge_prev[2] = {prev[0] - node[0], prev[1] - node[1]};
    const double length_next = sqrt(edge_next[0] * edge_next[0] +
                                    edge_next[1] * edge_next[1]);
    const double length_prev = sqrt(edge_prev[0] * edge_prev[0] +
                                    edge_prev[1] * edge_prev[1]);
    const double normal_next[2] = {
        edge_next[1] / length_next, -edge_next[0] / length_next
    };
    const double normal_prev[2] = {
        edge_prev[1] / length_prev, -edge_prev[0] / length_prev
    };
    const double directed[2] = {
        0.5 * length_next * normal_next[0] +
            0.5 * length_prev * normal_prev[0],
        0.5 * length_next * normal_next[1] +
            0.5 * length_prev * normal_prev[1]
    };

    *length = sqrt(directed[0] * directed[0] + directed[1] * directed[1]);
    normal[0] = directed[0] / *length;
    normal[1] = directed[1] / *length;
}

static void assertOriginalCorner(const double node[2], const double next[2],
                                 const double prev[2], double directed[2])
{
    double length, normal[2], original_length, original_normal[2];
    assert(calCornerGeometry(node, next, prev, &length, normal));
    originalCornerGeometry(node, next, prev, &original_length, original_normal);
    assertClose(length, original_length);
    assertClose(normal[0], original_normal[0]);
    assertClose(normal[1], original_normal[1]);
    assertClose(hypot(normal[0], normal[1]), 1.0);
    directed[0] = length * normal[0];
    directed[1] = length * normal[1];
}

static void testRectangularGrid(void)
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
    const size_t START_CELL_X = NUM_GHOST_LAYERS;
    const size_t END_CELL_X = START_CELL_X + NUM_CELLS_X - 1;
    const size_t START_CELL_Y = NUM_GHOST_LAYERS;
    const size_t END_CELL_Y = START_CELL_Y + NUM_CELLS_Y - 1;
    const size_t START_NODE_X = NUM_GHOST_LAYERS;
    const size_t END_NODE_X = START_NODE_X + NUM_NODES_X - 1;
    const size_t START_NODE_Y = NUM_GHOST_LAYERS;
    const size_t END_NODE_Y = START_NODE_Y + NUM_NODES_Y - 1;
    const double DX = 0.5;
    const double DY = 1.25;
    double coor_node[NUM_TOTAL_NODES_Y][NUM_TOTAL_NODES_X][2];
    double center_cell[NUM_TOTAL_CELLS_Y][NUM_TOTAL_CELLS_X][2];
    double area_cell[NUM_TOTAL_CELLS_Y][NUM_TOTAL_CELLS_X];

    for (size_t node_y = 0; node_y < NUM_TOTAL_NODES_Y; ++node_y) {
        for (size_t node_x = 0; node_x < NUM_TOTAL_NODES_X; ++node_x) {
            coor_node[node_y][node_x][0] =
                ((double)node_x - (double)START_NODE_X) * DX;
            coor_node[node_y][node_x][1] =
                ((double)node_y - (double)START_NODE_Y) * DY;
        }
    }

    /* Compute both physical and ghost cells directly from their four nodes. */
    for (size_t cell_y = 0; cell_y < NUM_TOTAL_CELLS_Y; ++cell_y) {
        for (size_t cell_x = 0; cell_x < NUM_TOTAL_CELLS_X; ++cell_x) {
            assert(cell_x + 1 < NUM_TOTAL_NODES_X);
            assert(cell_y + 1 < NUM_TOTAL_NODES_Y);
            area_cell[cell_y][cell_x] = calCellArea(
                coor_node[cell_y][cell_x],
                coor_node[cell_y][cell_x + 1],
                coor_node[cell_y + 1][cell_x + 1],
                coor_node[cell_y + 1][cell_x]);
            calCellCenter(coor_node[cell_y][cell_x],
                          coor_node[cell_y][cell_x + 1],
                          coor_node[cell_y + 1][cell_x + 1],
                          coor_node[cell_y + 1][cell_x],
                          center_cell[cell_y][cell_x]);
            assertClose(area_cell[cell_y][cell_x], DX * DY);
            assertClose(center_cell[cell_y][cell_x][0],
                        coor_node[cell_y][cell_x][0] + 0.5 * DX);
            assertClose(center_cell[cell_y][cell_x][1],
                        coor_node[cell_y][cell_x][1] + 0.5 * DY);
        }
    }

    size_t physical_cells = 0;
    double physical_area = 0.0;
    for (size_t cell_y = START_CELL_Y; cell_y <= END_CELL_Y; ++cell_y) {
        for (size_t cell_x = START_CELL_X; cell_x <= END_CELL_X; ++cell_x) {
            ++physical_cells;
            physical_area += area_cell[cell_y][cell_x];
        }
    }
    assert(physical_cells == NUM_CELLS_X * NUM_CELLS_Y);
    assertClose(physical_area, (double)physical_cells * DX * DY);

    size_t physical_nodes = 0;
    for (size_t node_y = START_NODE_Y; node_y <= END_NODE_Y; ++node_y) {
        for (size_t node_x = START_NODE_X; node_x <= END_NODE_X; ++node_x) {
            ++physical_nodes;
            assert(node_x > 0 && node_x < NUM_TOTAL_CELLS_X);
            assert(node_y > 0 && node_y < NUM_TOTAL_CELLS_Y);

            /* All four surrounding cells exist, including at boundary nodes. */
            assertClose(center_cell[node_y - 1][node_x - 1][0],
                        coor_node[node_y][node_x][0] - 0.5 * DX);
            assertClose(center_cell[node_y - 1][node_x - 1][1],
                        coor_node[node_y][node_x][1] - 0.5 * DY);
            assertClose(center_cell[node_y - 1][node_x][0],
                        coor_node[node_y][node_x][0] + 0.5 * DX);
            assertClose(center_cell[node_y - 1][node_x][1],
                        coor_node[node_y][node_x][1] - 0.5 * DY);
            assertClose(center_cell[node_y][node_x][0],
                        coor_node[node_y][node_x][0] + 0.5 * DX);
            assertClose(center_cell[node_y][node_x][1],
                        coor_node[node_y][node_x][1] + 0.5 * DY);
            assertClose(center_cell[node_y][node_x - 1][0],
                        coor_node[node_y][node_x][0] - 0.5 * DX);
            assertClose(center_cell[node_y][node_x - 1][1],
                        coor_node[node_y][node_x][1] + 0.5 * DY);
        }
    }
    assert(physical_nodes == NUM_NODES_X * NUM_NODES_Y);
    assertClose(coor_node[START_NODE_Y][START_NODE_X][0], 0.0);
    assertClose(coor_node[START_NODE_Y][START_NODE_X][1], 0.0);
    assertClose(coor_node[END_NODE_Y][END_NODE_X][0], (double)NUM_CELLS_X * DX);
    assertClose(coor_node[END_NODE_Y][END_NODE_X][1], (double)NUM_CELLS_Y * DY);
}

static void testInwardNormal(void)
{
    const double sw[2] = {0.0, 0.0};
    const double se[2] = {1.0, 0.0};
    const double ne[2] = {1.0, 1.0};
    const double nw[2] = {0.0, 1.0};
    double directed_sw[2], directed_se[2], directed_ne[2], directed_nw[2];

    assertOriginalCorner(sw, se, nw, directed_sw);
    assertOriginalCorner(se, ne, sw, directed_se);
    assertOriginalCorner(ne, nw, se, directed_ne);
    assertOriginalCorner(nw, sw, ne, directed_nw);
    assertClose(directed_sw[0], 0.5);
    assertClose(directed_sw[1], 0.5);
    assertClose(directed_se[0], -0.5);
    assertClose(directed_se[1], 0.5);
    assertClose(directed_ne[0], -0.5);
    assertClose(directed_ne[1], -0.5);
    assertClose(directed_nw[0], 0.5);
    assertClose(directed_nw[1], -0.5);
}

static void testSkewedGeometryAndTranslation(void)
{
    const double sw[2] = {0.25, -0.5};
    const double se[2] = {2.5, 0.0};
    const double ne[2] = {3.0, 1.75};
    const double nw[2] = {-0.25, 1.25};
    const double translation[2] = {10.0, -7.0};
    const double uniform_velocity[2] = {0.3, -0.7};
    double moved_sw[2], moved_se[2], moved_ne[2], moved_nw[2];
    double directed_sw[2], directed_se[2], directed_ne[2], directed_nw[2];
    double moved_directed_sw[2], moved_directed_se[2];
    double moved_directed_ne[2], moved_directed_nw[2];
    double center[2], moved_center[2];

    assertOriginalCorner(sw, se, nw, directed_sw);
    assertOriginalCorner(se, ne, sw, directed_se);
    assertOriginalCorner(ne, nw, se, directed_ne);
    assertOriginalCorner(nw, sw, ne, directed_nw);
    assertClose(calCellArea(sw, se, ne, nw), 4.8125);
    assertClose(calCellArea(sw, nw, ne, se), -4.8125);
    calCellCenter(sw, se, ne, nw, center);
    assertClose(center[0], 1.375);
    assertClose(center[1], 0.625);

    double uniform_area_rate = 0.0;
    for (size_t component = 0; component < 2; ++component) {
        const double closed_corner_sum = directed_sw[component] +
            directed_se[component] + directed_ne[component] + directed_nw[component];
        assertClose(closed_corner_sum, 0.0);
        uniform_area_rate += uniform_velocity[component] * closed_corner_sum;
        moved_sw[component] = sw[component] + translation[component];
        moved_se[component] = se[component] + translation[component];
        moved_ne[component] = ne[component] + translation[component];
        moved_nw[component] = nw[component] + translation[component];
    }
    assertClose(uniform_area_rate, 0.0);
    assertClose(calCellArea(moved_sw, moved_se, moved_ne, moved_nw), 4.8125);
    calCellCenter(moved_sw, moved_se, moved_ne, moved_nw, moved_center);
    assertOriginalCorner(moved_sw, moved_se, moved_nw, moved_directed_sw);
    assertOriginalCorner(moved_se, moved_ne, moved_sw, moved_directed_se);
    assertOriginalCorner(moved_ne, moved_nw, moved_se, moved_directed_ne);
    assertOriginalCorner(moved_nw, moved_sw, moved_ne, moved_directed_nw);
    for (size_t component = 0; component < 2; ++component) {
        assertClose(moved_center[component], center[component] + translation[component]);
        assertClose(moved_directed_sw[component], directed_sw[component]);
        assertClose(moved_directed_se[component], directed_se[component]);
        assertClose(moved_directed_ne[component], directed_ne[component]);
        assertClose(moved_directed_nw[component], directed_nw[component]);
    }
}

static void testDegenerateGeometry(void)
{
    const double node[2] = {0.0, 0.0};
    const double next[2] = {1.0, 0.0};
    const double prev[2] = {0.0, 1.0};
    const double nan_node[2] = {NAN, 0.0};
    const double infinite_node[2] = {0.0, INFINITY};
    const double line_ne[2] = {2.0, 0.0};
    const double line_nw[2] = {3.0, 0.0};
    double length = -1.0;
    double normal[2] = {-1.0, -1.0};

    assert(!calCornerGeometry(node, node, prev, &length, normal));
    assertClose(length, 0.0);
    assertClose(normal[0], 0.0);
    assertClose(normal[1], 0.0);
    assert(!calCornerGeometry(node, next, node, &length, normal));
    assert(!calCornerGeometry(node, next, next, &length, normal));
    assert(!calCornerGeometry(nan_node, next, prev, &length, normal));
    assert(!calCornerGeometry(node, infinite_node, prev, &length, normal));
    assert(!calCornerGeometry(NULL, next, prev, &length, normal));
    assertClose(length, 0.0);
    assertClose(normal[0], 0.0);
    assertClose(normal[1], 0.0);
    assertClose(calCellArea(node, next, line_ne, line_nw), 0.0);
}

int main(void)
{
    testRectangularGrid();
    testInwardNormal();
    testSkewedGeometryAndTranslation();
    testDegenerateGeometry();
    puts("mesh_calculation_test: all assertions passed");
    return 0;
}
