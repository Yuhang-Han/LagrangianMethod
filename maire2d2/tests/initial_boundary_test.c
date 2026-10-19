#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "../initial_bound_condition.h"

int main(void)
{
    /* A rectangular padded mesh: 6 x 5 nodes, with one ghost layer. */
    const double cell_h_0 = 0.2;
    const double cell_h_1 = 0.3;
    const struct {
        size_t node_x, node_y;
        double x, y;
    } samples[] = {
        {0, 0, 0.0, 0.0},
        {1, 1, 0.2, 0.3},
        {4, 3, 0.8, 0.9},
        {5, 4, 1.0, 1.2},
        {0, 2, 0.0, 0.6},
        {5, 1, 1.0, 0.3},
        {2, 0, 0.4, 0.0},
        {3, 4, 0.6, 1.2}
    };

    for (size_t sample = 0; sample < sizeof samples / sizeof samples[0]; ++sample) {
        double coor_node[2] = {-99.0, -99.0};
        initializeNodeCoor(samples[sample].node_x, samples[sample].node_y,
                           cell_h_0, cell_h_1, coor_node);
        assert(fabs(coor_node[0] - samples[sample].x) < 1e-12);
        assert(fabs(coor_node[1] - samples[sample].y) < 1e-12);

        coor_node[0] = 99.0;
        coor_node[1] = -99.0;
        boundarilizeCoor(samples[sample].node_x, samples[sample].node_y,
                        cell_h_0, cell_h_1, coor_node);
        assert(fabs(coor_node[0] - samples[sample].x) < 1e-12);
        assert(fabs(coor_node[1] - samples[sample].y) < 1e-12);
    }

    puts("initial_boundary_test: PASS");
    return 0;
}
