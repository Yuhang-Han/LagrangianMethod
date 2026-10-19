#include <math.h>
#include <stddef.h>

#include "mesh_calculation.h"

int calCornerGeometry(const double coor_node[2],
                      const double coor_next[2],
                      const double coor_prev[2],
                      double *length_corner,
                      double normal_corner[2])
{
    if (length_corner == NULL || normal_corner == NULL) {
        return 0;
    }

    *length_corner = 0.0;
    normal_corner[0] = 0.0;
    normal_corner[1] = 0.0;

    if (coor_node == NULL || coor_next == NULL || coor_prev == NULL) {
        return 0;
    }

    for (size_t component = 0; component < 2; ++component) {
        if (!isfinite(coor_node[component]) ||
            !isfinite(coor_next[component]) ||
            !isfinite(coor_prev[component])) {
            return 0;
        }
    }

    const double length_next = hypot(coor_next[0] - coor_node[0],
                                     coor_next[1] - coor_node[1]);
    const double length_prev = hypot(coor_prev[0] - coor_node[0],
                                     coor_prev[1] - coor_node[1]);
    if (!isfinite(length_next) || !isfinite(length_prev) ||
        length_next <= 0.0 || length_prev <= 0.0) {
        return 0;
    }

    /*
     * Keep the existing inward normal convention. For vertices ordered CCW,
     * the former Pl was the next vertex and Pr was the previous vertex.
     * The sum of their half-edge normals is
     *     C = 0.5 * Rcw(coor_prev - coor_next),
     * where Rcw(x, y) = (y, -x). The positive length and unit normal satisfy
     *     length_corner * normal_corner = C.
     */
    const double directed_corner[2] = {
        0.5 * (coor_prev[1] - coor_next[1]),
        0.5 * (coor_next[0] - coor_prev[0])
    };
    const double length = hypot(directed_corner[0], directed_corner[1]);
    if (!isfinite(length) || length <= 0.0) {
        return 0;
    }

    *length_corner = length;
    normal_corner[0] = directed_corner[0] / length;
    normal_corner[1] = directed_corner[1] / length;
    return 1;
}

double calCellArea(const double coor_sw[2], const double coor_se[2],
                   const double coor_ne[2], const double coor_nw[2])
{
    /* Shoelace area, translated to SW to avoid dependence on the origin. */
    const double vec_se[2] = {
        coor_se[0] - coor_sw[0], coor_se[1] - coor_sw[1]
    };
    const double vec_ne[2] = {
        coor_ne[0] - coor_sw[0], coor_ne[1] - coor_sw[1]
    };
    const double vec_nw[2] = {
        coor_nw[0] - coor_sw[0], coor_nw[1] - coor_sw[1]
    };

    return 0.5 * (vec_se[0] * vec_ne[1] - vec_se[1] * vec_ne[0] +
                  vec_ne[0] * vec_nw[1] - vec_ne[1] * vec_nw[0]);
}

void calCellCenter(const double coor_sw[2], const double coor_se[2],
                   const double coor_ne[2], const double coor_nw[2],
                   double center_cell[2])
{
    center_cell[0] = 0.25 * (coor_sw[0] + coor_se[0] +
                            coor_ne[0] + coor_nw[0]);
    center_cell[1] = 0.25 * (coor_sw[1] + coor_se[1] +
                            coor_ne[1] + coor_nw[1]);
}
