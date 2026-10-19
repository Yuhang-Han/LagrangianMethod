#ifndef MESH_CALCULATION_H
#define MESH_CALCULATION_H

/*
 * Arithmetic for a moving quadrilateral mesh. Callers access mesh arrays
 * directly as [node_y][node_x] and pass the required coordinate vectors.
 * Every coordinate vector stores Cartesian components as [x, y].
 */

/*
 * Calculate corner length and unit normal in one operation.
 * next and prev are the next and previous vertices in CCW cell order.
 * The normal points inward, preserving the existing GLACE implementation:
 *     length_corner * normal_corner = 0.5 * Rcw(coor_prev - coor_next).
 * Return 1 on success. Return 0 for nonfinite coordinates, a zero adjacent
 * edge, a zero directed corner vector, or a null pointer. On failure,
 * both outputs are zero when both output pointers are valid.
 */
int calCornerGeometry(const double coor_node[2],
                      const double coor_next[2],
                      const double coor_prev[2],
                      double *length_corner,
                      double normal_corner[2]);

/* Signed polygon area; positive for a simple, nondegenerate CCW quadrilateral. */
double calCellArea(const double coor_sw[2], const double coor_se[2],
                   const double coor_ne[2], const double coor_nw[2]);

/*
 * Arithmetic average of the four vertices. This equals the area centroid
 * for the initial Cartesian cells, but not for an arbitrary deformed cell.
 */
void calCellCenter(const double coor_sw[2], const double coor_se[2],
                   const double coor_ne[2], const double coor_nw[2],
                   double center_cell[2]);

#endif
