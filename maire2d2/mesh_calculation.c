#include <stdio.h>

#include "linear_algebra_2d.h"
#include "mesh_calculation.h"

/*
 * The node must be a vertex of cell (cell_x, cell_y); direction is +1 for
 * the next CCW vertex or -1 for the previous vertex. Output pointers must
 * be valid. As before, the caller is responsible for array bounds.
 */
void calNeighVertexIdx(size_t node_x, size_t node_y,
                      size_t cell_x, size_t cell_y, int direction,
                      size_t *neigh_node_x, size_t *neigh_node_y)
{
    if (node_x == cell_x && node_y == cell_y) {
        /* Southwest vertex. */
        *neigh_node_x = direction == 1 ? cell_x + 1 : cell_x;
        *neigh_node_y = direction == 1 ? cell_y : cell_y + 1;
    } else if (node_x == cell_x + 1 && node_y == cell_y) {
        /* Southeast vertex. */
        *neigh_node_x = direction == 1 ? cell_x + 1 : cell_x;
        *neigh_node_y = direction == 1 ? cell_y + 1 : cell_y;
    } else if (node_x == cell_x + 1 && node_y == cell_y + 1) {
        /* Northeast vertex. */
        *neigh_node_x = direction == 1 ? cell_x : cell_x + 1;
        *neigh_node_y = direction == 1 ? cell_y + 1 : cell_y;
    } else {
        /* Northwest vertex, under the documented incident-node precondition. */
        *neigh_node_x = direction == 1 ? cell_x : cell_x + 1;
        *neigh_node_y = direction == 1 ? cell_y : cell_y + 1;
    }
}


double calCornerLength(double coor_p[2], double coor_pl[2], double coor_pr[2])
{
	double vec_plp[2], vec_ppr[2];
	vecSubVec2(coor_p, coor_pl, vec_plp);
	vecSubVec2(coor_pr, coor_p, vec_ppr);

	double length_plp = calNorm2(vec_plp);
	double length_ppr = calNorm2(vec_ppr);

	double normal_plp[2], normal_ppr[2];
	calRHNormal2(vec_plp, normal_plp);
	calRHNormal2(vec_ppr, normal_ppr);

	double half_directed_normal_plp[2], half_directed_normal_ppr[2];
	numMulVec2(length_plp/2.0, normal_plp, half_directed_normal_plp);
	numMulVec2(length_ppr/2.0, normal_ppr, half_directed_normal_ppr);

	double directed_length_corner[2];
	vecAddVec2(half_directed_normal_plp, half_directed_normal_ppr,
				directed_length_corner);

	double length_corner = calNorm2(directed_length_corner);

	return length_corner;
}



void calCornerNormal(double coor_p[2], double coor_pl[2], double coor_pr[2],
						double normal_corner[2])
{
	double vec_plp[2], vec_ppr[2];
	vecSubVec2(coor_p, coor_pl, vec_plp);
	vecSubVec2(coor_pr, coor_p, vec_ppr);

	double length_plp = calNorm2(vec_plp);
	double length_ppr = calNorm2(vec_ppr);

	double normal_plp[2], normal_ppr[2];
	calRHNormal2(vec_plp, normal_plp);
	calRHNormal2(vec_ppr, normal_ppr);

	double half_directed_normal_plp[2], half_directed_normal_ppr[2];
	numMulVec2(length_plp/2.0, normal_plp, half_directed_normal_plp);
	numMulVec2(length_ppr/2.0, normal_ppr, half_directed_normal_ppr);

	double directed_length_corner[2];
	vecAddVec2(half_directed_normal_plp, half_directed_normal_ppr,
				directed_length_corner);

	double length_corner = calNorm2(directed_length_corner);

	numMulVec2(1.0/length_corner, directed_length_corner,
				normal_corner);
}
