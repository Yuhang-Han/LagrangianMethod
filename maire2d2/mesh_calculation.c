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


double calCellVolume(const double coor_sw[2], const double coor_se[2],
                     const double coor_ne[2], const double coor_nw[2])
{
    /* Half the cross product of the diagonals is the signed polygon area.
     * Coordinate differences avoid cancellation of large absolute products
     * when a small cell is translated far from the origin.
     */
    const double diagonal_0_x = coor_ne[0] - coor_sw[0];
    const double diagonal_0_y = coor_ne[1] - coor_sw[1];
    const double diagonal_1_x = coor_nw[0] - coor_se[0];
    const double diagonal_1_y = coor_nw[1] - coor_se[1];

    return 0.5 * (diagonal_0_x * diagonal_1_y
                - diagonal_0_y * diagonal_1_x);
}


double calCornerLength(double coor_p[2], double coor_pl[2], double coor_pr[2])
{
	double vec_ppl[2], vec_prp[2];
	vecSubVec2(coor_pl, coor_p, vec_ppl);
	vecSubVec2(coor_p, coor_pr, vec_prp);

	double length_ppl = calNorm2(vec_ppl);
	double length_prp = calNorm2(vec_prp);

	double normal_ppl[2], normal_prp[2];
	calRHNormal2(vec_ppl, normal_ppl);
	calRHNormal2(vec_prp, normal_prp);

	double half_directed_normal_ppl[2], half_directed_normal_prp[2];
	numMulVec2(length_ppl/2.0, normal_ppl, half_directed_normal_ppl);
	numMulVec2(length_prp/2.0, normal_prp, half_directed_normal_prp);

	double directed_length_corner[2];
	vecAddVec2(half_directed_normal_ppl, half_directed_normal_prp,
				directed_length_corner);

	double length_corner = calNorm2(directed_length_corner);

	return length_corner;
}



// CCW: pl means next, pr means prev
void calCornerNormal(double coor_p[2], double coor_pl[2], double coor_pr[2],
						double normal_corner[2])
{
	double vec_ppl[2], vec_prp[2];
	vecSubVec2(coor_pl, coor_p, vec_ppl);
	vecSubVec2(coor_p, coor_pr, vec_prp);

	double length_ppl = calNorm2(vec_ppl);
	double length_prp = calNorm2(vec_prp);

	double normal_ppl[2], normal_prp[2];
	calRHNormal2(vec_ppl, normal_ppl);
	calRHNormal2(vec_prp, normal_prp);

	double half_directed_normal_ppl[2], half_directed_normal_prp[2];
	numMulVec2(length_ppl/2.0, normal_ppl, half_directed_normal_ppl);
	numMulVec2(length_prp/2.0, normal_prp, half_directed_normal_prp);

	double directed_length_corner[2];
	vecAddVec2(half_directed_normal_ppl, half_directed_normal_prp,
				directed_length_corner);

	double length_corner = calNorm2(directed_length_corner);

	numMulVec2(1.0/length_corner, directed_length_corner,
				normal_corner);
}
