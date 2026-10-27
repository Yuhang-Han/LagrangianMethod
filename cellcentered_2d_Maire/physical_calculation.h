#ifndef PHYSICAL_CALCULATION_H
#define PHYSICAL_CALCULATION_H

#include "linear_algebra_2d.h"
#include "mesh_calculation.h"
#include <math.h>

double calImpendance(double density_cell, double pressure_cell, double GAMMA);

typedef enum {
	CORNER_SCHEME_GLACE,
	CORNER_SCHEME_EUCCLHYD
} CornerScheme;

void calMatCorner_glace(double length_corner, double normal_corner[2], double impendance,
					double matrix_corner[][2]);

/*
 * First-order acoustic EUCCLHYD (Maire HDR, Eq. 3.84, p. 90).
 * coor_prev, coor_p, coor_next follow the cell's counterclockwise boundary.
 * The two incident edges must have nonzero length; impendance is rho * c.
 * Each edge contributes its HALF length times z * (n tensor n).
 * The output is overwritten and is positive definite for positive z and
 * non-collinear edges.
 */
void calMatCorner_eucclhyd(const double coor_p[2], const double coor_next[2],
						 const double coor_prev[2], double impendance,
						 double matrix_corner[][2]);

/* Use the same scheme for the nodal solve and the corner-force calculation. */
void calMatCorner(CornerScheme scheme, double coor_p[2], double coor_next[2],
				  double coor_prev[2], double impendance, double matrix_corner[][2]);

/* F_pc = P_c * l_pc * n_pc - M_pc * (U_p - U_c), Maire HDR Eq. 3.83. */
void calCornerForce(double corner_force[2], double matrix_corner[][2],
					double velocity_cell[2], double pressure_cell,
					double length_corner,   double normal_corner[2],
					double velocity_node[2]);

double calCellDensity(double volume_cell , double mass_cell);
double calCellSpInterior(double sp_energy_cell, double sp_momentum_cell[2]);
double calCellPressure(double density_cell, double internal_energy, double GAMMA);

#endif
