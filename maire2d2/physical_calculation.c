#include "linear_algebra_2d.h"
#include "mesh_calculation.h"
#include "physical_calculation.h"
#include <math.h>




double calImpendance(double density_cell, double pressure_cell, double GAMMA)
{
	double acous_speed = sqrt(GAMMA * pressure_cell / density_cell);
	double impendance = acous_speed * density_cell;

	return impendance;
}


void calMatCorner_glace(double length_corner, double normal_corner[2], 
						double impendance, double matrix_corner[][2])
{
	double tensor_production[2][2];
	vecCrossVec2(normal_corner, normal_corner, 
		tensor_production);

	numMulMat2(impendance * length_corner, tensor_production,
		matrix_corner);
}


void calMatCorner_eucclhyd(const double coor_p[2], const double coor_next[2],
						 const double coor_prev[2], double impendance,
						 double matrix_corner[][2])
{
	/* Directed CCW edges: prev -> p and p -> next. */
	const double edge[2][2] = {
		{coor_p[0] - coor_prev[0], coor_p[1] - coor_prev[1]},
		{coor_next[0] - coor_p[0], coor_next[1] - coor_p[1]}
	};
	double m_xx = 0.0, m_xy = 0.0, m_yy = 0.0;

	for (size_t side = 0; side < 2; ++side) {
		const double length = hypot(edge[side][0], edge[side][1]);
		const double n_x = edge[side][1] / length;
		const double n_y = -edge[side][0] / length;
		const double weight = 0.5 * impendance * length;

		/* M_pc = z_c * (l^- n^- tensor n^- + l^+ n^+ tensor n^+). */
		m_xx += weight * n_x * n_x;
		m_xy += weight * n_x * n_y;
		m_yy += weight * n_y * n_y;
	}

	matrix_corner[0][0] = m_xx;
	matrix_corner[0][1] = m_xy;
	matrix_corner[1][0] = m_xy;
	matrix_corner[1][1] = m_yy;
}


void calMatCorner(CornerScheme scheme, double coor_p[2], double coor_next[2],
				  double coor_prev[2], double impendance, double matrix_corner[][2])
{
	switch (scheme) {
	case CORNER_SCHEME_GLACE: {
		const double length_corner = calCornerLength(coor_p, coor_next, coor_prev);
		double normal_corner[2];
		calCornerNormal(coor_p, coor_next, coor_prev, normal_corner);
		calMatCorner_glace(length_corner, normal_corner, impendance, matrix_corner);
		break;
	}
	case CORNER_SCHEME_EUCCLHYD:
		calMatCorner_eucclhyd(coor_p, coor_next, coor_prev, impendance, matrix_corner);
		break;
	default:
		/* An invalid selector must not silently choose a different method. */
		matrix_corner[0][0] = matrix_corner[0][1] = NAN;
		matrix_corner[1][0] = matrix_corner[1][1] = NAN;
	}
}


void calCornerForce(double corner_force[2], double matrix_corner[][2],
					double velocity_cell[2], double pressure_cell,
					double length_corner,   double normal_corner[2],
					double velocity_node[2])
{
	double pressure_force[2];
	numMulVec2(length_corner * pressure_cell, normal_corner, 
			pressure_force);

	double num_viscosity[2];
	double relative_velocity[2];
	vecSubVec2(velocity_node, velocity_cell, relative_velocity);
	matMulVec2(matrix_corner, relative_velocity, num_viscosity);

	vecSubVec2(pressure_force, num_viscosity, corner_force);
}




double calCellDensity(double volume_cell , double mass_cell)
{
	return mass_cell / volume_cell;
}

double calCellSpInterior(double sp_energy_cell, double sp_momentum_cell[2])
{
	return sp_energy_cell - 0.5 * pow(calNorm2(sp_momentum_cell), 2);
}

double calCellPressure(double density_cell, double internal_energy, double GAMMA)
{
	return density_cell * (GAMMA-1) * internal_energy;
}
