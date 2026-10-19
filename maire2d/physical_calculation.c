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
					

void calCornerForce(double corner_force[2], double impendance,
					double velociy_cell[2], double pressure_cell,
					double length_corner,   double normal_corner[2],
					double velociy_node[2])
{
	double gredient_pressure[2];
	numMulVec2(length_corner * pressure_cell, normal_corner, 
			gredient_pressure);

	double num_viscosity[2];
	double relative_velocity[2];
	vecSubVec2(velociy_node, velociy_cell, relative_velocity);

	double matrix_corner[2][2];
// 这里可以加一个选项
	calMatCorner_glace(length_corner, normal_corner, impendance, matrix_corner);

	matMulVec2(matrix_corner, relative_velocity, num_viscosity);

	vecSubVec2(gredient_pressure, num_viscosity, corner_force);
}
