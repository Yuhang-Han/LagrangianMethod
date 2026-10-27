#ifndef PHYSICAL_CALCULATION_H
#define PHYSICAL_CALCULATION_H

#include "linear_algebra_2d.h"
#include "mesh_calculation.h"
#include <math.h>

double calImpendance(double density_cell, double pressure_cell, double GAMMA);

void calMatCorner_glace(double length_corner, double normal_corner[2], double impendance,
					double matrix_corner[][2]);

void calCornerForce(double corner_force[2], double impendance,
					double velociy_cell[2], double pressure_cell,
					double length_corner,   double normal_corner[2],
					double velociy_node[2]);


#endif
