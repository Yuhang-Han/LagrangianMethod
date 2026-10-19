#include "linear_algebra_2d.h"
#include <math.h>
void attVec2(double vector1[2], double vector2[2])
{
	vector1[0] += vector2[0];
	vector1[1] += vector2[1];
}


void attMat2(double matrix1[][2], double matrix2[][2])
{
	matrix1[0][0] += matrix2[0][0];
	matrix1[0][1] += matrix2[0][1];
	matrix1[1][0] += matrix2[1][0];
	matrix1[1][1] += matrix2[1][1];
}


void vecAddVec2(double vector1[2], double vector2[2], double result_vector[2])
{
	result_vector[0] = vector1[0] + vector2[0];
	result_vector[1] = vector1[1] + vector2[1];
}

void vecSubVec2(double vector1[2], const double vector2[2], double result_vector[2])
{
	result_vector[0] = vector1[0] - vector2[0];
	result_vector[1] = vector1[1] - vector2[1];
}


void numMulVec2(double num, double vector[2], double result_vector[2])
{
	result_vector[0] = num * vector[0];
	result_vector[1] = num * vector[1];
}

void numMulMat2(double num, double matrix[][2], double result_matrix[][2])
{
	result_matrix[0][0] = num * matrix[0][0];
	result_matrix[0][1] = num * matrix[0][1];
	result_matrix[1][0] = num * matrix[1][0];
	result_matrix[1][1] = num * matrix[1][1];
}


double vecMulVec2(double vector1[2], double vector2[2])
{
	return vector1[0] * vector2[0] + vector1[1] * vector2[1];
}

void matMulVec2(double matrix[][2], double vector[2], double result_vector[2])
{
	result_vector[0] = matrix[0][0] * vector[0] \
					 + matrix[0][1] * vector[1];
	result_vector[1] = matrix[1][0] * vector[0] \
					 + matrix[1][1] * vector[1];
}


void vecCrossVec2(double vector1[2], double vector2[2], double result_matrix[][2])
{
	result_matrix[0][0] = vector1[0] * vector2[0];
	result_matrix[0][1] = vector1[0] * vector2[1];
	result_matrix[1][0] = vector1[1] * vector2[0];
	result_matrix[1][1] = vector1[1] * vector2[1];
}


double calNorm2(double vector[2])
{
	return sqrt( vecMulVec2(vector, vector) );
}


void calRHNormal2(double vector[2], double normal[2])
{
	double normal_rh[2];
	normal_rh[0] = vector[1];
	normal_rh[1] = -vector[0];

	double norm = calNorm2(normal_rh);

	numMulVec2(1.0/norm, normal_rh, normal); 
}


void inverseMatrix2(double matrix[][2], double result_matrix[][2])
{
	double det = matrix[0][0] * matrix[1][1] \
				- matrix[0][1] * matrix[1][0];
	double matrix_tmp[2][2] = { {matrix[1][1], -matrix[0][1]},
								{-matrix[1][0], matrix[0][0]} };

	numMulMat2(1.0/det, matrix_tmp, result_matrix);
}
