/***************************
* 2D linear algebra
****************************/

#ifndef LINEAR_ALGEBRA_2D_H
#define LINEAR_ALGEBRA_2D_H


#include <stdio.h>


// attach
void attVec2(double vector1[2], double vector2[2]);
void attMat2(double matrix1[][2], double matrix2[][2]);

void vecAddVec2(double vector1[2], double vector2[2], double result_vector[2]);
void vecSubVec2(double vector1[2], const double vector2[2], double result_vector[2]);

void numMulVec2(double num, double vector[2], double result_vector[2]);
void numMulMat2(double num, double matrix[][2], double result_matrix[][2]);

double vecMulVec2(double vector1[2], double vector2[2]);
void matMulVec2(double matrix[][2], double vector[2], double result_vector[2]);

void vecCrossVec2(double vector1[2], double vector2[2], double result_matrix[][2]);

double calNorm2(double vector[2]);

// calculate right-hand unit normal vector 
void calRHNormal2(double vector[2], double normal[2]);

void inverseMatrix2(double matrix[][2], double result_matrix[][2]);


#endif
