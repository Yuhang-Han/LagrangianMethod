/***************************
* 2D linear algebra
****************************/


#include <stdio.h>

void attVec2(double vector1[2], double vector2[2]);
void attMat2(double matrix1[][2], double matrix2[][2]);

void vecAddVec2(double vector1[2], double vector2[2], double result_vector[2]);
void vecMinVec2(double vector1[2], double vector2[2], double result_vector[2]);

void numMulVec2(double num, double vector[2], double result_vector[2]);
void numMulMat2(double num, double matrix[][2], double result_matrix[][2]);

void vecMulVec2(double vector1[2], double vector2[2], double result_number);
void matMulVec2(double matrix[][2], double vector[2], double result_vector[2]);

void vecCrossVec2(double vector1[2], double vector2[2], double result_matrix[][2]);

void inverseMatrix2(double matrix[][2]);
