#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "../linear_algebra_2d.h"

int main()
{
	double vec0[2] = {1, 2};
	double vec1[2] = {2, 3};
	double vec2[2] = {5, 7};

	double result_vec[2];
	double result_mat[2][2];

	vecAddVec2(vec1, vec2, result_vec);
	printf("vec1 add vec2 : (%f, %f)\n", result_vec[0], result_vec[1]);

	vecSubVec2(result_vec, vec2, result_vec);
	printf("vec1 add vec2 sub vec2: (%f, %f)\n", result_vec[0], result_vec[1]);

	numMulVec2(3, vec2, result_vec);
	printf("3 mul vec2 : (%f, %f)\n", result_vec[0], result_vec[1]);

	double result = vecMulVec2(vec1, vec2);
	printf("vec1 mul vec2 : %f\n", result);

	vecCrossVec2(vec1, vec2, result_mat);
	printf("vec1 cross vec2 : (%f, %f; %f, %f)\n", 
				result_mat[0][0], result_mat[0][1],
				result_mat[1][0], result_mat[1][1]);

	result = calNorm2(vec2);
	printf("norm of vec2 : %f\n", result);

	calRHNormal2(vec2, result_vec);
	printf("rh normal of vec2 : (%f, %f)\n", result_vec[0], result_vec[1]);
	
	attVec2(vec0, vec2);
	printf("attach vec2 to vec0 get: (%f, %f)\n", vec0[0], vec0[1]);


	double mat0[2][2] = {{2, 3}, {5, 7}};
	double mat1[2][2] = {{9, 11}, {13, 17}};

	numMulMat2(3, mat0, result_mat);
	printf("3 mul mat0 : (%f, %f; %f, %f)\n", 
				result_mat[0][0], result_mat[0][1],
				result_mat[1][0], result_mat[1][1]);

	matMulVec2(mat0, vec2, result_vec);
	printf("mat0 mul vec2 : (%f, %f)\n", result_vec[0], result_vec[1]);

	attMat2(mat1, mat0);
	printf("attach mat0 to mat1 : (%f, %f; %f, %f)\n", 
				mat1[0][0], mat1[0][1],
				mat1[1][0], mat1[1][1]);

	double inverse_mat[2][2];
	inverseMatrix2(mat0, inverse_mat);
	for (size_t row = 0; row < 2; ++row) {
		for (size_t col = 0; col < 2; ++col) {
			double product = mat0[row][0] * inverse_mat[0][col]
			               + mat0[row][1] * inverse_mat[1][col];
			assert(fabs(product - (row == col ? 1.0 : 0.0)) < 1e-12);
		}
	}
	printf("inverse of mat0: (%f, %f; %f, %f)\n", 
				inverse_mat[0][0], inverse_mat[0][1],
				inverse_mat[1][0], inverse_mat[1][1]);

	 
}
