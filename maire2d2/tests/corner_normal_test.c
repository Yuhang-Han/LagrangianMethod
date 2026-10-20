#include <stdio.h>
#include "../mesh_calculation.h"

int main() 
{
	double coor_p[2]  = {1., 1.};
	double coor_pl[2] = {0., 1.};
	double coor_pr[2] = {1., 0.};

	double normal_corner[2];

	calCornerNormal(coor_p, coor_pl, coor_pr, normal_corner);

	printf("corner normal: (%f, %f)\n", normal_corner[0], normal_corner[1]);
}
