#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include <math.h>

int main(int argc, char *argv[]) 
{

// Read the mesh size -------------------
	if(argc != 3){
		printf("Usage: staggeredGrid time_gird space_grid\n");
		return 1;
	}

    char *end;
    
    unsigned long long value1 = strtoull(argv[1], &end, 10);
    if (*end != '\0' || value1 > SIZE_MAX) {
        fprintf(stderr, "invalid size: %s\n", argv[1]);
        return 1;
    }
    unsigned long long value2 = strtoull(argv[2], &end, 10);
    if (*end != '\0' || value2 > SIZE_MAX) {
        fprintf(stderr, "invalid size: %s\n", argv[2]);
        return 1;
    }

    const size_t M = (size_t)value1;
    const size_t N = (size_t)value2;

    printf("time grid size = %zu\n", M);
    printf("space grid size = %zu\n", N);

	if(M % 2 == 1 || N % 2 == 1) {
		printf("M and N should be even number.\n");
		return 1;
	}


	const double X = 1;
	const double T = 0.1;

	const double C = 1.0;
	const double GAMMA = 1.4;

	if(M % 2 == 1 || N % 2 == 1) {
		printf("M and N should be even number.\n");
		return 1;
	}

	const double TAU = T/M;
	const double H   = X/N;

	const double EPS = 0.5;

	const double xSeparate = (0.0 + X) / 2;


// Initialize ----------------------------
	// set the separate at l = N/2
	double xEuler[N+1];			// Eulerian coordiantes X(x)
	for(size_t i = 0; i < N/2+1; ++i) {
		xEuler[N/2 + i] = xSeparate + i*H;
		xEuler[N/2 - i] = xSeparate - i*H;
	}
/*
	printf("xEuler(0): ");
	for (size_t j = 0; j < N+1; j++) {
		printf("%f, ", xEuler[j]);
	}
	printf("\n");
*/
	double rho0[N+2];
	// p_{l-1/2}^n, q^{n-1/2}_{l-1/2}, U^{n-1/2}_l, V^n_{l-1/2}
	double buffer1[N+2], buffer2[N+2], buffer3[N+1], buffer4[N+2];
	double buffer5[N+2], buffer6[N+2], buffer7[N+1], buffer8[N+2];
	double *pn = buffer1; double *qn = buffer2; 
	double *Un = buffer3; double *Vn = buffer4;
	double *p  = buffer5; double *q  = buffer6; 
	double *U  = buffer7; double *V  = buffer8;

// Sod shock tube	
	for(size_t ix = 0; ix < N+2; ++ix) {
		double x = (ix-0.5) * H;
		if (x <= xSeparate) {
			rho0[ix] = 1;
			Vn[ix] = 1/rho0[ix];
			pn[ix] = 1;
			qn[ix] = 0;
		}
		else {
			rho0[ix] = 0.125;
			Vn[ix] = 1/rho0[ix];
			pn[ix] = 0.1;
			qn[ix] = 0;
		}

	}
	for(size_t i = 0; i < N+1; ++i){
		double U0 = 0;
		U[i] = U0 - TAU/((rho0[i]+rho0[i+1]) * H) * (pn[i+1] - pn[i]);
	}

/*
	for(size_t ix = 0; ix < N+2; ++ix) {
		double x = (ix-0.5) * H;
		rho0[ix] = 1 + EPS/(H*2*M_PI) \
					* (cos(2*M_PI*(x-0.5*H)) - cos(2*M_PI*(x+0.5*H)));
		Vn[ix] = 1/rho0[ix];
		pn[ix] = 1 + GAMMA*(rho0[ix]-1);
		qn[ix] = 0;

	}

	for(size_t i = 0; i < N+1; ++i){
		double x = i * H;
		double U0 = 1 + EPS*sqrt(GAMMA)*sin(2*M_PI*x);
		U[i] = U0 - TAU/((rho0[i]+rho0[i+1]) * H) * (pn[i+1] - pn[i]);
	}
*/


/*
	printf("U(1/2): ");
	for (size_t j = 0; j < N+1; j++) {
		printf("%.3f, ", U[j]);
	}
	printf("\n\n");
	
	printf("V(0): ");
	for (size_t j = 0; j < N+2; j++) {
		printf("%.6f, ", Vn[j]);
	}
	printf("\n\n");
	
	printf("p(0): ");
	for (size_t j = 0; j < N+2; j++) {
		printf("%.3f, ", pn[j]);
	}
	printf("\n\n");

*/



// Time advance ---------------------------
	for(size_t it = 1; it < M+1; ++it) {
	
		if(it > 1) {
			for(size_t ix = 0; ix < N+1; ++ix) {			
				U[ix] = Un[ix] - TAU/((rho0[ix]+rho0[ix+1])/2 * H) \
						* (pn[ix+1] + qn[ix+1] - pn[ix] - qn[ix]);
			}		
		}

		for(size_t ix = 0; ix < N+1; ++ix){
			xEuler[ix] = xEuler[ix] + TAU*U[ix];
		}

		for(size_t ix = 0; ix < N; ++ix){
			V[ix+1] = Vn[ix+1] + TAU/(rho0[ix+1] * H) \
					* (U[ix+1] - U[ix]);
		}
		V[0] = V[1];	V[N+1] = V[N];		// free boundary
		//V[0] = V[N];	V[N+1] = V[1];		// periodic boundary
		
		for(size_t ix = 0; ix < N; ++ix){
		/*
			q[ix+1] = -2*C*C/(Vn[ix+1] + V[ix+1]) \
					* (U[ix+1] -U[ix]) * fabs(U[ix+1] -U[ix]);
		*/
		// test q = 0
			q[ix+1] = 0;
		}
		q[0] = q[1];	q[N+1] = q[N];
		//q[0] = q[N];	q[N+1] = q[1];
		
		for(size_t ix = 0; ix < N; ++ix){
			p[ix+1] = 1/((GAMMA+1)*V[ix+1] - (GAMMA-1)*Vn[ix+1]) \
					* ( \
						((GAMMA+1)*Vn[ix+1] - (GAMMA-1)*V[ix+1]) * pn[ix+1] \
						- 2*(GAMMA-1) * q[ix+1] * (V[ix+1] - Vn[ix+1])
					);
		}
		p[0] = p[1];	p[N+1] = p[N];
		//p[0] = p[N]; 	p[N+1] = p[1];

		// Rotate the pointers
		double *tmp;
		tmp = Un;  Un = U;	U = tmp;
		tmp = Vn;  Vn = V;  V = tmp;
		tmp = pn;  pn = p;  p = tmp;
		tmp = qn;  qn = q;  q = tmp;

/*
		if(it % 1 == 0) {
			printf("U(%zu-1/2): ", it);
			for (size_t j = 0; j < N+1; j++) {
				printf("%.3f, ", Un[j]);
			}
			printf("\n\n");

			printf("rho(%zu): ", it);
			for (size_t j = 0; j < N+2; j++) {
				printf("%.6f, ", 1.0/Vn[j]);
			}
			printf("\n\n");

			printf("q(%zu): ", it);
			for (size_t j = 0; j < N+2; j++) {
				printf("%.3f, ", qn[j]);
			}
			printf("\n\n");
			
			printf("p(%zu): ", it);
			for (size_t j = 0; j < N+2; j++) {
				printf("%.3f, ", pn[j]);
			}
			printf("\n\n");


			printf("xEuler(%zu): ", it);
			for (size_t j = 0; j < N+1; j++) {
				printf("%.3f, ", xEuler[j]);
			}
			printf("\n\n");

		}
*/

	}


/******************
* Output
*******************/

	char filename_Rho[256];
	char filename_U[256];
	char filename_P[256];

	snprintf(filename_Rho, sizeof(filename_Rho), 
			"build/output/staggeredGrid_M%zu_N%zu_Rho.csv", M, N);
	snprintf(filename_U, sizeof(filename_U), 
			"build/output/staggeredGrid_M%zu_N%zu_U.csv", M, N);
	snprintf(filename_P, sizeof(filename_P), 
			"build/output/staggeredGrid_M%zu_N%zu_P.csv", M, N);

	FILE *fpRho = fopen(filename_Rho, "w");
	FILE *fpU = fopen(filename_U, "w");
	FILE *fpP = fopen(filename_P, "w");


	for (size_t iCell = 1; iCell < N+1; iCell++) {
		double xCenter = 0.5 * (xEuler[iCell-1] + xEuler[iCell]);
		
	    fprintf(fpRho, "%.16f,%.6f\n", xCenter, 1.0/Vn[iCell]);
	    fprintf(fpP, "%.16f,%.6f\n", xCenter, pn[iCell]);
	}
	for (size_t j = 0; j < N+1; j++) {
	    fprintf(fpU, "%.16f,%.6f\n", xEuler[j], Un[j]);
	}

	fclose(fpRho);
	fclose(fpU);
	fclose(fpP);
	
}
