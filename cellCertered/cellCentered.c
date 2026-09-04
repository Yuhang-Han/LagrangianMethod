#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include <math.h>

int main(int argc, char *argv[]) 
{

// Read the mesh size -------------------
	if(argc != 2){
		printf("Usage: staggeredGrid space_grid\n");
		return 1;
	}

    char *end;
    
    unsigned long long value1 = strtoull(argv[1], &end, 10);
    if (*end != '\0' || value1 > SIZE_MAX) {
        fprintf(stderr, "invalid size: %s\n", argv[1]);
        return 1;
    }

    const size_t N = (size_t)value1;

    printf("space grid size = %zu\n", N);

	if(N % 2 == 1) {
		printf("N should be even number.\n");
		return 1;
	}


	const double X = 1;
	const double T = 0.1;

	const double GAMMA = 1.4;

	const double H   = X/N;

	// const double EPS = 0.5;
	const double CFL = 0.3;

	const double xSeparate = (0.0 + X) / 2;


// Initialize ----------------------------

	// Eulerian coordiantes x(X)
	double *xEuler = malloc((N+1)*sizeof(double));		
	for(size_t i = 0; i < N+1; ++i) {
		xEuler[i] = i*H;
	}


	double *Vn = malloc((N+2)*sizeof(double));
	double *un = malloc((N+2)*sizeof(double));
	double *pn = malloc((N+2)*sizeof(double));
	double *En = malloc((N+2)*sizeof(double));
	double *V  = malloc((N+2)*sizeof(double));
	double *u  = malloc((N+2)*sizeof(double));
	double *p  = malloc((N+2)*sizeof(double));
	double *E  = malloc((N+2)*sizeof(double));

	double *m  = malloc((N+2)*sizeof(double));

	// sin() IV for cal orders
/*
	for(size_t iGhost = 0; iGhost < N+2; ++iGhost) {
		double xR = iGhost * H;
		double xL = xR - H;
		
		double rho = 1 + EPS/(H*2*M_PI) \
					* ( cos(2*M_PI*xL) - cos(2*M_PI*xR) );
		Vn[iGhost] = 1/rho;
		pn[iGhost] = 1 + GAMMA*(rho-1);
		un[iGhost] = EPS*sqrt(GAMMA)/(H*2*M_PI) \
					* ( cos(2*M_PI*xL) - cos(2*M_PI*xR) );
		En[iGhost] = 0.5 * pow(un[iGhost],2) + pn[iGhost]*Vn[iGhost] / (GAMMA-1);

		m[iGhost] = H * rho;
	}
*/

	// Sod shock tube	
	for(size_t iGhost = 0; iGhost < N+2; ++iGhost) {
		double xCenter = (iGhost-0.5) * H;
		if (xCenter <= xSeparate) {
			double rho= 1;
			Vn[iGhost] = 1/rho;
			un[iGhost] = 0;
			pn[iGhost] = 1;
			En[iGhost] = 0.5 * pow(un[iGhost],2) + pn[iGhost]*Vn[iGhost] / (GAMMA-1);
			m[iGhost] = H * rho;
		}
		else {
			double rho = 0.125;
			Vn[iGhost] = 1/rho;
			un[iGhost] = 0;
			pn[iGhost] = 0.1;
			En[iGhost] = 0.5 * pow(un[iGhost],2) + pn[iGhost]*Vn[iGhost] / (GAMMA-1);
			m[iGhost] = H * rho;
		}
	}

	
/*
	printf("U(0): ");
	for (size_t j = 0; j < N+2; j++) {
		printf("%.3f, ", un[j]);
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
	
	double pn_half[N+1], un_half[N+1];		// states on interfaces
	double an[N+2];							// local acoustic speed
	double t = 0;							// timer

	size_t it = 0;

	while (t < T) {
		for(size_t iGhost = 0; iGhost < N+2; ++iGhost) {
			an[iGhost] = sqrt(GAMMA*pn[iGhost]*Vn[iGhost]);
		}
		
		// Set time step length
		// const double tau = T / 10000;
		 ++it;
			

		double tau = T;			
		for(size_t iCell = 1; iCell < N+1; ++iCell){
			double tau_local = CFL * H / (fabs(un[iCell]) + an[iCell]);  
			//printf("tau_local: %f\n", tau_local);
			tau = fmin(tau_local, tau);
		}
		//printf("tau: %f\n", tau);

		if (t + tau > T) {
			tau = T - t;
			t = T;
		}
		else {
			t += tau;
		}

		// Solve interface states by acoustic approximatation
		for(size_t iFace = 0; iFace < N+1; ++iFace) {
			double zL = an[iFace]   / Vn[iFace];		// acoustic impedance
			double zR = an[iFace+1] / Vn[iFace+1];
			
			un_half[iFace] = (zL*un[iFace] + zR*un[iFace+1]) / (zL + zR) \
							-(pn[iFace+1] - pn[iFace]) / (zL + zR);
							
			pn_half[iFace] = (zL*pn[iFace] + zR*pn[iFace+1]) / (zL + zR) \
							-(un[iFace+1] - un[iFace]) * (zL*zR) / (zL + zR);
		}

		// Update x's coordinates
		for(size_t iFace = 0; iFace < N+1; ++iFace) {
			xEuler[iFace] = xEuler[iFace] + tau*un_half[iFace];
		}
		
		// Advance cell states
		for(size_t iCell = 1; iCell < N+1; ++iCell) {
			V[iCell] = Vn[iCell] \
						+ tau/m[iCell] * (un_half[iCell] - un_half[iCell-1]);
			u[iCell] = un[iCell] 
						- tau/m[iCell] * (pn_half[iCell] - pn_half[iCell-1]);
			E[iCell] = En[iCell]\
						- tau/m[iCell] * ( pn_half[iCell]*un_half[iCell]
											- pn_half[iCell-1]*un_half[iCell-1] );
			p[iCell] = (GAMMA-1) * (E[iCell] - 0.5*pow(u[iCell],2)) / V[iCell];
		}


		// Apply periodic boundary condition
		/*
		V[0] = V[N];	V[N+1] = V[1];
		u[0] = u[N]; 	u[N+1] = u[1];
		E[0] = E[N]; 	E[N+1] = E[1];
		p[0] = p[N]; 	p[N+1] = p[1];
		*/

		// Free boundary condition
		V[0] = V[1];	V[N+1] = V[N];
		u[0] = u[1]; 	u[N+1] = u[N];
		E[0] = E[1]; 	E[N+1] = E[N];
		p[0] = p[1]; 	p[N+1] = p[N];

		
		// Rotate the pointers
		double *tmp;
		
		tmp = Vn;  Vn = V;  V = tmp;
		tmp = un;  un = u;	u = tmp;
		tmp = En;  En = E;  E = tmp;
		tmp = pn;  pn = p;  p = tmp;

/*
		if(it % 1 == 0) {
			printf("U(%f): ", t);
			for (size_t j = 0; j < N+2; j++) {
				printf("%.3f, ", un[j]);
			}
			printf("\n\n");

			printf("rho(%f): ", t);
			for (size_t j = 0; j < N+2; j++) {
				printf("%.6f, ", 1.0/Vn[j]);
			}
			printf("\n\n");

			printf("E(%f): ", t);
			for (size_t j = 0; j < N+2; j++) {
				printf("%.3f, ", En[j]);
			}
			printf("\n\n");
			
			printf("p(%f): ", t);
			for (size_t j = 0; j < N+2; j++) {
				printf("%.3f, ", pn[j]);
			}
			printf("\n\n");


			printf("xEuler(%f): ", t);
			for (size_t j = 0; j < N+1; j++) {
				printf("%.3f, ", xEuler[j]);
			}
			printf("\n\n");

		}
*/
	}

	printf("time step numer: %zu \n", it);


		
/******************
* Output
*******************/

	char filename_Rho[256];
	char filename_U[256];
	char filename_P[256];

	snprintf(filename_Rho, sizeof(filename_Rho), 
			"build/output/cellCentered_N%zu_Rho.csv", N);
	snprintf(filename_U, sizeof(filename_U), 
			"build/output/cellCentered_N%zu_U.csv", N);
	snprintf(filename_P, sizeof(filename_P), 
			"build/output/cellCentered_N%zu_P.csv", N);

	FILE *fpRho = fopen(filename_Rho, "w");
	FILE *fpU = fopen(filename_U, "w");
	FILE *fpP = fopen(filename_P, "w");


	for (size_t iCell = 1; iCell < N+1; iCell++) {
		double xCenter = 0.5 * (xEuler[iCell-1] + xEuler[iCell]);
		
	    fprintf(fpRho, 	"%.16f,%.16f\n", xCenter, 1.0/Vn[iCell]);
	    fprintf(fpU, 	"%.16f,%.16f\n", xCenter, un[iCell]);
	    fprintf(fpP, 	"%.16f,%.16f\n", xCenter, pn[iCell]);
	}

	fclose(fpRho);
	fclose(fpU);
	fclose(fpP);
	

	free(xEuler);
	free(Vn);	free(un); 	free(pn);	free(En);
	free(V);	free(u); 	free(p);	free(E);
	free(m);
}
