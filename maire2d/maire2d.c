#include <stdio.h>

#include "linear_algebra_2d.h"
#include "mesh_calculation.h"
#include "physical_calculation.h"


int main() 
{

/*******************
* Model Coefficients
*******************/

	const double GAMMA = 1.4;

	const double T = 1.0;

	const double CFL = 0.2;

	const double L = 1.0;
	const double X0 = 2.0*L;
	const double X1 = 2.0*L;

	const size_t NUM_CELL_VERTEXS = 4;
	const size_t NUM_NODE_NEIGH_CELLS = 4;

	const size_t NUM_CELLS_1D = 10;
	const size_t NUM_CELLS = NUM_CELLS_1D * NUM_CELLS_1D;
	const size_t NUM_NODES = (NUM_CELLS_1D + 1) * (NUM_CELLS_1D + 1);

/*******************
* Memory Allocation
*******************/





/*******************
* Initialization
*******************/



/********************
* Time Advancing
*********************/

// Calculate time step length: tau
	double t = 1e-4;


// Update
	while (t < T) {
		// apply boundary condition
		

		// update node velocity
		for(size_t iP = 0; iP < NUM_NODES; ++iP) {
			double node_matrix[2][2] = {0, };
			double pressure_gradient[2] = {0, 0};
			double ave_velocity_node[2] = {0, 0};
			
			for(size_t iCP = 0; iCP < NUM_NODE_NEIGH_CELLS; ++iCP) {
				// transform local cell idx to global
				size_t idxC = localToGlobalCellIdx(iP, iCP);
				// calculate corner length and normal by neigh nodes
				double length_corner  = calCornerLength(idxC, iP, meshinfo);
				double normal_corner[2] = calCornerNormal(idxC, iP, meshinfo);

				double* corner_matrix = calCornerMatrix();

				// 改成先数乘再加
				attachVector2(pressure_gradient, 
					numMulVector2(length_corner * pressure_cell[idxC], normal_corner));

				attachVector2(ave_velocity_node, 
					matMulVector(corner_matrix, sp_momentum[idxC]));

				attachMatrix2(node_matrix, corner_matrix);
			}
			
			double* node_matrix_inverse = inverseMatrix2(node_matrix);

			velocity_node[iP] = 
				matMulVector2(node_matrix_inverse, 
							attachVector2(pressure_gradient, ave_velocity_node));
		}
		
	
		// update node position
		for(size_t iP = 0; iP < NUM_NODES; ++iP) {
			// update a 2d vector: x_p += u_p * tau
			attachVector(coor_node[iP], velocity_node[iP], tau);
		}


		// update cell specific volume
		for(size_t iC = 0; iC < NUM_CELLS; ++iC) {
			double volume_flux = 0;
			// sum on cell vertexs
			for(size_t iPC = 0; iPC < NUM_CELL_VERTEXS; ++iPC) {
				// transform local node idx to global
				size_t idxP = localToGlobalNodeIdx(iC, iPC);
				double length_corner  = calCornerLength(iC, idxP, meshinfo);
				double normal_corner[2] = calCornerNormal(iC, idxP, meshinfo);
				
				sp_volume_flux += 
					length_corner * vectorMul2(normal_corner, velocity_node[idxP]);
			}
			sp_volume[iC] += tau * sp_volume_flux;
		}


		// update cell specific momentum
		for(size_t iC = 0; iC < NUM_CELLS; ++iC) {
			double resul_force = 0;
			for(size_t iPC = 0; iPC < NUM_CELL_VERTEX; ++iPC) {
				resul_force += calCornerForce(iC, iP, meshinfo, state);
			}
			attachVector(sp_momentum[iC], -resul_force, tau);
		}


		// update cell specific total energy
		for(size_t iC = 0; iC < NUM_CELLS; ++iC) {
			double expan_power = 0;
			for(size_t iPC = 0; iPC < NUM_CELL_VERTEX; ++iPC) {
				double corner_force = calCornerForce(iC, iP, meshinfo, state);
				size_t idxP = localToGlobalNodeIdx(iC, iPC);
				expan_power += vectorMul2(corner_force, velocity_node[idxP]);
			}
			sp_energy[iC] += - tau * expan_power;
		}

		
	}


}
