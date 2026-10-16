#include <stdio.h>

#include "linear_algebra_2d.h"
#include "mesh_calculation.h"
#include "physical_calculation.h"
#include "initial_bound_condition.h"

// 初始化, 可参考SGH
// 处理边界条件, 区分虚拟单元

// 结果输出

int main() 
{

/*******************
* Model Coefficients
*******************/

	const double GAMMA = 1.4;

	const double T = 1.0;

	const double CFL = 0.2;

	const double L0 = 1.0;
	const double L1 = 1.0;

	const size_t NUM_CELL_VERTEXS = 4;
	const size_t NUM_NODE_NEIGH_CELLS = 4;

	const size_t NUM_CELLS_X = 10;
	const size_t NUM_CELLS_Y = 10;
	const size_t NUM_NODES_X = NUM_CELLS_X + 1;
	const size_t NUM_NODES_Y = NUM_CELLS_Y + 1;
	const size_t NUM_CELLS = NUM_CELLS_X * NUM_CELLS_Y;
	const size_t NUM_NODES = NUM_NODES_X * NUM_NODES_Y;

/*******************
* Memory Allocation
*******************/
	double coor_node[NUM_NODES][2];
	double velocity_node[NUM_NODES][2];

	double mass_cell[NUM_CELLS];
	double density_cell[NUM_CELLS];
	double pressure_cell[NUM_CELLS];
	double volume_cell[NUM_CELLS];
	double sp_energy_cell[NUM_CELLS];
	double sp_momentum_cell[NUM_CELLS][2];



/*******************
* Initialization
*******************/

	// (uniform) Grid geometry
	double cell_h_0 = L0 / NUM_CELLS_X;
	double cell_h_1 = L1 / NUM_CELLS_Y;

	for(size_t iP = 0; iP < NUM_NODES; ++iP) {
		initializeNodeCoor(iP, NUM_NODES_X, cell_h_0, cell_h_1, coor_node[iP]);
	}

	
	// assign density, velocity, pressure, and get the others
	for(size_t iC = 0; iC < NUM_CELLS; ++iC) {

		volume_cell[iC]  = initializeCellVolume(cell_h_0, cell_h_1);
		density_cell[iC] = initializeCellDensity();
		mass_cell[iC] 	 = initializeCellMass(volume_cell[iC], density_cell[iC]);

		initializeCellSpMomentum(sp_momentum_cell[iC]);

		pressure_cell[iC] = initializeCellPresure(coor_node[2], GAMMA);

		sp_energy_cell[iC] = initializeCellSpEnergy(density_cell[iC],
								sp_momentum_cell[iC], pressure_cell[iC], GAMMA);
	}

	// we don't initialize veolicty_node


/********************
* Time Advancing
*********************/

	double t = 0.0;

	while (t < T) {
		// Apply boundary condition

		
		
		// Calculate time step length: tau
			double tau = 1e-4;


		// Update node velocity
		for(size_t iP = 0; iP < NUM_NODES; ++iP) {
			double pressure_gradient[2] = {0,};
			double ave_velocity_node[2] = {0,};

			double matrix_node[2][2] = {0,};
			
			for(size_t iCP = 0; iCP < NUM_NODE_NEIGH_CELLS; ++iCP) {
				// transform local cell idx to global
				size_t idxC = localToGlobalCellIdx(iP, iCP, NUM_NODES_X, NUM_CELLS_X);
				// calculate corner length and normal by neigh nodes

				size_t idxPl = calNeighVertexIdx(iP, idxC,  1, NUM_NODES_X, NUM_CELLS_X);
				size_t idxPr = calNeighVertexIdx(iP, idxC, -1, NUM_NODES_X, NUM_CELLS_X);

				double length_corner  = calCornerLength(coor_node[iP], 
											coor_node[idxPl], coor_node[idxPr]);
	
				double normal_corner[2];
				calCornerNormal(coor_node[iP], coor_node[idxPl], coor_node[idxPr],
								normal_corner);

				double impendance = calImpendance(density_cell[idxC], 
												  pressure_cell[idxC], GAMMA);
				double matrix_corner[2][2];
				calMatCorner_glace(length_corner, normal_corner, impendance,
									matrix_corner);

				double contrib_pressure_gradient[2];
				numMulVec2(length_corner * pressure_cell[idxC], normal_corner,
							contrib_pressure_gradient);
				attVec2(pressure_gradient, contrib_pressure_gradient);

				double contrib_ave_velocity_node[2];
				matMulVec2(matrix_corner, sp_momentum_cell[idxC],
							contrib_ave_velocity_node);
				attVec2(ave_velocity_node, contrib_ave_velocity_node);

				attMat2(matrix_node, matrix_corner);
			}

			
			double node_matrix_inverse[2][2];
			inverseMatrix2(matrix_node, node_matrix_inverse);

			double corrected_force[2];
			vecAddVec2(pressure_gradient, ave_velocity_node, corrected_force);
			matMulVec2(node_matrix_inverse, corrected_force, velocity_node[iP]);
		}
		
	
		// Update node position
		for(size_t iP = 0; iP < NUM_NODES; ++iP) {
			// update a 2d vector: x_p += u_p * tau
			double position_update[2];
			numMulVec2(tau, velocity_node[iP], position_update);
			attVec2(coor_node[iP], position_update);
		}


		// Update cell specific volume
		//		 	   specific momentum
		//			   specific total energy
		for(size_t iC = 0; iC < NUM_CELLS; ++iC) {
			double volume_flux = 0;
			double resul_force[2] = {0,};
			double expan_power = 0;

			
			// sum contributes on cell vertexs
			for(size_t iPC = 0; iPC < NUM_CELL_VERTEXS; ++iPC) {
				// transform local node idx to global
				size_t idxP = localToGlobalNodeIdx(iC, iPC, NUM_NODES_X, NUM_CELLS_X);
				size_t idxPl = calNeighVertexIdx(idxP, iC,  1, NUM_NODES_X, NUM_CELLS_X);
				size_t idxPr = calNeighVertexIdx(idxP, iC, -1, NUM_NODES_X, NUM_CELLS_X);

				double length_corner  = calCornerLength(coor_node[idxP], 
											coor_node[idxPl], coor_node[idxPr]);
														
				double normal_corner[2];
				calCornerNormal(coor_node[idxP], coor_node[idxPl], coor_node[idxPr], 
								normal_corner);

				double impendance = calImpendance(density_cell[iC], 
												pressure_cell[iC], GAMMA);
				double corner_force[2];
				calCornerForce(corner_force, impendance,
								sp_momentum_cell[iC], pressure_cell[iC], 
								length_corner, normal_corner,
								velocity_node[idxP]);

				/**********************************/

				double velocity_component = 
								vecMulVec2(normal_corner, velocity_node[idxP]);
				volume_flux += length_corner * velocity_component;

				attVec2(resul_force, corner_force);

				expan_power += vecMulVec2(corner_force, velocity_node[idxP]);
			}

			// update
			volume_cell[iC] += tau * volume_flux / mass_cell[iC];

			double contrib_sp_momentum[2];
			numMulVec2(-tau/mass_cell[iC], resul_force, contrib_sp_momentum);
			attVec2(sp_momentum_cell[iC], contrib_sp_momentum);

			sp_energy_cell[iC] += - tau * expan_power / mass_cell[iC];
		}
	}


}
