#include <stdio.h>

#include "linear_algebra_2d.h"
#include "mesh_calculation.h"
#include "physical_calculation.h"

// 先写函数再声明变量
// 矩阵换成二维数组, 不用指针

// 顺一遍主函数, 匹配函数定义, 补充没定义的函数

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

	const size_t NUM_CELLS_X = 10;
	const size_t NUM_NODES_X = NUM_CELLS_X + 1;
	const size_t NUM_CELLS = NUM_CELLS_X * NUM_CELLS_X;
	const size_t NUM_NODES = NUM_NODES_X * NUM_NODES_X;

/*******************
* Memory Allocation
*******************/




/*******************
* Initialization
*******************/



/********************
* Time Advancing
*********************/



// Update
	while (t < T) {
		// apply boundary condition
		
		// Calculate time step length: tau
			double tau = 1e-4;


		// update node velocity
		for(size_t iP = 0; iP < NUM_NODES; ++iP) {
			double pressure_gradient[2] = {0,};
			double ave_velocity_node[2] = {0,};

			double matrix_node[2][2] = {0,};
			
			for(size_t iCP = 0; iCP < NUM_NODE_NEIGH_CELLS; ++iCP) {
				// transform local cell idx to global
				size_t idxC = localToGlobalCellIdx(iP, iCP,
									NUM_NODES_X, NUM_CELLS_X);
				// calculate corner length and normal by neigh nodes

				size_t idxPl = calNeighVertexIdx(iP, idxC,  1, 
												NUM_NODES_X, NUM_CELLS_X);
				size_t idxPr = calNeighVertexIdx(iP, idxC, -1, 
												NUM_NODES_X, NUM_CELLS_X);
				double length_corner  = calCornerLength(coor_node[iP], 
											coor_node[idxPl], coor_node[idxPr]);
	
				double normal_corner[2];
				calCornerNormal(coor_node[iP], coor_pl, coor_pr, normal_corner);

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
		
	
		// update node position
		for(size_t iP = 0; iP < NUM_NODES; ++iP) {
			// update a 2d vector: x_p += u_p * tau
			double position_update[2];
			numMulVec2(tau, velocity_node[iP], position_update);
			attVec2(coor_node[iP], position_update);
		}


		// update cell specific volume
		for(size_t iC = 0; iC < NUM_CELLS; ++iC) {
			double sp_volume_flux = 0;
			// sum on cell vertexs
			for(size_t iPC = 0; iPC < NUM_CELL_VERTEXS; ++iPC) {
				// transform local node idx to global
				size_t idxP = localToGlobalNodeIdx(iC, iPC, 
										NUM_NODES_X, NUM_CELLS_X);
										
				size_t idxPl = calNeighVertexIdx(idxP, iC,  1, 
												NUM_NODES_X, NUM_CELLS_X);
				size_t idxPr = calNeighVertexIdx(idxP, iC, -1, 
												NUM_NODES_X, NUM_CELLS_X);
				double length_corner  = calCornerLength(coor_node[idxP], 
											coor_node[idxPl], coor_node[idxPr]);
														
				double normal_corner[2];
				calCornerNormal(coor_node[idxP], coor_pl, coor_pr, normal_corner);

				double velocity_component = 
						vecMulVec2(normal_corner, velocity_node[idxP]);
				sp_volume_flux += length_corner * velocity_component;
			}
			sp_volume_cell[iC] += tau * sp_volume_flux;
		}


		// update cell specific momentum
		for(size_t iC = 0; iC < NUM_CELLS; ++iC) {
			double resul_force[2] = {0,};
			for(size_t iPC = 0; iPC < NUM_CELL_VERTEXS; ++iPC) {

				size_t idxP = localToGlobalNodeIdx(iC, iPC, 
										NUM_NODES_X, NUM_CELLS_X);
										
				size_t idxPl = calNeighVertexIdx(idxP, iC,  1, 
												NUM_NODES_X, NUM_CELLS_X);
				size_t idxPr = calNeighVertexIdx(idxP, iC, -1, 
												NUM_NODES_X, NUM_CELLS_X);
				double length_corner  = calCornerLength(coor_node[idxP], 
											coor_node[idxPl], coor_node[idxPr]);
														
				double normal_corner[2];
				calCornerNormal(coor_node[idxP], coor_pl, coor_pr, normal_corner);

				double impendance = calImpendance(density_cell[iC], 
												pressure_cell[iC], GAMMA);
				double corner_force[2];
				calCornerForce(corner_force, sp_momentum_cell[iC], 
								pressure_cell[iC], impendance,
								length_corner, normal_corner,
								velociy_node[idxP]);
									
				attVec2(resul_force, corner_force);
				
			}

			double contrib_sp_momentum[2];
			numMulVec2(-tau, resul_force, contrib_sp_momentum);
			attVec2(sp_momentum_cell[iC], contrib_sp_momentum);
		}


		// update cell specific total energy
		for(size_t iC = 0; iC < NUM_CELLS; ++iC) {
			double expan_power = 0;
			for(size_t iPC = 0; iPC < NUM_CELL_VERTEXS; ++iPC) {
				size_t idxP = localToGlobalNodeIdx(iC, iPC, 
										NUM_NODES_X, NUM_CELLS_X);
										
				size_t idxPl = calNeighVertexIdx(idxP, iC,  1, 
												NUM_NODES_X, NUM_CELLS_X);
				size_t idxPr = calNeighVertexIdx(idxP, iC, -1, 
												NUM_NODES_X, NUM_CELLS_X);
				double length_corner  = calCornerLength(coor_node[idxP], 
											coor_node[idxPl], coor_node[idxPr]);
														
				double normal_corner[2];
				calCornerNormal(coor_node[idxP], coor_pl, coor_pr, normal_corner);

				double impendance = calImpendance(density_cell[iC], 
												pressure_cell[iC], GAMMA);
				double corner_force[2];
				calCornerForce(corner_force, sp_momentum_cell[iC], 
								pressure_cell[iC], impendance,
								length_corner, normal_corner,
								velociy_node[idxP]);
				
				expan_power += vecMulVec2(corner_force, velocity_node[idxP]);
			}
			sp_energy_cell[iC] += - tau * expan_power;
		}

		
	}


}
