#include <stdio.h>

#include "linear_algebra_2d.h"
#include "mesh_calculation.h"
#include "physical_calculation.h"
#include "initial_bound_condition.h"
#include "data_output.h"

// 顶不住了, 还是用二维编号吧

int main() 
{

/*******************
* Model Coefficients
*******************/

	const size_t OUTPUT_FREQUNCY = 100;

	const double GAMMA = 1.4;

	const double T = 1.0;

	const double CFL = 0.2;

	const double L0 = 1.0;
	const double L1 = 1.0;

	const size_t NUM_CELL_VERTEXS = 4;
	const size_t NUM_NODE_NEIGH_CELLS = 4;

	const size_t NUM_CELLS_X = 5;
	const size_t NUM_CELLS_Y = 5;
	const size_t NUM_CELLS = NUM_CELLS_X * NUM_CELLS_Y;
	
	const size_t NUM_NODES_X = NUM_CELLS_X + 1;
	const size_t NUM_NODES_Y = NUM_CELLS_Y + 1;
	const size_t NUM_NODES = NUM_NODES_X * NUM_NODES_Y;

	const size_t NUM_GHOST_CELLS_X = NUM_CELLS_X + 2;
	const size_t NUM_GHOST_CELLS_Y = NUM_CELLS_Y + 2;
	const size_t NUM_GHOST_CELLS = NUM_GHOST_CELLS_X * NUM_GHOST_CELLS_Y;

	const size_t NUM_GHOST_NODES_X = NUM_GHOST_CELLS_X + 1;
	const size_t NUM_GHOST_NODES_Y = NUM_GHOST_CELLS_Y + 1;
	const size_t NUM_GHOST_NODES = NUM_GHOST_NODES_X * NUM_GHOST_NODES_Y;

	const size_t NUM_BOUND_CELLS = NUM_GHOST_CELLS - NUM_CELLS;
	const size_t NUM_BOUND_NODES = NUM_GHOST_NODES - NUM_NODES;

printf("\n Computation coefficients:\n");
printf("cells_x: %zu, cells_y: %zu, cells: %zu\n", 
				NUM_CELLS_X, NUM_CELLS_Y, NUM_CELLS);
printf("nodes_x: %zu, nodes_y: %zu, nodes: %zu\n",
				NUM_NODES_X, NUM_NODES_Y, NUM_NODES);
printf("ghost cells_x: %zu, cells_y: %zu, cells: %zu\n", 
				NUM_GHOST_CELLS_X, NUM_GHOST_CELLS_Y, NUM_GHOST_CELLS);
printf("ghost nodes_x: %zu, nodes_y: %zu, nodes: %zu\n",
				NUM_GHOST_NODES_X, NUM_GHOST_NODES_Y, NUM_GHOST_NODES);
printf("bound cells: %zu\n", NUM_BOUND_CELLS);
printf("bound nodes: %zu\n", NUM_BOUND_NODES);


	// make physical and ghost node/cell dictionary
	size_t NODE_DIC[NUM_NODES];
	size_t CELL_DIC[NUM_CELLS];
	size_t BOUND_CELL_DIC[NUM_BOUND_CELLS];
	size_t BOUND_NODE_DIC[NUM_BOUND_NODES];

	categorizeNodes(NUM_GHOST_NODES_X, NUM_GHOST_NODES_Y, NUM_NODES, NUM_BOUND_NODES,
					NODE_DIC, BOUND_NODE_DIC);
	categorizeCells(NUM_GHOST_CELLS_X, NUM_GHOST_CELLS_Y, NUM_CELLS, NUM_BOUND_CELLS,
					CELL_DIC, BOUND_CELL_DIC);
	

/*******************
* Memory Allocation
*******************/
	double coor_node[NUM_GHOST_NODES][2];
	double velocity_node[NUM_GHOST_NODES][2];

	double mass_cell[NUM_GHOST_CELLS];
	double density_cell[NUM_GHOST_CELLS];
	double pressure_cell[NUM_GHOST_CELLS];
	double volume_cell[NUM_GHOST_CELLS];
	double sp_energy_cell[NUM_GHOST_CELLS];
	double sp_momentum_cell[NUM_GHOST_CELLS][2];



/*******************
* Initialization
*******************/

	// uniform grid geometry
	double cell_h_0 = L0 / NUM_CELLS_X;
	double cell_h_1 = L1 / NUM_CELLS_Y;

	for(size_t iP = 0; iP < NUM_NODES; ++iP) {
		size_t idxP = NODE_DIC[iP];
		initializeNodeCoor(idxP, NUM_GHOST_NODES_X, cell_h_0, cell_h_1, coor_node[idxP]);

		printf("node %zu (%f, %f)\n", idxP, coor_node[idxP][0], coor_node[idxP][1]);
	}
	// we don't need to initialize veolicty_node

	
	
	// assign density, velocity, pressure, and get the others
	for(size_t iC = 0; iC < NUM_CELLS; ++iC) {
		size_t idxC = CELL_DIC[iC];
		volume_cell[idxC]  = initializeCellVolume(cell_h_0, cell_h_1);
		density_cell[idxC] = initializeCellDensity();
		mass_cell[idxC]    = initializeCellMass(volume_cell[idxC], density_cell[idxC]);

		initializeCellSpMomentum(sp_momentum_cell[idxC]);

		pressure_cell[idxC] = initializeCellPresure(coor_node[2], GAMMA);

		sp_energy_cell[idxC] = initializeCellSpEnergy(density_cell[idxC],
								sp_momentum_cell[idxC], pressure_cell[idxC], GAMMA);
	}




/********************
* Time Advancing
*********************/

	double t = 0.0;

	size_t time_step = 0;

	while (t < T) {
		// Apply boundary condition
		for(size_t iBC = 0;  iBC < NUM_BOUND_CELLS; ++iBC) {
			size_t idxBC = BOUND_CELL_DIC[iBC];
			
			mass_cell[idxBC] 	  = boundarilizeMass(cell_h_0, cell_h_1);
			density_cell[idxBC]   = boundarilizeDensity();
			pressure_cell[idxBC]  = boundarilizePressure();
			volume_cell[idxBC] 	  = boundarilizeVolume(cell_h_0, cell_h_1);
			sp_energy_cell[idxBC] = boundarilizeSpEnergy(GAMMA);

			boundarilizeSpMomentum(sp_momentum_cell[idxBC]);
		}

		for(size_t iBP = 0; iBP < NUM_BOUND_NODES; ++iBP) {
			size_t idxBP = BOUND_NODE_DIC[iBP];

			boundarilizeCoor(idxBP, NUM_GHOST_NODES_X, cell_h_0, cell_h_1, coor_node[idxBP]);

			
		}
		

		// Calculate time step length: tau
			double tau = 1e-3;
			t += tau;
			++time_step;


		// Update node velocity
		for(size_t iP = 0; iP < NUM_NODES; ++iP) {
			size_t idxP = NODE_DIC[iP];
			double pressure_gradient[2] = {0,};
			double ave_velocity_node[2] = {0,};

			double matrix_node[2][2] = {0,};
			
			for(size_t iCP = 0; iCP < NUM_NODE_NEIGH_CELLS; ++iCP) {
				// transform local cell idx to global
				size_t idxC = localToGlobalCellIdx(idxP, iCP, 
									NUM_GHOST_NODES_X, NUM_GHOST_CELLS_X);
				//printf("trying to touch %zu th cell\n", idxC);
				// calculate corner length and normal by neigh nodes

				size_t idxPl = calNeighVertexIdx(idxP, idxC,  1, 
												NUM_GHOST_NODES_X, NUM_GHOST_CELLS_X);
				size_t idxPr = calNeighVertexIdx(idxP, idxC, -1, 
												NUM_GHOST_NODES_X, NUM_GHOST_CELLS_X);

	
				double length_corner  = calCornerLength(coor_node[idxP], 
											coor_node[idxPl], coor_node[idxPr]);
	
				double normal_corner[2];
				calCornerNormal(coor_node[idxP], coor_node[idxPl], coor_node[idxPr],
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
			matMulVec2(node_matrix_inverse, corrected_force, velocity_node[idxP]);
		}
		
	
		// Update node position
		for(size_t iP = 0; iP < NUM_NODES; ++iP) {
			size_t idxP = NODE_DIC[iP];
			// update a 2d vector: x_p += u_p * tau
			double position_update[2];
			numMulVec2(tau, velocity_node[idxP], position_update);
			attVec2(coor_node[idxP], position_update);
		}


		// Update cell specific volume
		//		 	   specific momentum
		//			   specific total energy
		for(size_t iC = 0; iC < NUM_CELLS; ++iC) {
			size_t idxC = CELL_DIC[iC];
			double volume_flux = 0;
			double resul_force[2] = {0,};
			double expan_power = 0;

			// sum contributes on cell vertexs
			for(size_t iPC = 0; iPC < NUM_CELL_VERTEXS; ++iPC) {
				// transform local node idx to global
				size_t idxP = localToGlobalNodeIdx(idxC, iPC, 
												NUM_GHOST_NODES_X, NUM_GHOST_CELLS_X);
				size_t idxPl = calNeighVertexIdx(idxP, idxC,  1, 
												NUM_GHOST_NODES_X, NUM_GHOST_CELLS_X);
				size_t idxPr = calNeighVertexIdx(idxP, idxC, -1, 
												NUM_GHOST_NODES_X, NUM_GHOST_CELLS_X);


				double length_corner  = calCornerLength(coor_node[idxP], 
											coor_node[idxPl], coor_node[idxPr]);
														
				double normal_corner[2];
				calCornerNormal(coor_node[idxP], coor_node[idxPl], coor_node[idxPr], 
								normal_corner);

				double impendance = calImpendance(density_cell[idxC], 
												pressure_cell[idxC], GAMMA);
				double corner_force[2];
				calCornerForce(corner_force, impendance,
								sp_momentum_cell[idxC], pressure_cell[idxC], 
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
			volume_cell[idxC] += tau * volume_flux / mass_cell[idxC];

			double contrib_sp_momentum[2];
			numMulVec2(-tau/mass_cell[idxC], resul_force, contrib_sp_momentum);
			attVec2(sp_momentum_cell[idxC], contrib_sp_momentum);

			sp_energy_cell[idxC] += - tau * expan_power / mass_cell[idxC];
		}

		if (time_step / OUTPUT_FREQUNCY == 0)
			putOutNodeCoor(t, NUM_NODES_X, NUM_NODES_Y, NUM_NODES, coor_node);
		//起点不对
	}


}
