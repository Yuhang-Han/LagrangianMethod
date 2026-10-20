#include <stdio.h>

#include "linear_algebra_2d.h"
#include "mesh_calculation.h"
#include "physical_calculation.h"
#include "initial_bound_condition.h"
#include "data_output.h"


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

	const size_t NUM_CELLS_X = 5;
	const size_t NUM_CELLS_Y = 5;
	const size_t NUM_CELLS = NUM_CELLS_X * NUM_CELLS_Y;
	
	const size_t NUM_NODES_X = NUM_CELLS_X + 1;
	const size_t NUM_NODES_Y = NUM_CELLS_Y + 1;
	const size_t NUM_NODES = NUM_NODES_X * NUM_NODES_Y;

	const size_t NUM_TOTAL_CELLS_X = NUM_CELLS_X + 2;
	const size_t NUM_TOTAL_CELLS_Y = NUM_CELLS_Y + 2;
	const size_t NUM_TOTAL_CELLS = NUM_TOTAL_CELLS_X * NUM_TOTAL_CELLS_Y;

	const size_t NUM_TOTAL_NODES_X = NUM_TOTAL_CELLS_X + 1;
	const size_t NUM_TOTAL_NODES_Y = NUM_TOTAL_CELLS_Y + 1;
	const size_t NUM_TOTAL_NODES = NUM_TOTAL_NODES_X * NUM_TOTAL_NODES_Y;

	const size_t NUM_BOUND_CELLS = NUM_TOTAL_CELLS - NUM_CELLS;
	const size_t NUM_BOUND_NODES = NUM_TOTAL_NODES - NUM_NODES;

printf("\n Computation coefficients:\n");
printf("cells_x: %zu, cells_y: %zu, cells: %zu\n", 
				NUM_CELLS_X, NUM_CELLS_Y, NUM_CELLS);
printf("nodes_x: %zu, nodes_y: %zu, nodes: %zu\n",
				NUM_NODES_X, NUM_NODES_Y, NUM_NODES);
printf("ghost cells_x: %zu, cells_y: %zu, cells: %zu\n", 
				NUM_TOTAL_CELLS_X, NUM_TOTAL_CELLS_Y, NUM_TOTAL_CELLS);
printf("ghost nodes_x: %zu, nodes_y: %zu, nodes: %zu\n",
				NUM_TOTAL_NODES_X, NUM_TOTAL_NODES_Y, NUM_TOTAL_NODES);
printf("bound cells: %zu\n", NUM_BOUND_CELLS);
printf("bound nodes: %zu\n", NUM_BOUND_NODES);


	// Inclusive physical and allocated index ranges.
	const size_t START_CELL_X = 1;
	const size_t END_CELL_X = NUM_CELLS_X;
	const size_t START_CELL_Y = 1;
	const size_t END_CELL_Y = NUM_CELLS_Y;
	const size_t START_NODE_X = 1;
	const size_t END_NODE_X = NUM_NODES_X;
	const size_t START_NODE_Y = 1;
	const size_t END_NODE_Y = NUM_NODES_Y;

	const size_t START_TOTAL_CELL_X = 0;
	const size_t END_TOTAL_CELL_X = NUM_TOTAL_CELLS_X - 1;
	const size_t START_TOTAL_CELL_Y = 0;
	const size_t END_TOTAL_CELL_Y = NUM_TOTAL_CELLS_Y - 1;
	const size_t START_TOTAL_NODE_X = 0;
	const size_t END_TOTAL_NODE_X = NUM_TOTAL_NODES_X - 1;
	const size_t START_TOTAL_NODE_Y = 0;
	const size_t END_TOTAL_NODE_Y = NUM_TOTAL_NODES_Y - 1;


/*******************
* Memory Allocation
*******************/
	double coor_node[NUM_TOTAL_NODES_Y][NUM_TOTAL_NODES_X][2];
	double velocity_node[NUM_TOTAL_NODES_Y][NUM_TOTAL_NODES_X][2];

	double mass_cell[NUM_TOTAL_CELLS_Y][NUM_TOTAL_CELLS_X];
	double density_cell[NUM_TOTAL_CELLS_Y][NUM_TOTAL_CELLS_X];
	double pressure_cell[NUM_TOTAL_CELLS_Y][NUM_TOTAL_CELLS_X];
	double volume_cell[NUM_TOTAL_CELLS_Y][NUM_TOTAL_CELLS_X];
	double sp_energy_cell[NUM_TOTAL_CELLS_Y][NUM_TOTAL_CELLS_X];
	double sp_momentum_cell[NUM_TOTAL_CELLS_Y][NUM_TOTAL_CELLS_X][2];



/*******************
* Initialization
*******************/

	// uniform grid geometry
	double cell_h_0 = L0 / NUM_CELLS_X;
	double cell_h_1 = L1 / NUM_CELLS_Y;

	for(size_t iPy = START_NODE_Y; iPy <= END_NODE_Y; ++iPy) {
		for(size_t iPx = START_NODE_X; iPx <= END_NODE_X; ++iPx) {
			initializeNodeCoor(iPx, iPy, cell_h_0, cell_h_1, coor_node[iPy][iPx]);

			printf("node (%zu, %zu) (%f, %f)\n", iPx, iPy, coor_node[iPy][iPx][0], coor_node[iPy][iPx][1]);
		}
	}
	// we don't need to initialize veolicty_node

	
	
	// assign density, velocity, pressure, and get the others
	for(size_t iCy = START_CELL_Y; iCy <= END_CELL_Y; ++iCy) {
		for(size_t iCx = START_CELL_X; iCx <= END_CELL_X; ++iCx) {
			volume_cell[iCy][iCx]  = initializeCellVolume(cell_h_0, cell_h_1);
			density_cell[iCy][iCx] = initializeCellDensity();
			mass_cell[iCy][iCx]    = initializeCellMass(volume_cell[iCy][iCx], density_cell[iCy][iCx]);

			initializeCellSpMomentum(sp_momentum_cell[iCy][iCx]);

			pressure_cell[iCy][iCx] = initializeCellPresure(coor_node[0][2], GAMMA);

			sp_energy_cell[iCy][iCx] = initializeCellSpEnergy(density_cell[iCy][iCx],
									sp_momentum_cell[iCy][iCx], pressure_cell[iCy][iCx], GAMMA);
		}
	}




/********************
* Time Advancing
*********************/

	double t = 0.0;

	size_t time_step = 0;

	while (t < T) {
		// Apply boundary condition
		for(size_t iCy = START_TOTAL_CELL_Y; iCy <= END_TOTAL_CELL_Y; ++iCy) {
			for(size_t iCx = START_TOTAL_CELL_X; iCx <= END_TOTAL_CELL_X; ++iCx) {
				if (iCx >= START_CELL_X && iCx <= END_CELL_X &&
				    iCy >= START_CELL_Y && iCy <= END_CELL_Y)
					continue;

			
				mass_cell[iCy][iCx] 	 = boundarilizeMass(cell_h_0, cell_h_1);
				density_cell[iCy][iCx]   = boundarilizeDensity();
				pressure_cell[iCy][iCx]  = boundarilizePressure();
				volume_cell[iCy][iCx]    = boundarilizeVolume(cell_h_0, cell_h_1);
				sp_energy_cell[iCy][iCx] = boundarilizeSpEnergy(GAMMA);

				boundarilizeSpMomentum(sp_momentum_cell[iCy][iCx]);
			}
		}

		for(size_t iPy = START_TOTAL_NODE_Y; iPy <= END_TOTAL_NODE_Y; ++iPy) {
			for(size_t iPx = START_TOTAL_NODE_X; iPx <= END_TOTAL_NODE_X; ++iPx) {
				if (iPx >= START_NODE_X && iPx <= END_NODE_X &&
				    iPy >= START_NODE_Y && iPy <= END_NODE_Y)
					continue;

				boundarilizeCoor(iPx, iPy, cell_h_0, cell_h_1, 
								 coor_node[iPy][iPx]);
			}
		}
		

		// Calculate time step length: tau
			double tau = 1e-3;
			t += tau;
			++time_step;


		// Update node velocity
		for(size_t iPy = START_NODE_Y; iPy <= END_NODE_Y; ++iPy) {
			for(size_t iPx = START_NODE_X; iPx <= END_NODE_X; ++iPx) {
				double pressure_gradient[2] = {0,};
				double ave_velocity_node[2] = {0,};

				double matrix_node[2][2] = {0,};
			
				for(size_t iCy = iPy - 1; iCy <= iPy; ++iCy) {
					for(size_t iCx = iPx - 1; iCx <= iPx; ++iCx) {
						size_t node_next_x, node_next_y;
						size_t node_prev_x, node_prev_y;
						calNeighVertexIdx(iPx, iPy, iCx, iCy, 1,
						                  &node_next_x, &node_next_y);
						calNeighVertexIdx(iPx, iPy, iCx, iCy, -1,
						                  &node_prev_x, &node_prev_y);

						double length_corner = \
							calCornerLength(coor_node[iPy][iPx],
											coor_node[node_next_y][node_next_x], 
											coor_node[node_prev_y][node_prev_x]);
	
						double normal_corner[2];
						calCornerNormal(coor_node[iPy][iPx], 
										coor_node[node_next_y][node_next_x], 
										coor_node[node_prev_y][node_prev_x],
										normal_corner);

						double impendance = \
							calImpendance(density_cell[iCy][iCx],
										  pressure_cell[iCy][iCx], GAMMA);
						double matrix_corner[2][2];
						calMatCorner_glace(length_corner, normal_corner, 
											impendance, matrix_corner);

						double contrib_pressure_gradient[2];
						numMulVec2(length_corner * pressure_cell[iCy][iCx], 
									normal_corner, contrib_pressure_gradient);
						attVec2(pressure_gradient, contrib_pressure_gradient);

						double contrib_ave_velocity_node[2];
						matMulVec2(matrix_corner, sp_momentum_cell[iCy][iCx],
									contrib_ave_velocity_node);
						attVec2(ave_velocity_node, contrib_ave_velocity_node);

						attMat2(matrix_node, matrix_corner);
					}
				}

			
				double node_matrix_inverse[2][2];
				inverseMatrix2(matrix_node, node_matrix_inverse);

				double corrected_force[2];
				vecAddVec2(pressure_gradient, ave_velocity_node, corrected_force);
				matMulVec2(node_matrix_inverse, corrected_force, 
							velocity_node[iPy][iPx]);
			}
		}
		
	
		// Update node position
		for(size_t iPy = START_NODE_Y; iPy <= END_NODE_Y; ++iPy) {
			for(size_t iPx = START_NODE_X; iPx <= END_NODE_X; ++iPx) {
				// update a 2d vector: x_p += u_p * tau
				double position_update[2];
				numMulVec2(tau, velocity_node[iPy][iPx], position_update);
				attVec2(coor_node[iPy][iPx], position_update);
			}
		}


		// Update cell volume
		//		 	   specific momentum
		//			   specific total energy
		for(size_t iCy = START_CELL_Y; iCy <= END_CELL_Y; ++iCy) {
			for(size_t iCx = START_CELL_X; iCx <= END_CELL_X; ++iCx) {
				double volume_flux = 0;
				double resul_force[2] = {0,};
				double expan_power = 0;

				// sum contributes on cell vertexs
				for(size_t iPy = iCy; iPy <= iCy + 1; ++iPy) {
					for(size_t iPx = iCx; iPx <= iCx + 1; ++iPx) {
						size_t node_next_x, node_next_y;
						size_t node_prev_x, node_prev_y;
						calNeighVertexIdx(iPx, iPy, iCx, iCy, 1,
						                  &node_next_x, &node_next_y);
						calNeighVertexIdx(iPx, iPy, iCx, iCy, -1,
						                  &node_prev_x, &node_prev_y);

						double length_corner  = \
								calCornerLength(coor_node[iPy][iPx],
												coor_node[node_next_y][node_next_x], 
												coor_node[node_prev_y][node_prev_x]);
														
						double normal_corner[2];
						calCornerNormal(coor_node[iPy][iPx], 
										coor_node[node_next_y][node_next_x], 
										coor_node[node_prev_y][node_prev_x],
										normal_corner);

						double impendance = 
								calImpendance(density_cell[iCy][iCx],
											  pressure_cell[iCy][iCx], 
											  GAMMA);
						double corner_force[2];
						calCornerForce(corner_force, impendance,
										sp_momentum_cell[iCy][iCx], 
										pressure_cell[iCy][iCx],
										length_corner, normal_corner,
										velocity_node[iPy][iPx]);

						/**********************************/

						double velocity_component =
									vecMulVec2(normal_corner, 
												velocity_node[iPy][iPx]);
						volume_flux += length_corner * velocity_component;

						attVec2(resul_force, corner_force);

						expan_power += vecMulVec2(corner_force, 
												velocity_node[iPy][iPx]);
					}
				}

				// update
				volume_cell[iCy][iCx] += tau * volume_flux \
												/ mass_cell[iCy][iCx];

				double contrib_sp_momentum[2];
				numMulVec2(-tau/mass_cell[iCy][iCx], 
							resul_force, contrib_sp_momentum);
				attVec2(sp_momentum_cell[iCy][iCx], contrib_sp_momentum);

				sp_energy_cell[iCy][iCx] += - tau * expan_power \
												/ mass_cell[iCy][iCx];
			}
		}

		if (time_step / OUTPUT_FREQUNCY == 0)
			putOutNodeCoor(t, NUM_TOTAL_NODES_X, NUM_TOTAL_NODES_Y,
						   START_NODE_X, END_NODE_X, START_NODE_Y, END_NODE_Y, 
						   coor_node);
		//起点不对
	}


}
