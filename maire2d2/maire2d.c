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
	const size_t END_TOTAL_CELL_X = NUM_GHOST_CELLS_X - 1;
	const size_t START_TOTAL_CELL_Y = 0;
	const size_t END_TOTAL_CELL_Y = NUM_GHOST_CELLS_Y - 1;
	const size_t START_TOTAL_NODE_X = 0;
	const size_t END_TOTAL_NODE_X = NUM_GHOST_NODES_X - 1;
	const size_t START_TOTAL_NODE_Y = 0;
	const size_t END_TOTAL_NODE_Y = NUM_GHOST_NODES_Y - 1;


/*******************
* Memory Allocation
*******************/
	double coor_node[NUM_GHOST_NODES_Y][NUM_GHOST_NODES_X][2];
	double velocity_node[NUM_GHOST_NODES_Y][NUM_GHOST_NODES_X][2];

	double mass_cell[NUM_GHOST_CELLS_Y][NUM_GHOST_CELLS_X];
	double density_cell[NUM_GHOST_CELLS_Y][NUM_GHOST_CELLS_X];
	double pressure_cell[NUM_GHOST_CELLS_Y][NUM_GHOST_CELLS_X];
	double volume_cell[NUM_GHOST_CELLS_Y][NUM_GHOST_CELLS_X];
	double sp_energy_cell[NUM_GHOST_CELLS_Y][NUM_GHOST_CELLS_X];
	double sp_momentum_cell[NUM_GHOST_CELLS_Y][NUM_GHOST_CELLS_X][2];



/*******************
* Initialization
*******************/

	// uniform grid geometry
	double cell_h_0 = L0 / NUM_CELLS_X;
	double cell_h_1 = L1 / NUM_CELLS_Y;

	for(size_t node_y = START_NODE_Y; node_y <= END_NODE_Y; ++node_y) {
		for(size_t node_x = START_NODE_X; node_x <= END_NODE_X; ++node_x) {
			initializeNodeCoor(node_x, node_y, cell_h_0, cell_h_1, coor_node[node_y][node_x]);

			printf("node (%zu, %zu) (%f, %f)\n", node_x, node_y, coor_node[node_y][node_x][0], coor_node[node_y][node_x][1]);
		}
	}
	// we don't need to initialize veolicty_node

	
	
	// assign density, velocity, pressure, and get the others
	for(size_t cell_y = START_CELL_Y; cell_y <= END_CELL_Y; ++cell_y) {
		for(size_t cell_x = START_CELL_X; cell_x <= END_CELL_X; ++cell_x) {
			volume_cell[cell_y][cell_x]  = initializeCellVolume(cell_h_0, cell_h_1);
			density_cell[cell_y][cell_x] = initializeCellDensity();
			mass_cell[cell_y][cell_x]    = initializeCellMass(volume_cell[cell_y][cell_x], density_cell[cell_y][cell_x]);

			initializeCellSpMomentum(sp_momentum_cell[cell_y][cell_x]);

			pressure_cell[cell_y][cell_x] = initializeCellPresure(coor_node[0][2], GAMMA);

			sp_energy_cell[cell_y][cell_x] = initializeCellSpEnergy(density_cell[cell_y][cell_x],
									sp_momentum_cell[cell_y][cell_x], pressure_cell[cell_y][cell_x], GAMMA);
		}
	}




/********************
* Time Advancing
*********************/

	double t = 0.0;

	size_t time_step = 0;

	while (t < T) {
		// Apply boundary condition
		for(size_t cell_y = START_TOTAL_CELL_Y; cell_y <= END_TOTAL_CELL_Y; ++cell_y) {
			for(size_t cell_x = START_TOTAL_CELL_X; cell_x <= END_TOTAL_CELL_X; ++cell_x) {
				if (cell_x >= START_CELL_X && cell_x <= END_CELL_X &&
				    cell_y >= START_CELL_Y && cell_y <= END_CELL_Y)
					continue;

			
				mass_cell[cell_y][cell_x] 	  = boundarilizeMass(cell_h_0, cell_h_1);
				density_cell[cell_y][cell_x]   = boundarilizeDensity();
				pressure_cell[cell_y][cell_x]  = boundarilizePressure();
				volume_cell[cell_y][cell_x] 	  = boundarilizeVolume(cell_h_0, cell_h_1);
				sp_energy_cell[cell_y][cell_x] = boundarilizeSpEnergy(GAMMA);

				boundarilizeSpMomentum(sp_momentum_cell[cell_y][cell_x]);
			}
		}

		for(size_t node_y = START_TOTAL_NODE_Y; node_y <= END_TOTAL_NODE_Y; ++node_y) {
			for(size_t node_x = START_TOTAL_NODE_X; node_x <= END_TOTAL_NODE_X; ++node_x) {
				if (node_x >= START_NODE_X && node_x <= END_NODE_X &&
				    node_y >= START_NODE_Y && node_y <= END_NODE_Y)
					continue;


				boundarilizeCoor(node_x, node_y, cell_h_0, cell_h_1, coor_node[node_y][node_x]);

			
			}
		}
		

		// Calculate time step length: tau
			double tau = 1e-3;
			t += tau;
			++time_step;


		// Update node velocity
		for(size_t node_y = START_NODE_Y; node_y <= END_NODE_Y; ++node_y) {
			for(size_t node_x = START_NODE_X; node_x <= END_NODE_X; ++node_x) {
				double pressure_gradient[2] = {0,};
				double ave_velocity_node[2] = {0,};

				double matrix_node[2][2] = {0,};
			
				for(size_t cell_y = node_y - 1; cell_y <= node_y; ++cell_y) {
					for(size_t cell_x = node_x - 1; cell_x <= node_x; ++cell_x) {
						size_t node_next_x, node_next_y;
						size_t node_prev_x, node_prev_y;
						calNeighVertexIdx(node_x, node_y, cell_x, cell_y, 1,
						                  &node_next_x, &node_next_y);
						calNeighVertexIdx(node_x, node_y, cell_x, cell_y, -1,
						                  &node_prev_x, &node_prev_y);

						double length_corner  = calCornerLength(coor_node[node_y][node_x],
													coor_node[node_next_y][node_next_x], coor_node[node_prev_y][node_prev_x]);
	
						double normal_corner[2];
						calCornerNormal(coor_node[node_y][node_x], coor_node[node_next_y][node_next_x], coor_node[node_prev_y][node_prev_x],
										normal_corner);

						double impendance = calImpendance(density_cell[cell_y][cell_x],
														  pressure_cell[cell_y][cell_x], GAMMA);
						double matrix_corner[2][2];
						calMatCorner_glace(length_corner, normal_corner, impendance,
											matrix_corner);

						double contrib_pressure_gradient[2];
						numMulVec2(length_corner * pressure_cell[cell_y][cell_x], normal_corner,
									contrib_pressure_gradient);
						attVec2(pressure_gradient, contrib_pressure_gradient);

						double contrib_ave_velocity_node[2];
						matMulVec2(matrix_corner, sp_momentum_cell[cell_y][cell_x],
									contrib_ave_velocity_node);
						attVec2(ave_velocity_node, contrib_ave_velocity_node);

						attMat2(matrix_node, matrix_corner);
					}
				}

			
				double node_matrix_inverse[2][2];
				inverseMatrix2(matrix_node, node_matrix_inverse);

				double corrected_force[2];
				vecAddVec2(pressure_gradient, ave_velocity_node, corrected_force);
				matMulVec2(node_matrix_inverse, corrected_force, velocity_node[node_y][node_x]);
			}
		}
		
	
		// Update node position
		for(size_t node_y = START_NODE_Y; node_y <= END_NODE_Y; ++node_y) {
			for(size_t node_x = START_NODE_X; node_x <= END_NODE_X; ++node_x) {
				// update a 2d vector: x_p += u_p * tau
				double position_update[2];
				numMulVec2(tau, velocity_node[node_y][node_x], position_update);
				attVec2(coor_node[node_y][node_x], position_update);
			}
		}


		// Update cell specific volume
		//		 	   specific momentum
		//			   specific total energy
		for(size_t cell_y = START_CELL_Y; cell_y <= END_CELL_Y; ++cell_y) {
			for(size_t cell_x = START_CELL_X; cell_x <= END_CELL_X; ++cell_x) {
				double volume_flux = 0;
				double resul_force[2] = {0,};
				double expan_power = 0;

				// sum contributes on cell vertexs
				for(size_t node_y = cell_y; node_y <= cell_y + 1; ++node_y) {
					for(size_t node_x = cell_x; node_x <= cell_x + 1; ++node_x) {
						size_t node_next_x, node_next_y;
						size_t node_prev_x, node_prev_y;
						calNeighVertexIdx(node_x, node_y, cell_x, cell_y, 1,
						                  &node_next_x, &node_next_y);
						calNeighVertexIdx(node_x, node_y, cell_x, cell_y, -1,
						                  &node_prev_x, &node_prev_y);

						double length_corner  = calCornerLength(coor_node[node_y][node_x],
													coor_node[node_next_y][node_next_x], coor_node[node_prev_y][node_prev_x]);
														
						double normal_corner[2];
						calCornerNormal(coor_node[node_y][node_x], coor_node[node_next_y][node_next_x], coor_node[node_prev_y][node_prev_x],
										normal_corner);

						double impendance = calImpendance(density_cell[cell_y][cell_x],
														pressure_cell[cell_y][cell_x], GAMMA);
						double corner_force[2];
						calCornerForce(corner_force, impendance,
										sp_momentum_cell[cell_y][cell_x], pressure_cell[cell_y][cell_x],
										length_corner, normal_corner,
										velocity_node[node_y][node_x]);

						/**********************************/

						double velocity_component =
										vecMulVec2(normal_corner, velocity_node[node_y][node_x]);
						volume_flux += length_corner * velocity_component;

						attVec2(resul_force, corner_force);

						expan_power += vecMulVec2(corner_force, velocity_node[node_y][node_x]);
					}
				}

				// update
				volume_cell[cell_y][cell_x] += tau * volume_flux / mass_cell[cell_y][cell_x];

				double contrib_sp_momentum[2];
				numMulVec2(-tau/mass_cell[cell_y][cell_x], resul_force, contrib_sp_momentum);
				attVec2(sp_momentum_cell[cell_y][cell_x], contrib_sp_momentum);

				sp_energy_cell[cell_y][cell_x] += - tau * expan_power / mass_cell[cell_y][cell_x];
			}
		}

		if (time_step / OUTPUT_FREQUNCY == 0)
			putOutNodeCoor(t, NUM_GHOST_NODES_X, NUM_GHOST_NODES_Y,
						   START_NODE_X, END_NODE_X, START_NODE_Y, END_NODE_Y, coor_node);
		//起点不对
	}


}
