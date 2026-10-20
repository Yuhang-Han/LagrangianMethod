
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "linear_algebra_2d.h"
#include "mesh_calculation.h"
#include "physical_calculation.h"
#include "initial_bound_condition.h"
#include "data_output.h"

double getDoubleMin(double a, double b) {
	if(a < b)
		return a;
	else 
		return b;
}


int main() 
{

/*******************
* Model Coefficients
*******************/
	const double EPS = 1e-10;

	const size_t OUTPUT_FREQUNCY = 10;

	const double GAMMA = SEDOV_GAMMA;

	const double T = SEDOV_FINAL_TIME;

	const double CFL_E = 0.25;
	const double CFL_V = 0.1;

	const double L0 = SEDOV_DOMAIN_LENGTH;
	const double L1 = SEDOV_DOMAIN_LENGTH;

	const size_t NUM_CELLS_X = SEDOV_NUM_CELLS;
	const size_t NUM_CELLS_Y = SEDOV_NUM_CELLS;
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
	printf("\n Initializing...\n");

	// uniform grid geometry
	double cell_h_0 = L0 / NUM_CELLS_X;
	double cell_h_1 = L1 / NUM_CELLS_Y;


	for(size_t iPy = START_NODE_Y; iPy <= END_NODE_Y; ++iPy) {
		for(size_t iPx = START_NODE_X; iPx <= END_NODE_X; ++iPx) {
			initializeNodeCoor(iPx, iPy, cell_h_0, cell_h_1, 
								START_NODE_X, START_NODE_Y, 
								coor_node[iPy][iPx]);

			//printf("node (%zu, %zu) (%f, %f)\n", iPx, iPy, coor_node[iPy][iPx][0], coor_node[iPy][iPx][1]);
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

			pressure_cell[iCy][iCx] = \
				initializeCellPresure(iCx, iCy, cell_h_0, cell_h_1,
								      NUM_TOTAL_NODES_X, NUM_TOTAL_NODES_Y,
								      coor_node, GAMMA);
			
			sp_energy_cell[iCy][iCx] = initializeCellSpEnergy(density_cell[iCy][iCx],
									sp_momentum_cell[iCy][iCx], pressure_cell[iCy][iCx], GAMMA);
		}
	}


	printf("Initialization is done\n");



/********************
* Time Advancing
*********************/

	printf("\n Time advancing...\n");

	double t = 0.0;

	size_t time_step = 0;

	while (t < T) {
	// Apply boundary condition
		applySedovCellBoundaries(NUM_CELLS_X, NUM_CELLS_Y,
				cell_h_0, cell_h_1, GAMMA, mass_cell, density_cell,
				pressure_cell, volume_cell, sp_energy_cell, sp_momentum_cell);
		applySedovNodeBoundaries(NUM_CELLS_X, NUM_CELLS_Y, coor_node);

		
		
	// Calculate node velocity
		//printf("calculating node velocity...\n");
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
				applySedovVelocityBoundary(iPx, iPy, START_NODE_X, START_NODE_Y,
										   velocity_node[iPy][iPx]);
			}
		}

		//printf("calculating time step length...\n");
	// Time step criterion on the volume variation
		double tau_e = T - t;
		double tau_v = T - t;

		/* Only physical cells have all four nodal velocities computed above. */
		for(size_t iCy = START_CELL_Y; iCy <= END_CELL_Y; ++iCy) {
			for(size_t iCx = START_CELL_X; iCx <= END_CELL_X; ++iCx) {

				double impendance = calImpendance(density_cell[iCy][iCx],
												  pressure_cell[iCy][iCx],
												  GAMMA);

				double width_cell_approx = L0 + L1;
				double volume_flux = 0.;
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

						width_cell_approx = getDoubleMin(length_corner, 
														 width_cell_approx);
					}
				}
				double sound_speed = impendance / density_cell[iCy][iCx];
				double tau_e_cell = CFL_E * width_cell_approx / sound_speed;
				tau_e = getDoubleMin(tau_e_cell, tau_e);
				
				
				double tau_v_cell = CFL_V * volume_cell[iCy][iCx] 
									/ (fabs(volume_flux) + EPS);
				/*					
				printf("tau_v, tau_v_cell: %f, %f\n", tau_v, tau_v_cell);
				printf("cell volume, volume flux: %f, %f\n", volume_cell[iCy][iCx],
															 volume_flux);
				*/
				tau_v = getDoubleMin(tau_v_cell, tau_v);

				//printf("tau_e: %f, tau_v: %f\n", tau_e, tau_v);
			}
		}


		double tau = tau_e;
		//double tau = getDoubleMin(tau_e, tau_v);
		//double tau = getDoubleMin(tau_e, 1e-3);
		t += tau;
		++time_step;

	// Update cell specific momentum
	//			   specific total energy
	// and sync    density
	//			   pressure
		//printf("updating cell state...\n");
		for(size_t iCy = START_CELL_Y; iCy <= END_CELL_Y; ++iCy) {
			for(size_t iCx = START_CELL_X; iCx <= END_CELL_X; ++iCx) {
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
						attVec2(resul_force, corner_force);

						expan_power += vecMulVec2(corner_force, 
												velocity_node[iPy][iPx]);
					}
				}

				// update
				double contrib_sp_momentum[2];
				numMulVec2(-tau/mass_cell[iCy][iCx], 
							resul_force, contrib_sp_momentum);
				attVec2(sp_momentum_cell[iCy][iCx], contrib_sp_momentum);

				sp_energy_cell[iCy][iCx] += - tau * expan_power \
												/ mass_cell[iCy][iCx];

			}
		}

	// Update node position
		//printf("updating node position...\n");
		for(size_t iPy = START_NODE_Y; iPy <= END_NODE_Y; ++iPy) {
			for(size_t iPx = START_NODE_X; iPx <= END_NODE_X; ++iPx) {
				// update a 2d vector: x_p += u_p * tau
				double position_update[2];
				numMulVec2(tau, velocity_node[iPy][iPx], position_update);
				attVec2(coor_node[iPy][iPx], position_update);
			}
		}

	// Synchronize cell volume with the moved mesh (geometric conservation).
		for(size_t iCy = START_CELL_Y; iCy <= END_CELL_Y; ++iCy) {
			for(size_t iCx = START_CELL_X; iCx <= END_CELL_X; ++iCx) {
				volume_cell[iCy][iCx] = calCellVolume(
						coor_node[iCy][iCx], coor_node[iCy][iCx + 1],
						coor_node[iCy + 1][iCx + 1], coor_node[iCy + 1][iCx]);

				if (!isfinite(volume_cell[iCy][iCx]) || volume_cell[iCy][iCx] <= 0.0) {
					fprintf(stderr,
							"Invalid cell volume at step %zu, t=%.17g, cell=(%zu,%zu): %.17g\n",
							time_step, t, iCx, iCy, volume_cell[iCy][iCx]);
					return EXIT_FAILURE;
				}
			}
		}
		
	// Sync cell density
	//			 pressure
		//printf("updating cell state...\n");
		for(size_t iCy = START_CELL_Y; iCy <= END_CELL_Y; ++iCy) {
			for(size_t iCx = START_CELL_X; iCx <= END_CELL_X; ++iCx) {

				// sync density and pressure
				density_cell[iCy][iCx] = calCellDensity(volume_cell[iCy][iCx],
														mass_cell[iCy][iCx]);
				double sp_interior = calCellSpInterior(sp_energy_cell[iCy][iCx],
													  sp_momentum_cell[iCy][iCx]);
				pressure_cell[iCy][iCx] = calCellPressure(density_cell[iCy][iCx],
														  sp_interior, GAMMA);
			}
		}

	
	// Put out the data
		if (time_step % OUTPUT_FREQUNCY == 0) {
			printf("Put out data at time: %.16f\n", t);
			putOutNodeCoor(t, NUM_TOTAL_NODES_X, NUM_TOTAL_NODES_Y,
						   START_NODE_X, END_NODE_X, START_NODE_Y, END_NODE_Y, 
						   coor_node);
		}

		//起点不对
	}

	printf("Time advance is done\n");
	printf("time_end: %f\n", t);
	printf("Put out data at time: %f\n", t);
	putOutNodeCoor(t, NUM_TOTAL_NODES_X, NUM_TOTAL_NODES_Y,
				   START_NODE_X, END_NODE_X, START_NODE_Y, END_NODE_Y, 
				   coor_node);


}
