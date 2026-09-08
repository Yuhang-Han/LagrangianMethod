#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>

typedef struct {
    // topology info
    size_t (*vertices)[4];
    size_t (*edges)[4];
    size_t (*corners)[4];
	// geometry info
    size_t *center;
    // physical info
    double *mass;
    double *density;
    double *internal;
} Grid;

typedef struct {
    // topology info
    size_t *center;
    size_t (*corners)[4];
    // physical info
    double *mass;
    double (*velocity)[2];
} DualGrid;

typedef struct {
	size_t *vertices;
	
    double *mass;
    double *pressure;
    double *force;
} Corners;

typedef struct {
    size_t (*vertices)[2];
    size_t *mid;
    size_t *boundary;
} Edges;

typedef struct {
    size_t *id;
	size_t (*neighbor_corners)[4];
} Nodes;

int main()
{
	const size_t N = 240;

	double (*coor)[2] = malloc(4*N * sizeof(*coor));
	double (*coor_mid)[2] = malloc(4*N * sizeof(*coor_mid));
	double (*coor_center)[2] = malloc(4*N * sizeof(*coor_center));

	Grid grid;
	grid.vertices = malloc(N * sizeof(*grid.vertices));
	grid.edges    = malloc(N * sizeof(*grid.edges));
	grid.corners  = malloc(N * sizeof(*grid.corners));
	grid.center   = malloc(N * sizeof(double));
	grid.mass     = malloc(N * sizeof(double));
	grid.density  = malloc(N * sizeof(double));
	grid.internal = malloc(N * sizeof(double));

	DualGrid dualgrid;
	dualgrid.center   = malloc(N * sizeof(size_t));
	dualgrid.corners  = malloc(N * sizeof(*dualgrid.corners));
	dualgrid.mass     = malloc(N * sizeof(double));
	dualgrid.velocity = malloc(N * sizeof(*dualgrid.velocity));

	Corners corners;
	corners.vertices = malloc(N * sizeof(*corners.vertices));
	corners.mass     = malloc(4*N * sizeof(double));
	corners.pressure = malloc(4*N * sizeof(double));
	corners.force    = malloc(4*N * sizeof(double));

	Edges edges0;
	edges0.vertices = malloc(4*N * sizeof(*edges0.vertices));
	edges0.mid	    = malloc(N * sizeof(double));
	edges0.boundary = malloc(N * sizeof(size_t));



	Nodes nodes;
	nodes.id 			   = malloc(N * sizeof(size_t));
	nodes.neighbor_corners = malloc(N * sizeof(nodes.neighbor_corners));


/****************************
* Read the mesh
****************************/


    FILE *fp = fopen("structured.msh", "r");

    if (!fp) {
        perror("structured.msh");
        return 1;
    }

    char line[256];

    size_t numNodes = 0;
    size_t numCells = 0;
    size_t numedges0 = 0;

    while (fgets(line, sizeof(line), fp)) {

        /* --------------------
           Read nodes
           -------------------- */
        if (strncmp(line, "$Nodes", 6) == 0) {

            fscanf(fp, "%zu", &numNodes);

            for (size_t i = 0; i < numNodes; i++) {

                size_t id;
                double x, y, z;

                fscanf(fp, "%zu %lf %lf %lf",
                       &id, &x, &y, &z);

                coor[id - 1][0] = x;
                coor[id - 1][1] = y;
            }
        }

        /* --------------------
           Read elements
           -------------------- */
        if (strncmp(line, "$Elements", 9) == 0) {

            size_t numElements;

            fscanf(fp, "%zu", &numElements);

            for (size_t i = 0; i < numElements; i++) {

                size_t id;
                size_t type;
                size_t numTags;

                fscanf(fp, "%zu %zu %zu",
                       &id, &type, &numTags);

				// get tags
				size_t physical_tag = 0;
				//size_t geometry_tag = 0;
				
                for (size_t j = 0; j < numTags; j++) {
                    size_t tag;
                    fscanf(fp, "%zu", &tag);

				    if (j == 0)
				        physical_tag = tag;

				    if (j == 1){
				        //geometry_tag = tag;
				    }
                }

                switch (type) {
                	case 15:	// posize_t element, ignore
	                	size_t n;
	                	fscanf(fp, "%zu", &n);
	                	break;
	                	
					case 1:		// boundary element
	                	size_t bn1, bn2;
	                	fscanf(fp, "%zu %zu", &bn1, &bn2);
	                	
						edges0.vertices[numedges0][0] = \
								(bn1 < bn2 ? (bn1-1) : (bn2-1));	// min
						edges0.vertices[numedges0][1] = \
								(bn1 > bn2 ? (bn1-1) : (bn2-1));	// max

						edges0.boundary[numedges0] = physical_tag;

						numedges0++;
						
						break;
						
                	case 3:
	                	size_t qn1, qn2, qn3, qn4;

	                	fscanf(fp, "%zu %zu %zu %zu",
	                	       &qn1, &qn2, &qn3, &qn4);

	                	grid.vertices[numCells][0] = qn1 - 1;
	                	grid.vertices[numCells][1] = qn2 - 1;
	                	grid.vertices[numCells][2] = qn3 - 1;
	                	grid.vertices[numCells][3] = qn4 - 1;

	                	numCells++;
	                	
	                	break;
	                	
	                default:
	                	fgets(line, sizeof(line), fp);
                }
            }
        }
    }

    fclose(fp);


    printf("Nodes: %zu\n", numNodes);
    printf("Quads: %zu\n", numCells);
    printf("edges0: %zu\n", numedges0);

    printf("\nNode coordinates:\n");

    for (size_t i = 0; i < numNodes; i++) {
        printf("%zu : %g %g\n",
               i,
               coor[i][0],
               coor[i][1]);
    }

    printf("\nQuad connectivity:\n");

    for (size_t e = 0; e < numCells; e++) {

        printf("%zu : %zu %zu %zu %zu\n",
               e,
               grid.vertices[e][0],
               grid.vertices[e][1],
               grid.vertices[e][2],
               grid.vertices[e][3]);
    }

    printf("\nEdge0 connectivity:\n");

    for (size_t e = 0; e < numedges0; e++) {

        printf("%zu : %zu %zu\n",
               e,
               edges0.vertices[e][0],
               edges0.vertices[e][1]);
    }
    

/****************************
* Construct interior edges
****************************/

	// Create egdes
	size_t cEdge0 = numedges0;

	for (size_t iCell = 0; iCell < numCells; ++ iCell) {
		for (size_t iCellEdge = 0; iCellEdge < 4; ++iCellEdge) {
			size_t cell_edge_id0 = iCellEdge;
			size_t cell_edge_id1 = (iCellEdge + 1) % 4;

			size_t id0 = grid.vertices[iCell][cell_edge_id0];
			size_t id1 = grid.vertices[iCell][cell_edge_id1];

			edges0.vertices[cEdge0][0] = (id0 < id1 ? id0 : id1);		// min
			edges0.vertices[cEdge0][1] = (id0 > id1 ? id0 : id1);		// max

			edges0.boundary[cEdge0] = 0;		// interior edge

			++cEdge0;
		}
	}

	numedges0 = cEdge0;


	// Merge identical edges
	Edges edges;
	edges.vertices = malloc(4*N * sizeof(*edges.vertices));
	edges.mid	   = malloc(N * sizeof(double));
	edges.boundary = malloc(N * sizeof(size_t));
	// initialize
	for(size_t i = 0; i < 4*N; ++i) {
		edges.vertices[i][0] = 0;
		edges.vertices[i][1] = 0;
	}

	size_t cEdge = 0;
	for(size_t iEdge0 = 0; iEdge0 < numedges0; ++iEdge0) {
	
		size_t repFlag = 0;
		
		for(size_t iEdge = 0; iEdge < cEdge; ++iEdge) {
			if( (edges0.vertices[iEdge0][0] == edges.vertices[iEdge][0]) 
				&& (edges0.vertices[iEdge0][1] == edges.vertices[iEdge][1]) )
			{
				repFlag = 1;
			}

		}

		if(repFlag == 0) {
			edges.vertices[cEdge][0] = edges0.vertices[iEdge0][0];
			edges.vertices[cEdge][1] = edges0.vertices[iEdge0][1];

			edges.boundary[cEdge] = edges0.boundary[iEdge0];

			++cEdge;
		}
	}

	size_t numEdges = cEdge;


	free(edges0.vertices);
	free(edges0.mid);
	free(edges0.boundary);
	

    printf("\nEdge connectivity:\n");
    printf("Edges: %zu \n", numEdges);

    for (size_t e = 0; e < numEdges; e++) {

        printf("%zu : %zu %zu\n",
               e,
               edges.vertices[e][0],
               edges.vertices[e][1]);
    }
    

/****************************
* Construct geometry
****************************/

	// midpoint
	size_t cMid = 0;
	for (size_t iEdge = 0; iEdge < numEdges; ++ iEdge) {
		double leftx  = coor[ edges.vertices[iEdge][0] ][0];
		double lefty  = coor[ edges.vertices[iEdge][0] ][1];
		double rightx = coor[ edges.vertices[iEdge][1] ][0];
		double righty = coor[ edges.vertices[iEdge][1] ][1];

		double midx = 1/2.0 * (leftx + rightx);
		double midy = 1/2.0 * (lefty + righty);

		coor_mid[cMid][0] = midx;
		coor_mid[cMid][1] = midy;
		
		edges.mid[iEdge] = cMid;
		
		++cMid;
	}

/*
    printf("\nEdge midpoint:\n");
    for (size_t e = 0; e < numEdges; e++) {
        printf("%zu : left(%f,%f) rifht(%f,%f) mid(%f,%f)\n",
               e,
               coor[ edges.vertices[e][0] ][0],
               coor[ edges.vertices[e][0] ][1],
               coor[ edges.vertices[e][1] ][0],
               coor[ edges.vertices[e][1] ][1],
               coor_mid[ edges.mid[e] ][0],
               coor_mid[ edges.mid[e] ][1] );
    }
*/

	// cell center
	size_t cCenter = 0;
	for (size_t iCell = 0; iCell < numCells; ++ iCell) {
		double coor1x = coor[ grid.vertices[iCell][0] ][0];
		double coor1y = coor[ grid.vertices[iCell][0] ][1];
		double coor2x = coor[ grid.vertices[iCell][1] ][0];
		double coor2y = coor[ grid.vertices[iCell][1] ][1];
		double coor3x = coor[ grid.vertices[iCell][2] ][0];
		double coor3y = coor[ grid.vertices[iCell][2] ][1];
		double coor4x = coor[ grid.vertices[iCell][3] ][0];
		double coor4y = coor[ grid.vertices[iCell][3] ][1];
		
		double centerx = 1/4.0 * (coor1x + coor2x + coor3x + coor4x);
		double centery = 1/4.0 * (coor1y + coor2y + coor3y + coor4y);

		coor_center[cCenter][0] = centerx;
		coor_center[cCenter][1] = centery;

		grid.center[iCell] = cCenter;

		++cCenter;
	}

/*
	printf("\n Cell center: \n");
    for (size_t iCell = 0; iCell < numCells	; iCell++) {
        printf("%zu : v1(%.2f,%.2f) v2(%.2f,%.2f) v3(%.2f,%.2f) v4(%.2f,%.2f) center(%.2f,%.2f)\n",
               	iCell,
				coor[ grid.vertices[iCell][0] ][0],
				coor[ grid.vertices[iCell][0] ][1],
				coor[ grid.vertices[iCell][1] ][0],
				coor[ grid.vertices[iCell][1] ][1],
				coor[ grid.vertices[iCell][2] ][0],
				coor[ grid.vertices[iCell][2] ][1],
				coor[ grid.vertices[iCell][3] ][0],
				coor[ grid.vertices[iCell][3] ][1],
               	coor_center[ grid.center[iCell] ][0],
               	coor_center[ grid.center[iCell] ][1] );
    }
*/

/****************************
* Construct topology
****************************/

	// cell--edge
	for (size_t iCell = 0; iCell < numCells; ++ iCell) {

		size_t iCellEdge = 0;

		for (size_t iEdge = 0; iEdge < numEdges; ++ iEdge) {
		
			size_t public_nodes = 0;
			
			for (size_t iEdgeNode = 0; iEdgeNode < 2; ++iEdgeNode){
				for (size_t iCellNode = 0; iCellNode < 4; ++iCellNode) {
				
					if(grid.vertices[iCell][iCellNode] 
						== edges.vertices[iEdge][iEdgeNode]) {
							++public_nodes;
					}
				}
			}

			if(public_nodes == 2) {
				grid.edges[iCell][iCellEdge] = iEdge;
				++ iCellEdge;
			} 
			else if(public_nodes > 2) {
				printf("Wrong in constructing cell--edge topology\n");
				return 1;
			}
		}
		
	}


















/****************************
* Free the memory
****************************/


	free(coor);	free(coor_mid);	free(coor_center);

	free(grid.vertices);
	free(grid.edges);
	free(grid.corners);
	free(grid.center);
	free(grid.mass);
	free(grid.density);
	free(grid.internal);

	free(dualgrid.center);
	free(dualgrid.corners);
	free(dualgrid.mass);
	free(dualgrid.velocity);

	free(corners.vertices);
	free(corners.mass);
	free(corners.pressure);
	free(corners.force);

	free(edges.vertices);
	free(edges.mid);
	free(edges.boundary);

	free(nodes.id);
	free(nodes.neighbor_corners);
	

}
