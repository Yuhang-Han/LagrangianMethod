#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>

typedef struct {
    // topology info
    size_t (*vertices)[4];
    size_t (*children)[4];
    // physical info
    double *mass;
    double *density;
    double *internal;
} Grid;

typedef struct {
    // topology info
    size_t *center;
    size_t (*children)[4];
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
    size_t node[2];
    size_t boundary;
} Edge;



int main()
{
	const size_t N = 240;

	double (*coordinates)[2] = malloc(4*N * sizeof(*coordinates));

	Grid grid;
	grid.vertices = malloc(N * sizeof(*grid.vertices));
	grid.children = malloc(N * sizeof(*grid.children));
	grid.mass     = malloc(N * sizeof(double));
	grid.density  = malloc(N * sizeof(double));
	grid.internal = malloc(N * sizeof(double));

	DualGrid dualgrid;
	dualgrid.center   = malloc(N * sizeof(size_t));
	dualgrid.children = malloc(N * sizeof(*dualgrid.children));
	dualgrid.mass     = malloc(N * sizeof(double));
	dualgrid.velocity = malloc(N * sizeof(*dualgrid.velocity));

	Corners corners;
	corners.vertices = malloc(N * sizeof(*corners.vertices));
	corners.mass     = malloc(4*N * sizeof(double));
	corners.pressure = malloc(4*N * sizeof(double));
	corners.force    = malloc(4*N * sizeof(double));

	Edge *edges = NULL;
	edges = malloc(N * sizeof(Edge));


/****************************
* Read the mesh
****************************/


    FILE *fp = fopen("distorted.msh", "r");

    if (!fp) {
        perror("distorted.msh");
        return 1;
    }

    char line[256];

    size_t numNodes = 0;
    size_t numQuads = 0;
    size_t numEdges = 0;

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

                coordinates[id - 1][0] = x;
                coordinates[id - 1][1] = y;
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
	                	
						edges[numEdges].node[0] = bn1 - 1;
						edges[numEdges].node[1] = bn2 - 1;

						edges[numEdges].boundary = physical_tag;

						numEdges++;
						
						break;
						
                	case 3:
	                	size_t qn1, qn2, qn3, qn4;

	                	fscanf(fp, "%zu %zu %zu %zu",
	                	       &qn1, &qn2, &qn3, &qn4);

	                	grid.vertices[numQuads][0] = qn1 - 1;
	                	grid.vertices[numQuads][1] = qn2 - 1;
	                	grid.vertices[numQuads][2] = qn3 - 1;
	                	grid.vertices[numQuads][3] = qn4 - 1;

	                	numQuads++;
	                	
	                	break;
	                	
	                default:
	                	fgets(line, sizeof(line), fp);
                }
            }
        }
    }

    fclose(fp);

    printf("Nodes: %zu\n", numNodes);
    printf("Quads: %zu\n", numQuads);
    printf("Edges: %zu\n", numEdges);

    printf("\nNode coordinates:\n");

    for (size_t i = 0; i < numNodes; i++) {
        printf("%zu : %g %g\n",
               i,
               coordinates[i][0],
               coordinates[i][1]);
    }

    printf("\nQuad connectivity:\n");

    for (size_t e = 0; e < numQuads; e++) {

        printf("%zu : %zu %zu %zu %zu\n",
               e,
               grid.vertices[numQuads][0],
               grid.vertices[numQuads][1],
               grid.vertices[numQuads][2],
               grid.vertices[numQuads][3]);
    }

    printf("\nEdge connectivity:\n");

    for (size_t e = 0; e < numEdges; e++) {

        printf("%zu : %zu %zu\n",
               e,
               edges[e].node[0],
               edges[e].node[1]);
    }


/****************************
* Construct topology
****************************/

	size_t iCorner = 0;
	size_t numExNodes = 0;
	for(size_t iCell = 0; iCell < numQuads; ++iCell) {

		double coor1x = coordinates[ grid.vertices[iCell][0] ][0];
		double coor1y = coordinates[ grid.vertices[iCell][0] ][1];
		double coor2x = coordinates[ grid.vertices[iCell][1] ][0];
		double coor2y = coordinates[ grid.vertices[iCell][1] ][1];
		double coor3x = coordinates[ grid.vertices[iCell][2] ][0];
		double coor3y = coordinates[ grid.vertices[iCell][2] ][1];
		double coor4x = coordinates[ grid.vertices[iCell][3] ][0];
		double coor4y = coordinates[ grid.vertices[iCell][3] ][1];
		
		double cellCenterx = 1/4 * (coor1x + coor2x + coor3x + coor4x);
		double cellCentery = 1/4 * (coor1y + coor2y + coor3y + coor4y);

		double cellBound1Midpointx = 1/2 * (coor1x + coor2x);
		double cellBound1Midpointy = 1/2 * (coor1y + coor2y);
		double cellBound2Midpointx = 1/2 * (coor2x + coor3x);
		double cellBound2Midpointy = 1/2 * (coor2y + coor3y);
		double cellBound3Midpointx = 1/2 * (coor3x + coor4x);
		double cellBound3Midpointy = 1/2 * (coor3y + coor4y);
		double cellBound4Midpointx = 1/2 * (coor4x + coor1x);
		double cellBound4Midpointy = 1/2 * (coor4y + coor1y);

		coordinates[numNodes+numExNodes][0] = cellCenterx;
		coordinates[numNodes+numExNodes][1] = cellCentery;

		coordinates[numNodes+numExNodes+1][0] = cellBound1Midpointx;
		coordinates[numNodes+numExNodes+1][1] = cellBound1Midpointy;
		coordinates[numNodes+numExNodes+2][0] = cellBound2Midpointx;
		coordinates[numNodes+numExNodes+2][1] = cellBound2Midpointy;
		coordinates[numNodes+numExNodes+3][0] = cellBound3Midpointx;
		coordinates[numNodes+numExNodes+3][1] = cellBound3Midpointy;
		coordinates[numNodes+numExNodes+4][0] = cellBound4Midpointx;
		coordinates[numNodes+numExNodes+4][1] = cellBound4Midpointy;

		corners.vertices[iCorner][0] = numNodes+numExNodes;
		corners.vertices[iCorner][1] = numNodes+numExNodes+1;
		corners.vertices[iCorner][2] = gird.vertices[1];
		corners.vertices[iCorner][3] = numNodes+numExNodes+2;

		corners.vertices[iCorner+1][0] = numNodes+numExNodes;
		corners.vertices[iCorner+1][1] = numNodes+numExNodes+2;
		corners.vertices[iCorner+1][2] = gird.vertices[2];
		corners.vertices[iCorner+1][3] = numNodes+numExNodes+3;
		
		corners.vertices[iCorner][0] = numNodes+numExNodes;
		corners.vertices[iCorner][1] = numNodes+numExNodes+3;
		corners.vertices[iCorner][2] = gird.vertices[3];
		corners.vertices[iCorner][3] = numNodes+numExNodes+4;
		
		corners.vertices[iCorner][0] = numNodes+numExNodes;
		corners.vertices[iCorner][1] = numNodes+numExNodes+4;
		corners.vertices[iCorner][2] = gird.vertices[0];
		corners.vertices[iCorner][3] = numNodes+numExNodes+1;

		numExNodes += 5;
		iCorner += 4;


		
	}

















    


	free(grid.vertices);
	free(grid.children);
	free(grid.mass);
	free(grid.density);
	free(grid.internal);

	free(dualgrid.center);
	free(dualgrid.children);
	free(dualgrid.mass);
	free(dualgrid.velocity);

	free(corners.vertices);
	free(corners.mass);
	free(corners.pressure);
	free(corners.force);
	
	free(coordinates);
    free(edges);
	
}
