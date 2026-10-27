#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    double x;
    double y;
} Node;

typedef struct {
    int node[4];
} Quad;

typedef struct {
    int node[2];
    int boundary;
} Edge;



// structure of .msh

// Node:
// 4 -0.3 1.8 0
// node id, x, y, z

// Element: 
// 57 3 2 0 1 83 87 89 88
// element ID, element type, number of tags, physical tag, geometrical tag
// node IDs


int main(void)
{
    FILE *fp = fopen("distorted.msh", "r");

    if (!fp) {
        perror("distorted.msh");
        return 1;
    }

    char line[256];

    Node *nodes = NULL;
    Quad *quads = NULL;
    Edge *edges = NULL;

    int numNodes = 0;
    int numQuads = 0;
    int numEdges = 0;

    while (fgets(line, sizeof(line), fp)) {

        /* --------------------
           Read nodes
           -------------------- */
        if (strncmp(line, "$Nodes", 6) == 0) {

            fscanf(fp, "%d", &numNodes);

            nodes = malloc(numNodes * sizeof(Node));

            for (int i = 0; i < numNodes; i++) {

                int id;
                double x, y, z;

                fscanf(fp, "%d %lf %lf %lf",
                       &id, &x, &y, &z);

                nodes[id - 1].x = x;
                nodes[id - 1].y = y;
            }
        }

        /* --------------------
           Read elements
           -------------------- */
        if (strncmp(line, "$Elements", 9) == 0) {

            int numElements;

            fscanf(fp, "%d", &numElements);

            /* Maximum possible allocation.
               We shrink conceptually by numQuads. */
            quads = malloc(numElements * sizeof(Quad));
            edges = malloc(numElements * sizeof(Edge));

            for (int i = 0; i < numElements; i++) {

                int id;
                int type;
                int numTags;

                fscanf(fp, "%d %d %d",
                       &id, &type, &numTags);



				// get tags
				int physical_tag = 0;
				//int geometry_tag = 0;
				
                for (int j = 0; j < numTags; j++) {
                    int tag;
                    fscanf(fp, "%d", &tag);

				    if (j == 0)
				        physical_tag = tag;

				    if (j == 1){
				        //geometry_tag = tag;
				    }
                }

                switch (type) {
                	case 15:	// point element, ignore
	                	int n;
	                	fscanf(fp, "%d", &n);
	                	break;
	                	
					case 1:		// boundary element
	                	int bn1, bn2;
	                	fscanf(fp, "%d %d", &bn1, &bn2);
	                	
						edges[numEdges].node[0] = bn1 - 1;
						edges[numEdges].node[1] = bn2 - 1;

						edges[numEdges].boundary = physical_tag;

						numEdges++;
						
						break;
						
                	case 3:
	                	int qn1, qn2, qn3, qn4;

	                	fscanf(fp, "%d %d %d %d",
	                	       &qn1, &qn2, &qn3, &qn4);

	                	quads[numQuads].node[0] = qn1 - 1;
	                	quads[numQuads].node[1] = qn2 - 1;
	                	quads[numQuads].node[2] = qn3 - 1;
	                	quads[numQuads].node[3] = qn4 - 1;

	                	numQuads++;
	                	
	                	break;
	                	
	                default:
	                	fgets(line, sizeof(line), fp);
                }
            }
        }
    }

    fclose(fp);

    printf("Nodes: %d\n", numNodes);
    printf("Quads: %d\n", numQuads);
    printf("Edges: %d\n", numEdges);

    printf("\nNode coordinates:\n");

    for (int i = 0; i < numNodes; i++) {
        printf("%d : %g %g\n",
               i,
               nodes[i].x,
               nodes[i].y);
    }

    printf("\nQuad connectivity:\n");

    for (int e = 0; e < numQuads; e++) {

        printf("%d : %d %d %d %d\n",
               e,
               quads[e].node[0],
               quads[e].node[1],
               quads[e].node[2],
               quads[e].node[3]);
    }

    printf("\nEdge connectivity:\n");

    for (int e = 0; e < numEdges; e++) {

        printf("%d : %d %d\n",
               e,
               edges[e].node[0],
               edges[e].node[1]);
    }


    free(nodes);
    free(quads);
    free(edges);

    return 0;
}
