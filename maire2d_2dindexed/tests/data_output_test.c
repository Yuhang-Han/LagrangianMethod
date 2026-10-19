#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "../data_output.h"

static void assertFileAbsent(const char *filename)
{
    FILE *fp = fopen(filename, "r");
    if (fp != NULL) {
        fclose(fp);
        assert(0 && "invalid output input unexpectedly created a file");
    }
}

static void testPhysicalWindow(void)
{
    const size_t NUM_CELLS_X = 3;
    const size_t NUM_CELLS_Y = 2;
    const size_t NUM_GHOST_LAYERS = 2;
    const size_t NUM_TOTAL_NODES_X = NUM_CELLS_X + 1 + 2 * NUM_GHOST_LAYERS;
    const size_t NUM_TOTAL_NODES_Y = NUM_CELLS_Y + 1 + 2 * NUM_GHOST_LAYERS;
    const size_t START_NODE_X = NUM_GHOST_LAYERS;
    const size_t END_NODE_X = START_NODE_X + NUM_CELLS_X;
    const size_t START_NODE_Y = NUM_GHOST_LAYERS;
    const size_t END_NODE_Y = START_NODE_Y + NUM_CELLS_Y;
    double coor_node[NUM_TOTAL_NODES_Y][NUM_TOTAL_NODES_X][2];

    for (size_t node_y = 0; node_y < NUM_TOTAL_NODES_Y; ++node_y) {
        for (size_t node_x = 0; node_x < NUM_TOTAL_NODES_X; ++node_x) {
            coor_node[node_y][node_x][0] = NAN;
            coor_node[node_y][node_x][1] = NAN;
        }
    }
    for (size_t node_y = START_NODE_Y; node_y <= END_NODE_Y; ++node_y) {
        for (size_t node_x = START_NODE_X; node_x <= END_NODE_X; ++node_x) {
            coor_node[node_y][node_x][0] = 10.0 * node_y + 0.125 * node_x;
            coor_node[node_y][node_x][1] = -20.0 * node_x + 0.25 * node_y;
        }
    }

    /* NaN ghosts must neither invalidate physical output nor appear in it. */
    const double time = 314159.25;
    const char *filename = "build/mesh_t_314159.25.csv";
    assert(putOutNodeCoor(time, NUM_TOTAL_NODES_X, NUM_TOTAL_NODES_Y,
                         START_NODE_X, END_NODE_X,
                         START_NODE_Y, END_NODE_Y, coor_node));
    FILE *fp = fopen(filename, "r");
    assert(fp != NULL);
    char line[256];
    assert(fgets(line, sizeof line, fp) != NULL);
    assert(strcmp(line, "t,node,i,j,x,y\n") == 0);
    size_t expected_record = 0;
    for (size_t node_y = START_NODE_Y; node_y <= END_NODE_Y; ++node_y) {
        for (size_t node_x = START_NODE_X; node_x <= END_NODE_X; ++node_x) {
            double output_time, x, y;
            size_t node_record, i, j;
            char newline;
            assert(fgets(line, sizeof line, fp) != NULL);
            assert(sscanf(line, "%lf,%zu,%zu,%zu,%lf,%lf%c",
                          &output_time, &node_record, &i, &j,
                          &x, &y, &newline) == 7);
            assert(newline == '\n');
            assert(output_time == time);
            assert(node_record == expected_record);
            assert(i == node_x && j == node_y);
            assert(x == coor_node[node_y][node_x][0]);
            assert(y == coor_node[node_y][node_x][1]);
            ++expected_record;
        }
        if (node_y < END_NODE_Y) {
            assert(fgets(line, sizeof line, fp) != NULL);
            assert(strcmp(line, "\n") == 0);
        }
    }
    assert(expected_record == (NUM_CELLS_X + 1) * (NUM_CELLS_Y + 1));
    assert(fgets(line, sizeof line, fp) == NULL);
    assert(feof(fp));
    assert(fclose(fp) == 0);
    assert(remove(filename) == 0);

    const char *invalid_filename = "build/mesh_t_314160.5.csv";
    assertFileAbsent(invalid_filename);
    assert(!putOutNodeCoor(314160.5, NUM_TOTAL_NODES_X, NUM_TOTAL_NODES_Y,
                          END_NODE_X, START_NODE_X,
                          START_NODE_Y, END_NODE_Y, coor_node));
    assert(!putOutNodeCoor(314160.5, NUM_TOTAL_NODES_X, NUM_TOTAL_NODES_Y,
                          START_NODE_X, NUM_TOTAL_NODES_X,
                          START_NODE_Y, END_NODE_Y, coor_node));
    assert(!putOutNodeCoor(314160.5, NUM_TOTAL_NODES_X, NUM_TOTAL_NODES_Y,
                          START_NODE_X, END_NODE_X,
                          START_NODE_Y, NUM_TOTAL_NODES_Y, coor_node));
    assert(!putOutNodeCoor(NAN, NUM_TOTAL_NODES_X, NUM_TOTAL_NODES_Y,
                          START_NODE_X, END_NODE_X,
                          START_NODE_Y, END_NODE_Y, coor_node));
    assertFileAbsent(invalid_filename);

    /* The last physical node is included in preflight validation. */
    coor_node[END_NODE_Y][END_NODE_X][1] = INFINITY;
    assert(!putOutNodeCoor(314160.5, NUM_TOTAL_NODES_X, NUM_TOTAL_NODES_Y,
                          START_NODE_X, END_NODE_X,
                          START_NODE_Y, END_NODE_Y, coor_node));
    assertFileAbsent(invalid_filename);
}

int main(void)
{
    testPhysicalWindow();
    puts("data_output_test: PASS");
    return 0;
}
