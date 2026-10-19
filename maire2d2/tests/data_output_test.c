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
        assert(0 && "invalid input unexpectedly created an output file");
    }
}

int main(void)
{
    const size_t NUM_GHOST_NODES_X = 8;
    const size_t NUM_GHOST_NODES_Y = 7;
    const size_t START_NODE_X = 2, END_NODE_X = 5;
    const size_t START_NODE_Y = 2, END_NODE_Y = 4;
    double coor_node[NUM_GHOST_NODES_Y][NUM_GHOST_NODES_X][2];

    for (size_t j = 0; j < NUM_GHOST_NODES_Y; ++j) {
        for (size_t i = 0; i < NUM_GHOST_NODES_X; ++i) {
            coor_node[j][i][0] = NAN;
            coor_node[j][i][1] = NAN;
        }
    }
    for (size_t j = START_NODE_Y; j <= END_NODE_Y; ++j) {
        for (size_t i = START_NODE_X; i <= END_NODE_X; ++i) {
            coor_node[j][i][0] = 10.0 * j + 0.125 * i;
            coor_node[j][i][1] = -20.0 * i + 0.25 * j;
        }
    }

    const double time = 314159.25;
    const char *filename = "build/mesh_t_314159.25.csv";
    putOutNodeCoor(time, NUM_GHOST_NODES_X, NUM_GHOST_NODES_Y,
                   START_NODE_X, END_NODE_X, START_NODE_Y, END_NODE_Y,
                   coor_node);
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
            assert(node_record == expected_record++);
            assert(i == node_x && j == node_y);
            assert(x == coor_node[node_y][node_x][0]);
            assert(y == coor_node[node_y][node_x][1]);
        }
        if (node_y < END_NODE_Y) {
            assert(fgets(line, sizeof line, fp) != NULL);
            assert(strcmp(line, "\n") == 0);
        }
    }
    assert(expected_record == 12);
    assert(fgets(line, sizeof line, fp) == NULL);
    assert(feof(fp));
    assert(fclose(fp) == 0);
    assert(remove(filename) == 0);

    const char *invalid_filename = "build/mesh_t_314160.5.csv";
    assertFileAbsent(invalid_filename);
    putOutNodeCoor(314160.5, NUM_GHOST_NODES_X, NUM_GHOST_NODES_Y,
                   END_NODE_X, START_NODE_X, START_NODE_Y, END_NODE_Y,
                   coor_node);
    assertFileAbsent(invalid_filename);
    putOutNodeCoor(314160.5, NUM_GHOST_NODES_X, NUM_GHOST_NODES_Y,
                   START_NODE_X, NUM_GHOST_NODES_X, START_NODE_Y, END_NODE_Y,
                   coor_node);
    assertFileAbsent(invalid_filename);
    putOutNodeCoor(314160.5, NUM_GHOST_NODES_X, NUM_GHOST_NODES_Y,
                   START_NODE_X, END_NODE_X, START_NODE_Y, NUM_GHOST_NODES_Y,
                   coor_node);
    assertFileAbsent(invalid_filename);

    coor_node[END_NODE_Y][END_NODE_X][1] = INFINITY;
    putOutNodeCoor(314160.5, NUM_GHOST_NODES_X, NUM_GHOST_NODES_Y,
                   START_NODE_X, END_NODE_X, START_NODE_Y, END_NODE_Y,
                   coor_node);
    assertFileAbsent(invalid_filename);

    puts("data_output_test: PASS");
    return 0;
}
