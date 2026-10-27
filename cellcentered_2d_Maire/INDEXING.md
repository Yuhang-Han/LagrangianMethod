# Indexing-only 2D refactor

This revision supersedes the earlier broad refactor. It starts from the original
maire2d source at upstream commit
3285b4f9f3bfa05eff836ad1a3cfebbf1c396681 and preserves the original main() structure.
The original backup is maire2d-backup-3285b4f.zip.

## Scope

The production changes are limited to seven files: maire2d.c and the source/header
pairs for mesh_calculation, initial_bound_condition, and data_output.

main() keeps its original sections, stack allocation, inline calculations,
coefficient values, and execution order. Its edits adapt array declarations,
index ranges, traversals, subscripts, and index-dependent function calls.
There is no simulation wrapper or extracted nodal/cell update function.

The physical and linear-algebra modules, compile/debug scripts, and existing
linearAlgebra_test.c match the original files byte-for-byte. calCornerLength(),
calCornerNormal(), and all scalar initial/boundary functions also retain their
original definitions and declarations.

## Names and inclusive ranges

C arrays use field[y][x], with vector component last:

~~~c
double coor_node[NUM_GHOST_NODES_Y][NUM_GHOST_NODES_X][2];
double pressure_cell[NUM_GHOST_CELLS_Y][NUM_GHOST_CELLS_X];

coor_node[node_y][node_x][0];
pressure_cell[cell_y][cell_x];
~~~

NUM_CELLS_X/Y and NUM_NODES_X/Y count physical entities. The original
NUM_GHOST_CELLS_X/Y and NUM_GHOST_NODES_X/Y names are retained: they count the
entire allocated extent, including physical entities and one ghost layer per side.

Bounds are explicit and inclusive:

~~~c
const size_t START_CELL_X = 1;
const size_t END_CELL_X = NUM_CELLS_X;
const size_t START_CELL_Y = 1;
const size_t END_CELL_Y = NUM_CELLS_Y;
const size_t START_NODE_X = 1;
const size_t END_NODE_X = NUM_NODES_X;
const size_t START_NODE_Y = 1;
const size_t END_NODE_Y = NUM_NODES_Y;
~~~

START_TOTAL_CELL_X/Y and START_TOTAL_NODE_X/Y are zero. Their END_TOTAL
counterparts are the corresponding allocated count minus one.

For the default 5 x 5 physical mesh:

| Entity | Physical indices per axis | Ghost indices per axis | Allocated shape |
| --- | --- | --- | --- |
| Cells | 1 through 5 | 0 and 6 | 7 x 7 |
| Nodes | 1 through 6 | 0 and 7 | 8 x 8 |

node_x/node_y and cell_x/cell_y are actual array indices. The dictionaries,
local/global conversions, and numbered-corner lookup table are removed.

The four cells around a node are traversed directly with cell_y from node_y-1
to node_y and cell_x from node_x-1 to node_x. The four nodes of a cell are
traversed with node_y from cell_y to cell_y+1 and node_x from cell_x to cell_x+1.
Both traversals keep the original accumulation order: southwest, southeast,
northwest, northeast. Ghost traversals cover the allocated rectangle and skip
the physical rectangle, including all four exterior corners.

## Functions adapted for indexing

| Function | Change |
| --- | --- |
| calNeighVertexIdx | Accepts actual node and cell x/y indices and writes the neighboring node x/y indices. Direction +1 selects the next CCW vertex and -1 selects the previous vertex. |
| initializeNodeCoor | Accepts node_x/node_y directly and retains the original coordinate arithmetic. |
| boundarilizeCoor | Supplies the previously declared but missing boundary-coordinate function using the same constant Cartesian coordinate formula. |
| putOutNodeCoor | Accepts the padded 2D array and inclusive physical-node bounds. Retains void return, CSV schema, error behavior, and row separators. |

The neighbor helper has the same valid-input preconditions as the original:
its node belongs to the cell, direction is +1 or -1, output pointers are valid,
and the caller ensures that array indices are in bounds. It introduces no
separate logical index or corner numbering.

The coordinate formula is deliberately unchanged:

~~~c
coor_node[0] = 0. + node_x * cell_h_0;
coor_node[1] = 0. + node_y * cell_h_1;
~~~

boundarilizeCoor() calls initializeNodeCoor() to provide constant coordinates.
It does not introduce a moving-mesh extrapolation boundary policy.
CSV i/j are the actual array subscripts; the node column is only a sequential
output record number, retained for file-format compatibility.

## Preserved solver behavior

The original time increment remains before the nodal solve. The fixed step,
inline volume/momentum/energy formulas, normal convention, and division-based
output condition are unchanged. No extra initial/final outputs or geometric
stop conditions are added.

The original pressure argument coor_node[2] is translated literally to
coor_node[0][2]. This retains its existing uninitialized-ghost input at that
point in main(); it does not substitute a cell center or repair pressure
initialization. The original coordinate offset is likewise preserved.
These legacy behaviors are outside this indexing-only revision.

## Build and focused verification

From maire2d:

~~~sh
mkdir -p build
./compile
./tests/run_tests.sh
~~~

The original compile script compiles and links the program. It does not execute
main(). The focused test runner executes only:

1. calNeighVertexIdx: all corners and both directions, plus rectangular-mesh
   boundary-node connections to ghost cells.
2. initializeNodeCoor and boundarilizeCoor: explicit physical and ghost
   coordinate samples on a rectangular mesh.
3. putOutNodeCoor: the physical output window, padded stride, actual indices,
   row separators, ghost exclusion, and invalid inputs.

No solver integration or Sedov run is part of this verification. The original
algebra test is preserved but not run, and the obsolete flattened-index test is
replaced by the neighbor-index test. The test runner uses GNU C11, strict
warnings, and undefined-behavior/bounds sanitizers; GNU C keeps the original
M_PI usage without changing unrelated scalar code.

Verification result: compilation/linking succeeded, and all three focused test
executables passed. A read-only comparison against the original confirmed the
main section/order and unchanged numerical expressions. The solver was not run.

## Apply the revised patch

Use maire2d-indexing-only.zip as a complete replacement source directory, or
apply maire2d-indexing-only.patch to a checkout of the stated original commit:

~~~sh
git apply --check /path/to/maire2d-indexing-only.patch
git apply /path/to/maire2d-indexing-only.patch
~~~

This patch is against the original backup, not the earlier broad refactor.
The source archive excludes generated build files. Remote GitHub files have
not been changed.
