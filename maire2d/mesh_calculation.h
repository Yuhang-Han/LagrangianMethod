#include "linear_algebra_2d.h"

size_t localToGlobalCellIdx(size_t idxP, size_t idxCP);
size_t localToGlobalNodeIdx(size_t idxC, size_t idxPC);

double calCornerLength(size_t idxC, size_t idxP, double* meshinfo);
double calCornerNormal(size_t idxC, size_t idxP, double* meshinfo);
double calCornerMatrix(size_t idxC, size_t idxP, double* meshinfo);
