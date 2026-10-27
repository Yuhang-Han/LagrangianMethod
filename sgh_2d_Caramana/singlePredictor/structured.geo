// Rectangle dimensions
Lx = 2.0;
Ly = 2.0;

// Number of mesh points along each direction
Nx = 6;
Ny = 6;

// Geometry points
Point(1) = {0,  0,  0};
Point(2) = {Lx, 0,  0};
Point(3) = {Lx, Ly, 0};
Point(4) = {0,  Ly, 0};

// Boundary lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Surface
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Structured mesh
Transfinite Curve {1, 3} = Nx;
Transfinite Curve {2, 4} = Ny;

Transfinite Surface {1};

// Convert triangles into quadrilaterals
Recombine Surface {1};

// Optional physical groups
Physical Curve("bottom") = {1};
Physical Curve("right")  = {2};
Physical Curve("top")    = {3};
Physical Curve("left")   = {4};

Physical Surface("domain") = {1};
