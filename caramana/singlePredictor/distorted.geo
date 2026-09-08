Point(1) = {0,   0,   0, 0.25};
Point(2) = {4,   0.2, 0, 0.25};
Point(3) = {3.5, 2.2, 0, 0.25};
Point(4) = {-0.3,1.8, 0, 0.25};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Physical Curve("bottom", 1) = {1};
Physical Curve("right",  2) = {2};
Physical Curve("top",    3) = {3};
Physical Curve("left",   4) = {4};

Curve Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

Physical Surface("domain", 10) = {1};

Mesh.RecombinationAlgorithm = 3;
Recombine Surface {1};
