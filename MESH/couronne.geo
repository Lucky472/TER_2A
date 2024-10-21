// Gmsh project created on Mon Oct 21 11:50:52 2024
SetFactory("OpenCASCADE");
//+
Circle(1) = {0, 0, 0, 2, 0, 2*Pi};
//+
Circle(2) = {0, 0, 0, 1, 0, 2*Pi};
//+
Curve Loop(1) = {1};
//+
Curve Loop(2) = {2};
//+
Plane Surface(1) = {1, 2};
//+
Physical Curve("Cercle_int", 10) = {2};
//+
Physical Curve("Cercle_ext", 11) = {1};
//+
Physical Surface("Couronne", 100) = {1};
