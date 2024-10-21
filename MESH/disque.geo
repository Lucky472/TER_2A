// Gmsh project created on Mon Oct 21 11:31:16 2024
SetFactory("OpenCASCADE");
//+
Circle(1) = {0, 0, 0, 2, 0, 2*Pi};
//+
Curve Loop(1) = {1};
//+
Plane Surface(1) = {1};
//+
Physical Curve("Cercle_ext", 11) = {1};
//+
Physical Surface("Cercle_ext", 100) = {1};
