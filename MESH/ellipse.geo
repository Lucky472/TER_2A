// Gmsh project created on Mon Oct 21 16:14:38 2024
SetFactory("OpenCASCADE");
//+
Ellipse(1) = {0, 0, 0, 2, 1, 0, 2*Pi};
//+
Circle(2) = {-0.75, 0, 0, 0.75, 0, 2*Pi};
//+
Circle(3) = {1, 0, 0, 0.5, 0, 2*Pi};
//+
Curve Loop(1) = {1};
//+
Curve Loop(2) = {2};
//+
Curve Loop(3) = {3};
//+
Plane Surface(1) = {1, 2, 3};
//+
Physical Curve("Cercle_int_gauche", 10) = {2};
//+
Physical Curve("Cercle_int_droit", 11) = {3};
//+
Physical Curve("Ellipse_trouee", 20) = {1};
//+
Physical Surface("Ellipse", 100) = {1};
