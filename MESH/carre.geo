// Gmsh project created on Mon Oct 21 11:07:57 2024
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {0, 1, 0, 1.0};
//+
Point(4) = {1, 1, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 4};
//+
Line(3) = {4, 3};
//+
Line(4) = {3, 1};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Physical Curve("Ouest", 10) = {4};
//+
Physical Curve("Est", 11) = {2};
//+
Physical Curve("NordSud", 20) = {3, 1};
//+
Physical Surface("Carré", 100) = {1};
