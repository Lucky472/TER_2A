// Gmsh project created on Sat Oct 15 15:28:09 2022
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Physical Line("ouest", 10) = {4};
//+
Physical Line("est", 11) = {2};
//+
Physical Line("nord-sud", 20) = {1, 3};
//+
Physical Surface("domaine", 30) = {1};
