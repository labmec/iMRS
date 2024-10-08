// Gmsh project created on Fri Jul  5 14:35:05 2024
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {345, 0, 0, 1.0};
//+
Point(3) = {345, 1000, 0, 1.0};
//+
Point(4) = {0, 1000, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {2, 3, 4, 1};
//+
Plane Surface(1) = {1};
//+
Physical Curve("bcl", 5) = {4};
//+
Physical Curve("bcr", 6) = {2};
//+
Physical Curve("bct", 7) = {3};
//+
Physical Curve("bcb", 8) = {1};
//+
Physical Surface("dom", 9) = {1};
//+
Transfinite Surface {1};
//+
Transfinite Curve {4, 2} = 5 Using Progression 1;
//+
Transfinite Curve {3, 1, 3, 1} = 3 Using Progression 1;
//+
Recombine Surface {1};
