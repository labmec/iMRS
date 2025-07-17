// Gmsh project created on Fri Jul  5 14:35:05 2024
SetFactory("OpenCASCADE");
rin = 30.15; //sandscreen tube radius
rout = 30.15+345; //external radius
h = 1000; //height
gap = 100; //gap to purge
//+
Point(1) = {rin, 0, 0, 1.0};
//+
Point(2) = {rout, 0, 0, 1.0};
//+
Point(3) = {rout, h-gap, 0, 1.0};
//+
Point(4) = {rout, h, 0, 1.0};
//+
Point(5) = {rin, h, 0, 1.0};
//+
Point(6) = {rin, h-gap, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 1};
//+
Curve Loop(1) = {1,2,3,4,5,6};
//+
Plane Surface(1) = {1};
//+
Physical Curve("bottom", 2) = {1};
//+
Physical Curve("diffusor", 3) = {2,3};
//+
Physical Curve("lid", 4) = {4};
//+
Physical Curve("gap", 5) = {5};
//+
Physical Curve("sandscreen", 6) = {6};
//+
Physical Surface("dom", 1) = {1};
//+
Transfinite Surface {1} = {1,2,4,5};
//+
Transfinite Curve {1,4} = 21 Using Progression 1;
//+
Transfinite Curve {2,6} = 19 Using Progression 1;
//+
Transfinite Curve {3,5} = 3 Using Progression 1;
//+
Recombine Surface {1};
