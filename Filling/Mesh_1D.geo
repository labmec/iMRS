// Gmsh project created on Fri Jul  5 14:35:05 2024
SetFactory("OpenCASCADE");
rin = 30.15; //sandscreen tube radius
rout = 30.15+345; //external radius
h = 1000; //height

ref =1;
xnels = 2^ref;
ynels = 2^ref; 

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


Curve Loop(1) = {1,2,3,4};

Plane Surface(1) = {1};

Physical Curve("inlet", 2) = {4};
//+
Physical Curve("outlet", 3) = {2};
//+
Physical Curve("noflux", 4) = {1,3};

Physical Surface("dom", 1) = {1};


Transfinite Surface {:};
//+
Transfinite Curve {1,3} = xnels + 1;
//+
Transfinite Curve {2,4} = ynels +1 ;
//+
Recombine Surface {:};
