//  Geo file generated by Discrete Fracture Network methods 
// Fracture #1 

// POINTS DEFINITION 

h = 100000;

Point(1) = {0,0,0,h};
Point(2) = {2,0,0,h};
Point(3) = {0,2,0,h};
Point(4) = {2,2,0,h};
Point(5) = {0,0,2,h};
Point(6) = {2,0,2,h};
Point(7) = {0,2,2,h};
Point(8) = {2,2,2,h};




//
// // FACES DEFINITION
//
// Curve Loop(15) = {-10,6,2,-9};
// Surface(15) = {15};
// Curve Loop(17) = {12,7,-4,-8};
// Surface(17) = {17};
// Curve Loop(28) = {36,23,4,20};
// Surface(28) = {28};
// Curve Loop(29) = {-36,21,2,22};
// Surface(29) = {29};
// Curve Loop(30) = {-37,-24,7,20};
// Surface(30) = {30};
// Curve Loop(31) = {37,21,-6,-25};
// Surface(31) = {31};
// Curve Loop(32) = {-38,27,8,-23};
// Surface(32) = {32};
// Curve Loop(33) = {38,-22,-9,26};
// Surface(33) = {33};
// Curve Loop(34) = {-39,27,12,24};
// Surface(34) = {34};
// Curve Loop(35) = {39,25,10,26};
// Surface(35) = {35};
// Curve Loop(40) = {-36,-37,-39,38};
// Surface(40) = {40};
//
//
//
// // VOLUMES DEFINITION
//
// Surface Loop(3) = {15,29,31,33,35,40};
// Volume(3) = {3};
// Surface Loop(4) = {17,28,30,32,34,40};
// Volume(4) = {4};
//
//
// // COARSE ELEMENTS GROUPING
//
// Physical Volume("c1",1) = {3,4};
//
//
//
//  // FRACTURES
//
// frac1[] = {40};
//
// Physical Surface("Fracture2",2) = {frac1[]};
//
//
//
// // BOUNDARY CONDITIONS
//
// // Physical Surface("bc1",1) = {15,17,28,29,30,31,32,33,34,35};
// // Physical Surface("inlet") = {31,30};
// // Physical Surface("outlet") = {32,33};
// // Physical Surface("noflux") = {28,29,15,35,34,17};
// Physical Surface("inlet") = {17};
// Physical Surface("outlet") = {15};
// Physical Surface("noflux") = {29,28,31,30,35,34,33,32};
//
//
// BCfrac0[] = {36,37,39,38};
// Physical Curve("BCfrac0",10) = {BCfrac0[]};
//
// // OPTIONS
//
// Coherence Mesh;
// Transfinite Curve{:} = 2;
// Transfinite Surface{:};
// Transfinite Volume{:};
// Recombine Surface{:};
// Recombine Volume{:};//+


// //+
// Line(1) = {3, 4};
// //+
// Line(2) = {4, 8};
// //+
// Line(3) = {8, 7};
// //+
// Line(4) = {7, 3};
// //+
// Line(5) = {7, 5};
// //+
// Line(6) = {5, 6};
// //+
// Line(7) = {6, 8};
// //+
// Line(8) = {6, 2};
// //+
// Line(9) = {2, 4};
// //+
// Line(10) = {1, 3};
// //+
// Line(11) = {5, 1};
// //+
// Line(12) = {1, 2};
// //+
// Curve Loop(1) = {3, 5, 6, 7};
// //+
// Plane Surface(1) = {1};
// //+
// Curve Loop(2) = {7, -2, -9, -8};
// //+
// Plane Surface(2) = {2};
// //+
// Curve Loop(3) = {9, -1, -10, 12};
// //+
// Plane Surface(3) = {3};
// //+
// Curve Loop(4) = {4, -10, -11, -5};
// //+
// Plane Surface(4) = {4};
// //+
// Curve Loop(5) = {1, 2, 3, 4};
// //+
// Plane Surface(5) = {5};
// //+
// Curve Loop(6) = {12, -8, -6, 11};
// //+
// Plane Surface(6) = {6};
// //+
// Surface Loop(1) = {4, 5, 3, 2, 1, 6};
// //+
// Volume(1) = {1};
// Coherence Mesh;
// Transfinite Curve{:} = 2;
// Transfinite Surface{:};
// Transfinite Volume{:};
// Recombine Surface{:};
// Recombine Volume{:};//+
//+
Line(1) = {3, 4};
//+
Line(2) = {4, 8};
//+
Line(3) = {7, 8};
//+
Line(4) = {3, 7};
//+
Line(5) = {5, 7};
//+
Line(6) = {5, 6};
//+
Line(7) = {6, 8};
//+
Line(8) = {1, 3};
//+
Line(9) = {2, 4};
//+
Line(10) = {1, 2};
//+
Line(11) = {1, 5};
//+
Line(12) = {2, 6};
//+
Curve Loop(1) = {1, 2, -3, -4};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {3, -7, -6, 5};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {10, 12, -6, -11};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {1, -9, -10, 8};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {2, -7, -12, 9};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {4, -5, -11, 8};
//+
Plane Surface(6) = {6};
//+
Physical Surface("inlet") = {2};
//+
Physical Surface("outlet") = {4};
//+
Physical Surface("noflux") = {1, 5, 6, 3};
//+
Surface Loop(1) = {6, 1, 4, 5, 2, 3};
//+
Volume(1) = {1};
//+
Physical Volume("k33") = {1};

Coherence Mesh;
Transfinite Curve{:} = 10;
Transfinite Surface{:};
Transfinite Volume{:};
Recombine Surface{:};
Recombine Volume{:};//+
