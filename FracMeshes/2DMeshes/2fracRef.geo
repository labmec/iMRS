//  Geo file generated by DFNMesh project
//	https://github.com/labmec/dfnMesh// POINTS DEFINITION 

h = 0;

Point(1) = {0,0,0, h};
Point(2) = {2,0,0, h};
Point(3) = {0,1,0, h};
Point(4) = {2,1,0, h};
Point(5) = {0,2,0, h};
Point(6) = {2,2,0, h};
Point(7) = {0,0,2, h};
Point(8) = {2,0,2, h};
Point(9) = {0,1,2, h};
Point(10) = {2,1,2, h};
Point(11) = {0,2,2, h};
Point(12) = {2,2,2, h};
Point(13) = {0,1,1, h};
Point(14) = {2,1,1, h};
Point(15) = {2,0,1, h};
Point(16) = {0,0,1, h};
Point(17) = {0,2,1, h};
Point(18) = {2,2,1, h};


// LINES DEFINITION 

Line(3) = {9,7};
Line(4) = {9,10};
Line(5) = {8,10};
Line(6) = {7,8};
Line(11) = {3,1};
Line(12) = {3,4};
Line(13) = {2,4};
Line(14) = {1,2};
Line(15) = {11,9};
Line(16) = {12,11};
Line(17) = {10,12};
Line(20) = {5,3};
Line(21) = {6,5};
Line(22) = {4,6};
Line(34) = {3,13};
Line(35) = {13,9};
Line(36) = {4,14};
Line(37) = {14,10};
Line(38) = {2,15};
Line(39) = {15,8};
Line(40) = {1,16};
Line(41) = {16,7};
Line(42) = {5,17};
Line(43) = {17,11};
Line(44) = {6,18};
Line(45) = {18,12};
Line(60) = {16,13};
Line(61) = {14,13};
Line(62) = {14,15};
Line(63) = {15,16};
Line(64) = {17,13};
Line(65) = {18,17};
Line(66) = {18,14};


// FACES DEFINITION 

Curve Loop(23) = {6,5,-4,3};
Surface(23) = {23};
Curve Loop(28) = {14,13,-12,11};
Surface(28) = {28};
Curve Loop(29) = {4,17,16,15};
Surface(29) = {29};
Curve Loop(33) = {12,22,21,20};
Surface(33) = {33};
Curve Loop(46) = {60,35,3,-41};
Surface(46) = {46};
Curve Loop(47) = {-60,-40,-11,34};
Surface(47) = {47};
Curve Loop(48) = {-61,37,-4,-35};
Surface(48) = {48};
Curve Loop(49) = {61,-34,12,36};
Surface(49) = {49};
Curve Loop(50) = {-62,37,-5,-39};
Surface(50) = {50};
Curve Loop(51) = {62,-38,13,36};
Surface(51) = {51};
Curve Loop(52) = {-63,39,-6,-41};
Surface(52) = {52};
Curve Loop(53) = {63,-40,14,38};
Surface(53) = {53};
Curve Loop(54) = {-64,43,15,-35};
Surface(54) = {54};
Curve Loop(55) = {64,-34,-20,42};
Surface(55) = {55};
Curve Loop(56) = {-65,45,16,-43};
Surface(56) = {56};
Curve Loop(57) = {65,-42,-21,44};
Surface(57) = {57};
Curve Loop(58) = {-66,45,-17,-37};
Surface(58) = {58};
Curve Loop(59) = {66,-36,22,44};
Surface(59) = {59};
Curve Loop(67) = {-60,-63,-62,61};
Surface(67) = {67};
Curve Loop(68) = {61,-64,-65,66};
Surface(68) = {68};



// VOLUMES DEFINITION 

Surface Loop(4) = {23,46,48,50,52,67};
Volume(4) = {4};
Surface Loop(5) = {28,47,49,51,53,67};
Volume(5) = {5};
Surface Loop(6) = {33,49,55,57,59,68};
Volume(6) = {6};
Surface Loop(7) = {29,48,54,56,58,68};
Volume(7) = {7};


// COARSE ELEMENTS GROUPING

Physical Volume("c1",1) = {4,5};
Physical Volume("c2",2) = {6,7};



 // FRACTURES

frac1[] = {67,68};
frac2[] = {48,49};

Physical Surface("Fracture12",12) = {frac1[], frac2[]};



// BOUNDARY CONDITIONS

Physical Surface("bc1",1) = {23,28,29,33,46,47,50,51,52,53,54,55,56,57,58,59};


BCfrac0_0[] = { 63};
Physical Curve("BCfrac0_0", 101) = {BCfrac0_0[]};

BCfrac0_1[] = { 62,66};
Physical Curve("BCfrac0_1", 102) = {BCfrac0_1[]};

BCfrac0_2[] = { 65};
Physical Curve("BCfrac0_2", 103) = {BCfrac0_2[]};

BCfrac0_3[] = { 64,60};
Physical Curve("BCfrac0_3", 104) = {BCfrac0_3[]};



BCfrac1_0[] = { 37,36};
Physical Curve("BCfrac1_0", 201) = {BCfrac1_0[]};

BCfrac1_1[] = { 12};
Physical Curve("BCfrac1_1", 202) = {BCfrac1_1[]};

BCfrac1_2[] = { 34,35};
Physical Curve("BCfrac1_2", 203) = {BCfrac1_2[]};

BCfrac1_3[] = { 4};
Physical Curve("BCfrac1_3", 204) = {BCfrac1_3[]};



// INTER-FRACTURE INTERSECTIONS

fracIntersection_0_1[] = { 61};
fracIntersection_1_0[] = fracIntersection_0_1;
Physical Curve("fracIntersection_0_1") = fracIntersection_0_1[];


// OPTIONS

Coherence Mesh;
 Transfinite Curve {:} = 2;
 Transfinite Surface{:};
 Transfinite Volume{:};
 Recombine Surface{:};
 Recombine Volume{:};