ht=250;

//Y-coordinate of the points of the two internal lines
yval1=1.5;
yval2=0.75;

//X-coordinate of the points of the two internal lines
valx1= 0.3333;
valx2= 0.6666;

//Number of elements in each partition of x
nelx1=1;
nelx2=2;
nelx3=3;

//Number of elements in each partition of y
nely1=1;
nely2=2;
nely3=3;

//Number of elements in each partition of z
nelz1=1;
nelz2=2;
nelz3=3;



Point(1)= {0., 2.25, 0., ht};
Point(2)= {valx1, 2.25, 0., ht};
Point(3)= {valx2, 2.25, 0., ht};
Point(4)= {1., 2.25, 0., ht};
Point(5)= {0., yval1, 0., ht};
Point(6)= {valx1, yval1, 0., ht};
Point(7)= {valx2, yval1, 0., ht};
Point(8)= {1., yval1, 0., ht};
Point(9)= {0., yval2, 0., ht};
Point(10)= {valx1, yval2, 0., ht};
Point(11)= {valx2, yval2, 0., ht};
Point(12)= {1., yval2, 0., ht};
Point(13)= {0., 0., 0., ht};
Point(14)= {valx1, 0., 0., ht};
Point(15)= {valx2, 0., 0., ht};
Point(16)= {1., 0., 0., ht};
Point(17)= {0., 2.25, 0.333333, ht};
Point(18)= {valx1, 2.25, 0.333333, ht};
Point(19)= {valx2, 2.25, 0.333333, ht};
Point(20)= {1., 2.25, 0.333333, ht};
Point(21)= {0., yval1, 0.333333, ht};
Point(22)= {valx1, yval1, 0.333333, ht};
Point(23)= {valx2, yval1, 0.333333, ht};
Point(24)= {1., yval1, 0.333333, ht};
Point(25)= {0., yval2, 0.333333, ht};
Point(26)= {valx1, yval2, 0.333333, ht};
Point(27)= {valx2, yval2, 0.333333, ht};
Point(28)= {1., yval2, 0.333333, ht};
Point(29)= {0., 0., 0.333333, ht};
Point(30)= {valx1, 0., 0.333333, ht};
Point(31)= {valx2, 0., 0.333333, ht};
Point(32)= {1., 0., 0.333333, ht};
Point(33)= {0., 2.25, 0.666667, ht};
Point(34)= {valx1, 2.25, 0.666667, ht};
Point(35)= {valx2, 2.25, 0.666667, ht};
Point(36)= {1., 2.25, 0.666667, ht};
Point(37)= {0., yval1, 0.666667, ht};
Point(38)= {valx1, yval1, 0.666667, ht};
Point(39)= {valx2, yval1, 0.666667, ht};
Point(40)= {1., yval1, 0.666667, ht};
Point(41)= {0., yval2, 0.666667, ht};
Point(42)= {valx1, yval2, 0.666667, ht};
Point(43)= {valx2, yval2, 0.666667, ht};
Point(44)= {1., yval2, 0.666667, ht};
Point(45)= {0., 0., 0.666667, ht};
Point(46)= {valx1, 0., 0.666667, ht};
Point(47)= {valx2, 0., 0.666667, ht};
Point(48)= {1., 0., 0.666667, ht};
Point(49)= {0., 2.25, 1., ht};
Point(50)= {valx1, 2.25, 1., ht};
Point(51)= {valx2, 2.25, 1., ht};
Point(52)= {1., 2.25, 1., ht};
Point(53)= {0., yval1, 1., ht};
Point(54)= {valx1, yval1, 1., ht};
Point(55)= {valx2, yval1, 1., ht};
Point(56)= {1., yval1, 1., ht};
Point(57)= {0., yval2, 1., ht};
Point(58)= {valx1, yval2, 1., ht};
Point(59)= {valx2, yval2, 1., ht};
Point(60)= {1., yval2, 1., ht};
Point(61)= {0., 0., 1., ht};
Point(62)= {valx1, 0., 1., ht};
Point(63)= {valx2, 0., 1., ht};
Point(64)= {1., 0., 1., ht};
Line(1)= {1, 2};
Line(2)= {2, 6};
Line(3)= {6, 5};
Line(4)= {5, 1};
Line(5)= {2, 3};
Line(6)= {3, 7};
Line(7)= {7, 6};
Line(8)= {3, 4};
Line(9)= {4, 8};
Line(10)= {8, 7};
Line(11)= {6, 10};
Line(12)= {10, 9};
Line(13)= {9, 5};
Line(14)= {7, 11};
Line(15)= {11, 10};
Line(16)= {8, 12};
Line(17)= {12, 11};
Line(18)= {10, 14};
Line(19)= {14, 13};
Line(20)= {13, 9};
Line(21)= {11, 15};
Line(22)= {15, 14};
Line(23)= {12, 16};
Line(24)= {16, 15};
Line(25)= {17, 18};
Line(26)= {18, 22};
Line(27)= {22, 21};
Line(28)= {21, 17};
Line(29)= {18, 19};
Line(30)= {19, 23};
Line(31)= {23, 22};
Line(32)= {19, 20};
Line(33)= {20, 24};
Line(34)= {24, 23};
Line(35)= {22, 26};
Line(36)= {26, 25};
Line(37)= {25, 21};
Line(38)= {23, 27};
Line(39)= {27, 26};
Line(40)= {24, 28};
Line(41)= {28, 27};
Line(42)= {26, 30};
Line(43)= {30, 29};
Line(44)= {29, 25};
Line(45)= {27, 31};
Line(46)= {31, 30};
Line(47)= {28, 32};
Line(48)= {32, 31};
Line(49)= {33, 34};
Line(50)= {34, 38};
Line(51)= {38, 37};
Line(52)= {37, 33};
Line(53)= {34, 35};
Line(54)= {35, 39};
Line(55)= {39, 38};
Line(56)= {35, 36};
Line(57)= {36, 40};
Line(58)= {40, 39};
Line(59)= {38, 42};
Line(60)= {42, 41};
Line(61)= {41, 37};
Line(62)= {39, 43};
Line(63)= {43, 42};
Line(64)= {40, 44};
Line(65)= {44, 43};
Line(66)= {42, 46};
Line(67)= {46, 45};
Line(68)= {45, 41};
Line(69)= {43, 47};
Line(70)= {47, 46};
Line(71)= {44, 48};
Line(72)= {48, 47};
Line(73)= {49, 50};
Line(74)= {50, 54};
Line(75)= {54, 53};
Line(76)= {53, 49};
Line(77)= {50, 51};
Line(78)= {51, 55};
Line(79)= {55, 54};
Line(80)= {51, 52};
Line(81)= {52, 56};
Line(82)= {56, 55};
Line(83)= {54, 58};
Line(84)= {58, 57};
Line(85)= {57, 53};
Line(86)= {55, 59};
Line(87)= {59, 58};
Line(88)= {56, 60};
Line(89)= {60, 59};
Line(90)= {58, 62};
Line(91)= {62, 61};
Line(92)= {61, 57};
Line(93)= {59, 63};
Line(94)= {63, 62};
Line(95)= {60, 64};
Line(96)= {64, 63};
Line Loop(1)= {1, 2,3,4};
Surface(1)= {1};
Line Loop(2)= {5, 6,7,-2};
Surface(2)= {2};
Line Loop(3)= {8, 9,10,-6};
Surface(3)= {3};
Line Loop(4)= {-3, 11,12,13};
Surface(4)= {4};
Line Loop(5)= {-7, 14,15,-11};
Surface(5)= {5};
Line Loop(6)= {-10, 16,17,-14};
Surface(6)= {6};
Line Loop(7)= {-12, 18,19,20};
Surface(7)= {7};
Line Loop(8)= {-15, 21,22,-18};
Surface(8)= {8};
Line Loop(9)= {-17, 23,24,-21};
Surface(9)= {9};
Line Loop(10)= {25, 26,27,28};
Surface(10)= {10};
Line Loop(11)= {29, 30,31,-26};
Surface(11)= {11};
Line Loop(12)= {32, 33,34,-30};
Surface(12)= {12};
Line Loop(13)= {-27, 35,36,37};
Surface(13)= {13};
Line Loop(14)= {-31, 38,39,-35};
Surface(14)= {14};
Line Loop(15)= {-34, 40,41,-38};
Surface(15)= {15};
Line Loop(16)= {-36, 42,43,44};
Surface(16)= {16};
Line Loop(17)= {-39, 45,46,-42};
Surface(17)= {17};
Line Loop(18)= {-41, 47,48,-45};
Surface(18)= {18};
Line Loop(19)= {49, 50,51,52};
Surface(19)= {19};
Line Loop(20)= {53, 54,55,-50};
Surface(20)= {20};
Line Loop(21)= {56, 57,58,-54};
Surface(21)= {21};
Line Loop(22)= {-51, 59,60,61};
Surface(22)= {22};
Line Loop(23)= {-55, 62,63,-59};
Surface(23)= {23};
Line Loop(24)= {-58, 64,65,-62};
Surface(24)= {24};
Line Loop(25)= {-60, 66,67,68};
Surface(25)= {25};
Line Loop(26)= {-63, 69,70,-66};
Surface(26)= {26};
Line Loop(27)= {-65, 71,72,-69};
Surface(27)= {27};
Line Loop(28)= {73, 74,75,76};
Surface(28)= {28};
Line Loop(29)= {77, 78,79,-74};
Surface(29)= {29};
Line Loop(30)= {80, 81,82,-78};
Surface(30)= {30};
Line Loop(31)= {-75, 83,84,85};
Surface(31)= {31};
Line Loop(32)= {-79, 86,87,-83};
Surface(32)= {32};
Line Loop(33)= {-82, 88,89,-86};
Surface(33)= {33};
Line Loop(34)= {-84, 90,91,92};
Surface(34)= {34};
Line Loop(35)= {-87, 93,94,-90};
Surface(35)= {35};
Line Loop(36)= {-89, 95,96,-93};
Surface(36)= {36};
Line(97)= {2, 18};
Line(98)= {17, 1};
Line(99)= {6, 22};
Line(100)= {5, 21};
Line(101)= {3, 19};
Line(102)= {7, 23};
Line(103)= {4, 20};
Line(104)= {8, 24};
Line(105)= {10, 26};
Line(106)= {9, 25};
Line(107)= {11, 27};
Line(108)= {12, 28};
Line(109)= {14, 30};
Line(110)= {13, 29};
Line(111)= {15, 31};
Line(112)= {16, 32};
Line(113)= {18, 34};
Line(114)= {33, 17};
Line(115)= {22, 38};
Line(116)= {21, 37};
Line(117)= {19, 35};
Line(118)= {23, 39};
Line(119)= {20, 36};
Line(120)= {24, 40};
Line(121)= {26, 42};
Line(122)= {25, 41};
Line(123)= {27, 43};
Line(124)= {28, 44};
Line(125)= {30, 46};
Line(126)= {29, 45};
Line(127)= {31, 47};
Line(128)= {32, 48};
Line(129)= {34, 50};
Line(130)= {49, 33};
Line(131)= {38, 54};
Line(132)= {37, 53};
Line(133)= {35, 51};
Line(134)= {39, 55};
Line(135)= {36, 52};
Line(136)= {40, 56};
Line(137)= {42, 58};
Line(138)= {41, 57};
Line(139)= {43, 59};
Line(140)= {44, 60};
Line(141)= {46, 62};
Line(142)= {45, 61};
Line(143)= {47, 63};
Line(144)= {48, 64};
Line Loop(37)= {1, 97,-25,98};
Surface(37)= {37};
Line Loop(38)= {2, 99,-26,-97};
Surface(38)= {38};
Line Loop(39)= {3, 100,-27,-99};
Surface(39)= {39};
Line Loop(40)= {4, -98,-28,-100};
Surface(40)= {40};
Line Loop(41)= {5, 101,-29,-97};
Surface(41)= {41};
Line Loop(42)= {6, 102,-30,-101};
Surface(42)= {42};
Line Loop(43)= {7, 99,-31,-102};
Surface(43)= {43};
Line Loop(44)= {8, 103,-32,-101};
Surface(44)= {44};
Line Loop(45)= {9, 104,-33,-103};
Surface(45)= {45};
Line Loop(46)= {10, 102,-34,-104};
Surface(46)= {46};
Line Loop(47)= {11, 105,-35,-99};
Surface(47)= {47};
Line Loop(48)= {12, 106,-36,-105};
Surface(48)= {48};
Line Loop(49)= {13, 100,-37,-106};
Surface(49)= {49};
Line Loop(50)= {14, 107,-38,-102};
Surface(50)= {50};
Line Loop(51)= {15, 105,-39,-107};
Surface(51)= {51};
Line Loop(52)= {16, 108,-40,-104};
Surface(52)= {52};
Line Loop(53)= {17, 107,-41,-108};
Surface(53)= {53};
Line Loop(54)= {18, 109,-42,-105};
Surface(54)= {54};
Line Loop(55)= {19, 110,-43,-109};
Surface(55)= {55};
Line Loop(56)= {20, 106,-44,-110};
Surface(56)= {56};
Line Loop(57)= {21, 111,-45,-107};
Surface(57)= {57};
Line Loop(58)= {22, 109,-46,-111};
Surface(58)= {58};
Line Loop(59)= {23, 112,-47,-108};
Surface(59)= {59};
Line Loop(60)= {24, 111,-48,-112};
Surface(60)= {60};
Line Loop(61)= {25, 113,-49,114};
Surface(61)= {61};
Line Loop(62)= {26, 115,-50,-113};
Surface(62)= {62};
Line Loop(63)= {27, 116,-51,-115};
Surface(63)= {63};
Line Loop(64)= {28, -114,-52,-116};
Surface(64)= {64};
Line Loop(65)= {29, 117,-53,-113};
Surface(65)= {65};
Line Loop(66)= {30, 118,-54,-117};
Surface(66)= {66};
Line Loop(67)= {31, 115,-55,-118};
Surface(67)= {67};
Line Loop(68)= {32, 119,-56,-117};
Surface(68)= {68};
Line Loop(69)= {33, 120,-57,-119};
Surface(69)= {69};
Line Loop(70)= {34, 118,-58,-120};
Surface(70)= {70};
Line Loop(71)= {35, 121,-59,-115};
Surface(71)= {71};
Line Loop(72)= {36, 122,-60,-121};
Surface(72)= {72};
Line Loop(73)= {37, 116,-61,-122};
Surface(73)= {73};
Line Loop(74)= {38, 123,-62,-118};
Surface(74)= {74};
Line Loop(75)= {39, 121,-63,-123};
Surface(75)= {75};
Line Loop(76)= {40, 124,-64,-120};
Surface(76)= {76};
Line Loop(77)= {41, 123,-65,-124};
Surface(77)= {77};
Line Loop(78)= {42, 125,-66,-121};
Surface(78)= {78};
Line Loop(79)= {43, 126,-67,-125};
Surface(79)= {79};
Line Loop(80)= {44, 122,-68,-126};
Surface(80)= {80};
Line Loop(81)= {45, 127,-69,-123};
Surface(81)= {81};
Line Loop(82)= {46, 125,-70,-127};
Surface(82)= {82};
Line Loop(83)= {47, 128,-71,-124};
Surface(83)= {83};
Line Loop(84)= {48, 127,-72,-128};
Surface(84)= {84};
Line Loop(85)= {49, 129,-73,130};
Surface(85)= {85};
Line Loop(86)= {50, 131,-74,-129};
Surface(86)= {86};
Line Loop(87)= {51, 132,-75,-131};
Surface(87)= {87};
Line Loop(88)= {52, -130,-76,-132};
Surface(88)= {88};
Line Loop(89)= {53, 133,-77,-129};
Surface(89)= {89};
Line Loop(90)= {54, 134,-78,-133};
Surface(90)= {90};
Line Loop(91)= {55, 131,-79,-134};
Surface(91)= {91};
Line Loop(92)= {56, 135,-80,-133};
Surface(92)= {92};
Line Loop(93)= {57, 136,-81,-135};
Surface(93)= {93};
Line Loop(94)= {58, 134,-82,-136};
Surface(94)= {94};
Line Loop(95)= {59, 137,-83,-131};
Surface(95)= {95};
Line Loop(96)= {60, 138,-84,-137};
Surface(96)= {96};
Line Loop(97)= {61, 132,-85,-138};
Surface(97)= {97};
Line Loop(98)= {62, 139,-86,-134};
Surface(98)= {98};
Line Loop(99)= {63, 137,-87,-139};
Surface(99)= {99};
Line Loop(100)= {64, 140,-88,-136};
Surface(100)= {100};
Line Loop(101)= {65, 139,-89,-140};
Surface(101)= {101};
Line Loop(102)= {66, 141,-90,-137};
Surface(102)= {102};
Line Loop(103)= {67, 142,-91,-141};
Surface(103)= {103};
Line Loop(104)= {68, 138,-92,-142};
Surface(104)= {104};
Line Loop(105)= {69, 143,-93,-139};
Surface(105)= {105};
Line Loop(106)= {70, 141,-94,-143};
Surface(106)= {106};
Line Loop(107)= {71, 144,-95,-140};
Surface(107)= {107};
Line Loop(108)= {72, 143,-96,-144};
Surface(108)= {108};
Surface Loop(1001)= {1, 10,37,38,39,40};
Volume(1001)= {1001};
Surface Loop(1002)= {2, 11,41,42,43,38};
Volume(1002)= {1002};
Surface Loop(1003)= {3, 12,44,45,46,42};
Volume(1003)= {1003};
Surface Loop(1004)= {4, 13,39,47,48,49};
Volume(1004)= {1004};
Surface Loop(1005)= {5, 14,43,50,51,47};
Volume(1005)= {1005};
Surface Loop(1006)= {6, 15,46,52,53,50};
Volume(1006)= {1006};
Surface Loop(1007)= {7, 16,48,54,55,56};
Volume(1007)= {1007};
Surface Loop(1008)= {8, 17,51,57,58,54};
Volume(1008)= {1008};
Surface Loop(1009)= {9, 18,53,59,60,57};
Volume(1009)= {1009};
Surface Loop(1010)= {10, 19,61,62,63,64};
Volume(1010)= {1010};
Surface Loop(1011)= {11, 20,65,66,67,62};
Volume(1011)= {1011};
Surface Loop(1012)= {12, 21,68,69,70,66};
Volume(1012)= {1012};
Surface Loop(1013)= {13, 22,63,71,72,73};
Volume(1013)= {1013};
Surface Loop(1014)= {14, 23,67,74,75,71};
Volume(1014)= {1014};
Surface Loop(1015)= {15, 24,70,76,77,74};
Volume(1015)= {1015};
Surface Loop(1016)= {16, 25,72,78,79,80};
Volume(1016)= {1016};
Surface Loop(1017)= {17, 26,75,81,82,78};
Volume(1017)= {1017};
Surface Loop(1018)= {18, 27,77,83,84,81};
Volume(1018)= {1018};
Surface Loop(1019)= {19, 28,85,86,87,88};
Volume(1019)= {1019};
Surface Loop(1020)= {20, 29,89,90,91,86};
Volume(1020)= {1020};
Surface Loop(1021)= {21, 30,92,93,94,90};
Volume(1021)= {1021};
Surface Loop(1022)= {22, 31,87,95,96,97};
Volume(1022)= {1022};
Surface Loop(1023)= {23, 32,91,98,99,95};
Volume(1023)= {1023};
Surface Loop(1024)= {24, 33,94,100,101,98};
Volume(1024)= {1024};
Surface Loop(1025)= {25, 34,96,102,103,104};
Volume(1025)= {1025};
Surface Loop(1026)= {26, 35,99,105,106,102};
Volume(1026)= {1026};
Surface Loop(1027)= {27, 36,101,107,108,105};
Volume(1027)= {1027};

Physical Volume("k33") = Volume{:};
Physical Surface("inlet") = {79, 82, 84};
Physical Surface("outlet") = {85, 89, 92, 37, 41, 44};
Physical Surface("noflux")= {1, 40, 2, 3, 45, 4, 49, 5, 6, 52, 7, 55, 56, 8, 58, 9, 60, 59, 61, 64, 65, 68, 69, 73, 76, 80, 83, 28, 88, 29,30, 93, 31, 97, 32, 33, 100, 34, 103, 104, 35, 106, 36, 108, 107};


Transfinite Curve{:}=2;
Transfinite Curve{20,44,68,92,18,42,66,90,21,45,69,93,23,47,71,95}=nely1+1;
Transfinite Curve{13,37,61,85,11,35,59,83,14,38,62,86,16,40,64,88}=nely2+1;
Transfinite Curve{4,28,52,76,2,26,50,74,6,30,54,78,9,33,57,81}=nely3+1;

Transfinite Curve{98,100,106,110,97,99,105,109,101,102,107,111,103,104,108,112}=nelx1+1;
Transfinite Curve{114,116, 122, 126,113,115,121,125,117,118,123,127,119,120,124,128}=nelx2+1;
Transfinite Curve{130,132,138,142,129,131,137,141,133,134,139,143,135,136,140,144}=nelx3+1;

Transfinite Curve{19,43,67,91,12,36,60,84,3,27,51,75,1,25,49,73}=nelz1+1;
Transfinite Curve{22,46,70,94,15,39,63,87,7,31,55,79,5,29,53,77}=nelz2+1;
Transfinite Curve{24,48,72,96,17,41,65,89,10,34,58,82,8,32,56,80}=nelz3+1;

Transfinite Surface{:};
Transfinite Volume{:};
Recombine Surface{:};
Recombine Volume{:};//+

