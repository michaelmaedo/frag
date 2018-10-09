// Gmsh project created on Mon Sep  3 18:18:27 2018
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Point(5) = {0.2, 1, 0, 1.0};
//+
Point(6) = {0.8, 1, 0, 1.0};
//+
Point(7) = {0.8, .2, 0, 1.0};
//+
Point(8) = {0.2, .2, 0, 1.0};
//+
Point(9) = {0, 0.2, 0, 1.0};
//+
Point(10) = {0.2, 0., 0, 1.0};
//+
Point(11) = {0.8, 0., 0, 1.0};
//+
Point(12) = {1, 0.2, 0, 1.0};
//+
Line(1) = {1, 10};
//+
Line(2) = {10, 11};
//+
Line(3) = {11, 2};
//+
Line(4) = {9, 8};
//+
Line(5) = {8, 7};
//+
Line(6) = {7, 12};
//+
Line(7) = {4, 5};
//+
Line(8) = {5, 6};
//+
Line(9) = {6, 3};
//+
Line(10) = {1, 9};
//+
Line(11) = {9, 4};
//+
Line(12) = {10, 8};
//+
Line(13) = {8, 5};
//+
Line(14) = {11, 7};
//+
Line(15) = {7, 6};
//+
Line(16) = {2, 12};
//+
Line(17) = {12, 3};
//+
Line Loop(1) = {13, -7, -11, 4};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {8, -15, -5, 13};
//+
Plane Surface(2) = {2};
//+
Line Loop(3) = {15, 9, -17, -6};
//+
Plane Surface(3) = {3};
//+
Line Loop(4) = {4, -12, -1, 10};
//+
Plane Surface(4) = {4};
//+
Line Loop(5) = {2, 14, -5, -12};
//+
Plane Surface(5) = {5};
//+
Line Loop(6) = {6, -16, -3, 14};
//+
Plane Surface(6) = {6};
//+
Physical Line("ux=0_uy=0.MBC") = {2, 3, 1};
//+
Physical Line("ux=0.MBC") = {11, 10, 17, 16};
//+
Physical Line("fy=-5.MBC") = {8};
//+
Physical Line("suction.FBC") = {8};
//+
Physical Surface("matfrag_h=0.001.dam") = {2};
//+
Physical Surface("mat2.dam") = {1, 3, 4, 5, 6};
//+
Transfinite Surface {1};
//+
Transfinite Surface {4};
//+
Transfinite Surface {5};
//+
Transfinite Surface {6};
//+
Transfinite Surface {3};
//+
Transfinite Surface {2} Alternated;
//+
Transfinite Line {11, 13, 15, 17} = 5 Using Progression 1;
//+
Transfinite Line {10, 10, 12, 14, 16} = 2 Using Progression 1;
//+
Transfinite Line {1, 4, 7, 3, 6, 9} = 2 Using Progression 1;
//+
Transfinite Line {8, 5, 2} = 4 Using Progression 1;
