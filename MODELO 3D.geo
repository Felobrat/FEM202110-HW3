// Gmsh project created on Thu May 06 09:56:02 2021
SetFactory("OpenCASCADE");

a = 2;
b = 4;
c = 16;
n = 1.0;
e = 0.45;
f = 0.50;
r = 1;
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, b, 0, 1.0};
//+
Point(3) = {a, b, 0, 1.0};
//+
Point(4) = {a, 0, 0, 1.0};
//+
Point(5) = {a, 0, f-e, 1.0};
//+
Point(6) = {a+c, 0, f-e, 1.0};
//+
Point(7) = {a+c, b, f-e, 1.0};
//+
Point(8) = {a, b, f-e, 1.0};
//+
Point(9) = {a+c, 0, 0, 1.0};
//+
Point(10) = {(2*a+c), 0, 0, 1.0};
//+
Point(11) = {(a+c), b, 0, 1.0};
//+
Point(12) = {(2*a+c), b, 0, 1.0};
//+
Line(1) = {1, 4};
//+
Line(2) = {4, 3};
//+
Line(3) = {3, 2};
//+
Line(4) = {2, 1};
//+
Line(5) = {8, 5};
//+
Line(6) = {5, 6};
//+
Line(7) = {6, 7};
//+
Line(8) = {7, 8};
//+
Line(9) = {9, 10};
//+
Line(10) = {10, 12};
//+
Line(11) = {12, 11};
//+
Line(12) = {11, 9};
//+
Point(13) = {(2*a+c)/2, b/2, f-e, 1.0};
//+
Circle(13) = {(2*a+c)/2, b/2, f-e, r, 0, 2*Pi};
//+
Curve Loop(1) = {6, 7, 8, 5};
//+
Curve Loop(2) = {13};
//+
Plane Surface(1) = {1, 2};
//+
Curve Loop(3) = {9, 10, 11, 12};
//+
Plane Surface(2) = {3};
//+
Curve Loop(4) = {1, 2, 3, 4};
//+
Plane Surface(3) = {4};
//+
//+
Extrude {0, 0, f} {
  Curve{9}; Curve{10}; Curve{11}; Curve{12}; Curve{1}; Curve{4}; Curve{3}; Curve{2}; 
}
//+
Extrude {0, 0, 0.4} {
  Curve{6}; Curve{5}; Curve{8}; Curve{7}; Curve{13}; 
}//+
Curve Loop(18) = {26, 24, 29, 28};
//+
Plane Surface(17) = {18};
//+
Curve Loop(19) = {20, 21, 16, 18};
//+
Plane Surface(18) = {19};
//+
//+
Curve Loop(20) = {36, 34, 32, 37};
//+
Curve Loop(21) = {39};
//+
Plane Surface(19) = {20, 21};
//+
Surface Loop(1) = {3, 8, 9, 10, 11, 17};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {18, 6, 5, 4, 7, 2};
//+
Volume(2) = {2};
//+
Surface Loop(3) = {19, 14, 15, 12, 13, 1, 16};
//+
Volume(3) = {3};
//+
Physical Volume(40) = {1};
//+
Physical Volume(41) = {3};
//+
Physical Volume(42) = {2};
//+
Transfinite Curve {39, 13} = 40 Using Progression 1;
