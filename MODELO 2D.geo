// Gmsh project created on Thu May 06 09:56:02 2021
SetFactory("OpenCASCADE");

a = 2;
b = 4;
c = 16;
n = 1.0;
e = 4.5;
f = 5.0;

//+
Point(1) = {0, 0, 0, n};
//+
Point(2) = {0, 4, 0, n};
//+
Point(3) = {a, 4, 0, n};
//+
Point(4) = {a, 0, 0, n};
//+
Point(5) = {a, 0, 0.0, n};
//+
Point(6) = {a+c, 0, 0.0, n};
//+
Point(7) = {a+c, b, 0.0, n};
//+
Point(8) = {a, b, 0.0, n};
//+
Point(9) = {a+c, 0, 0, n};
//+
Point(10) = {2*a+c, 0, 0, n};
//+
Point(11) = {a+c, b, 0, n};
//+
Point(12) = {2*a+c, b, 0, n};
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
Point(13) = {(2*a+c)/2, b/2, 0.0, 1.0};
//+
Circle(13) = {(2*a+c)/2, b/2, 0.0, 1, 0, 2*Pi};
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
Physical Curve("Fixed Displacements", 14) = {4};
//+
Physical Surface("Empotrado", 15) = {3};
//+
Physical Surface("Placa", 16) = {1};
//+
Physical Surface("Hook", 17) = {2};
//+
Physical Curve("Stress", 18) = {10};
//+
Transfinite Curve {13} = 40 Using Progression 1;
