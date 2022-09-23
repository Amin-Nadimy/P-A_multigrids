Point(1) = {0, 0, 0};
Point(2) = {1, 0, 0};
Point(3) = {1, 1, 0};
Point(4) = {0, 1, 0};
Point(5) = {0.75, 0.5, 0};
Point(6) = {0.25, 0.5, 0};
Point(7) = {0.5, 0.75, 0};
Point(8) = {0.5, 0.25, 0};
Point(9) = {0.676776695, 0.676776695, 0};
Point(10) = {0.676776695, 0.323223305, 0};
Point(11) = {0.323223305, 0.323223305, 0};
Point(12) = {0.323223305, 0.676776695, 0};
Point(13) = {0.722751631, 0.61, 0};
Point(14) = {0.722751631, 0.39, 0};
Point(15) = {0.277248369, 0.39, 0};
Point(16) = {0.277248369, 0.61, 0};
Point(17) = {0.61, 0.277248369, 0};
Point(18) = {0.39, 0.277248369, 0};
Point(19) = {0.61, 0.722751631, 0};
Point(20) = {0.39, 0.722751631, 0};
Line(1) = {7, 19};
Line(2) = {19, 9};
Line(3) = {9, 13};
Line(4) = {13, 5};
Line(5) = {5, 14};
Line(6) = {14, 10};
Line(7) = {10, 17};
Line(8) = {17, 8};
Line(9) = {8, 18};
Line(10) = {18, 11};
Line(11) = {11, 15};
Line(12) = {15, 6};
Line(13) = {6, 16};
Line(14) = {16, 12};
Line(15) = {12, 20};
Line(16) = {20, 7};
Line(17) = {4, 3};
Line(18) = {3, 2};
Line(19) = {2, 1};
Line(20) = {1, 4};
Mesh.MshFileVersion = 2;
//+
Curve Loop(1) = {13, 14, 15, 16, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {20, 17, 18, 19};
//+
Plane Surface(2) = {1, 2};
//+
Physical Surface(1) = {2};



Mesh.Algorithm = 1;
Mesh.MshFileVersion = 2;
Mesh.RecombinationAlgorithm = 0;
Mesh.CharacteristicLengthMin = 0.01;
Mesh.CharacteristicLengthMax = 0.1;



//+
Physical Surface(2) = {1};
