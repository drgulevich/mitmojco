/*********************************************************************
 *
 *  Generation of mesh for FFO with tapered ends 
 *
 *********************************************************************/

lceff = 0.11; // effective element size after refinement
lc = 2*lceff;

lambdaJ = 5.5;

// MJTL
L0 = 400./lambdaJ;
W0 = 16./lambdaJ;
S0 = 40.0/lambdaJ;
P0 = 2.0/lambdaJ;

//AJTL
L1 = 228./lambdaJ;
W1 = 4./lambdaJ;
S1 = 20.0/lambdaJ;
P1 = 2.0/lambdaJ;

Point(1) = {-0.5*L0, 0.5*(W0-P0), 0, lc};
Point(2) = {-0.5*L0+S0, 0,  0, lc};
Point(3) = {-0.5*W1, 0,  0, lc};
Point(4) = {0.5*W1, 0,  0, lc};
Point(5) = {0.5*L0-S0, 0, 0, lc};
Point(6) = {0.5*L0, 0.5*(W0-P0), 0, lc};
Point(7) = {0.5*L0, 0.5*(W0+P0), 0, lc};
Point(8) = {0.5*L0-S0, W0, 0, lc};
Point(9) = {0.5*W1, W0, 0, lc};
Point(10) = {0.5*W1, L1-S1, 0, lc};
Point(11) = {0.5*P1, L1, 0, lc};
Point(12) = {-0.5*P1, L1, 0, lc};
Point(13) = {-0.5*W1, L1-S1, 0, lc};
Point(14) = {-0.5*W1, W0, 0, lc};
Point(15) = {-0.5*L0+S0, W0,  0, lc};
Point(16) = {-0.5*L0, 0.5*(W0+P0), 0, lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,9};
Line(9) = {9,10};
Line(10) = {10,11};
Line(11) = {11,12};
Line(12) = {12,13};
Line(13) = {13,14};
Line(14) = {14,15};
Line(15) = {15,16};
Line(16) = {16,1};

Line(17) = {2,15}; // vertical internal line, left
Line(18) = {3,14}; // vertical internal line, middle left
Line(19) = {4,9}; // vertical internal line, middle right
Line(20) = {5,8}; // vertical internal line, right
Line(21) = {14,9}; // internal line, middle
Line(22) = {13,10}; // internal line, top

Line Loop(23) = {1,17,15,16};
Line Loop(24) = {2,18,14,-17};
Line Loop(25) = {3,19,-21,-18};
Line Loop(26) = {4,20,8,-19};
Line Loop(27) = {5,6,7,-20};
Line Loop(28) = {21,9,-22,13};
Line Loop(29) = {22,10,11,12};

Plane Surface(30) = {23};
Plane Surface(31) = {24};
Plane Surface(32) = {25};
Plane Surface(33) = {26};
Plane Surface(34) = {27};
Plane Surface(35) = {28};
Plane Surface(36) = {29};

Transfinite Line{2,-14,-4,8} = Ceil( (0.5*L0-S0-0.5*W1)/lc ) Using Progression 1;
Transfinite Line{-1,15,5,-7} = Ceil(S0/lc) Using Progression 1;
Transfinite Line{17,18,19,20} = Ceil( W0/lc ) Using Progression 1;
Transfinite Line{-16,6} = 3 Using Progression 1;
Transfinite Line{3,21,22} = Ceil(W1/lc) Using Progression 1; // central rectangle
Transfinite Line{9,13} = Ceil( (L1-S1-W0)/lc ) Using Progression 1;
Transfinite Line{10,-12} = Ceil(S1/lc) Using Progression 1;
Transfinite Line{11} = 3 Using Progression 1;
 
Transfinite Surface{30,31,32,33,34,35,36};
//Recombine Surface{14,15,16};

// some parameters for the meshing:
Mesh.Algorithm = 8;
//Recombine Surface{14,15,16};
Mesh.RecombineAll = 1;
//Mesh.CharacteristicLengthFactor = lc;
Mesh.SubdivisionAlgorithm = 1;
//Mesh.Smoothing = 1;
Mesh.Optimize

//Show "*";

Coherence;

Physical Line(0) = {9,13};
Physical Line(1) = {16};
Physical Line(2) = {6}; // MJTL radiation end
Physical Line(3) = {2,3,4}; // bottom of MJTL
Physical Line(4) = {1};
Physical Line(5) = {5};
Physical Line(6) = {10,12};
Physical Line(7) = {11}; // AJTL radiation end
Physical Line(8) = {8,14}; // top of MJTL
Physical Line(9) = {15};
Physical Line(10) = {7};

Physical Surface(0) = {30,31,32,33,34,35,36};
