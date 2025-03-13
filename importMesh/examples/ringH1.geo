//Inputs
gridsize = 1;
nelem = 5;

Point(1) = {0,0,0,gridsize};
Point(2) = {1,0,0,gridsize};
Point(3) = {-1,0,0,gridsize};
Point(4) = {5,0,0,gridsize};
Point(5) = {-5,0,0,gridsize};

Circle(11) = {2,1,3};
Circle(12) = {3,1,2};
Circle(13) = {4,1,5};
Circle(14) = {5,1,4};

Line(21) = {2,4};
Line(22) = {3,5};

Transfinite Line{11,12,13,14,21,22} = nelem;

Line Loop(31) = {21,13,-22,-11};
Line Loop(32) = {22,14,-21,-12};

Plane Surface(41) = {31};
Plane Surface(42) = {32};

Transfinite Surface{41};
Transfinite Surface{42};

//Define labels for surfaces and volumes
Physical Curve(1) = {11,12,13,14}; // "Dirichlet"
Physical Surface(1) = {41,42}; // "Domain"

//Mesh.ElementOrder =;
// High-order options
// 1: optimisation
// 2: elastic+optimisation
// 3: elastic
Mesh.HighOrderOptimize = 3; 

Mesh 2;

Save "ringH1unif.msh";
