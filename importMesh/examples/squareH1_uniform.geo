//Inputs
gridsize = 0.5;

Point(1) = {0,0,0,gridsize};
Point(2) = {1,0,0,gridsize};
Point(3) = {1,1,0,gridsize};
Point(4) = {0,1,0,gridsize};

Line(5) = {1,2};
Line(6) = {2,3};
Line(7) = {3,4};
Line(8) = {4,1};

Line Loop(9) = {5,6,7,8};

Plane Surface(10) = 9;

//Define labels for surfaces and volumes
Physical Curve(1) = {6,7,8}; // "Dirichlet"
Physical Curve(2) = {5}; // "Neumann"
Physical Surface(1) = {10}; // "Domain"

//Mesh.ElementOrder =;
// High-order options
// 1: optimisation
// 2: elastic+optimisation
// 3: elastic
Mesh.HighOrderOptimize = 3; 

Mesh 2;

Save "squareH1unif.msh";