//Inputs
boxdim = 1;
gridsize = 0.5;

//Create 2D square mesh
Point(1) = {0,0,0,gridsize};
Point(2) = {boxdim,0,0,gridsize};
Point(3) = {boxdim,boxdim,0,gridsize};
Point(4) = {0,boxdim,0,gridsize};

Line(5) = {1,2};
Line(6) = {2,3};
Line(7) = {3,4};
Line(8) = {4,1};

Line Loop(9) = {5,6,7,8};
Plane Surface(10) = 9;

Transfinite Line{5,6,7,8} = boxdim/gridsize+boxdim;
Transfinite Surface{10};

//Now make 3D by extrusion.
newEntities[] = 
Extrude { 0,0,boxdim }
{
	Surface{10};
	Layers{boxdim/gridsize};
};

//Define labels for surfaces and volumes
// front: newEntities[0]
// back: 10
// left: newEntities[5]
// right: newEntities[3]
// top: newEntities[4]
// bottom: newEntities[2]
Physical Surface(1) = {10,newEntities[0],newEntities[5],newEntities[3],newEntities[4]}; // "Dirichlet"
Physical Surface(2) = {newEntities[2]}; // "Neumann"
Physical Volume(1) = {newEntities[1]}; // "Domain"

//Mesh.ElementOrder =;
// High-order options
// 1: optimisation
// 2: elastic+optimisation
// 3: elastic
Mesh.HighOrderOptimize = 3; 

Mesh 3;

Save "cubeH1unif.msh";
