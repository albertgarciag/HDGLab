To generate a mesh of degree 1 using GMSH, run the following command via terminal:

gmsh file.geo


To generate a mesh of degree P>1 using GMSH, run the following command via terminal:

gmsh -degree P file.geo

P: Degree 2,..,5

High-order option to be specified in the .geo file via Mesh.HighOrderOptimize = M; 
// M=1: optimisation
// M=2: elastic+optimisation
// M=3: elastic


All information on GMSH is available at http://gmsh.info
