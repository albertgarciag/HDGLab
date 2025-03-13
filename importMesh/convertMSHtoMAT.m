%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Convert mesh from .msh format to .mat format
%  Supports meshes of triangles (TRI) in 2D and tetrahedra (TET) in 3D
%  up to polynomial degree 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2020 Ruben Sevilla, Matteo Giacomini and Antonio Huerta    
% Contact:      matteo.giacomini@upc.edu                                   
% ZCCE (SU)   - https://www.swansea.ac.uk/engineering/zcce/
% LaCaN (UPC) - https://www.lacan.upc.edu/                                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars, clc
close all
setpathMesh

optionNodes = 0; % GMSH utilises uniform nodal distributions

%% Option setup
fileName = 'squareH1P1unif';
nsd = 2;          % Number of spatial dimensions
pDegree = 1;      % Polynomial degree
isPlotMesh = 1;   % Boolean variable to plot the imported mesh
outputPath = 'meshFiles';

%% Construct the reference element
[refElem, refFace] = getRefData(pDegree, nsd, optionNodes);
refElem = refElem(1, pDegree);
refFace = refFace(pDegree, pDegree);

%% Store the coordinates of the nodes with duplication
[X, T, faces, matElem] = scanMeshFileMSH(fileName, pDegree, nsd);

%% Extract internal and external faces
intFaces = getInternalFaces(T, refElem, refFace);
extFaces = getBoundaryFaces(T, faces, refElem, refFace);

%% Build the mat structure of the mesh
mesh = buildMeshStruct(X, T, matElem, intFaces, extFaces, refElem, optionNodes);

%% Save mesh structure in .mat file
meshFile = sprintf('%s/%s', outputPath, fileName);
save(meshFile, 'mesh');

%% Visualise the imported mesh
if isPlotMesh
    visualiseMesh
end