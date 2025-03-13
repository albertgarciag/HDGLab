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
%
% mesh = buildMeshStruct(X, T, matElem, intFaces, extFaces, refElem, optionNodes)
%
% Inputs:
%
% X:              vector of the nodes
% T:              structure containing the element connectivity
% matElem:        matrix storing the ID of the region to which each 
%                 element belongs
% intFaces:       structure containing the internal faces (seen from left 
%                 and right) and their connectivity
% extFaces:       structure containing the external faces, the ID of the 
%                 contour and the imposed boundary condition
% refElem:        reference element 
% optionNodes:    flag for Fekete (=1) or uniform (=0) nodal distributions
%
% Outputs:
%
% mesh:           mesh structure
%

function mesh = buildMeshStruct(X, T, matElem, intFaces, extFaces, refElem, optionNodes)

%% Store mesh information
% Number of spatial dimensions
mesh.nsd = refElem.nsd;

% Number of elements
mesh.nOfElements = size(T,1);

% Number of nodes
mesh.nOfNodes = 0;

% ID of the region to which each element belongs (useful for multi-material problems)
mesh.matElem = matElem;

% Internal faces structure
mesh.intFaces = intFaces;

% Number of internal faces
mesh.nOfIntFaces = size(intFaces,1);

% Flag identifying nodal distribution
% 0: Uniform nodes
% 1: Fekete nodes
mesh.optionNodes = optionNodes;

% External faces structure
mesh.extFaces = extFaces;

% Number of external faces
mesh.nOfExtFaces = size(extFaces,1);

% Degree of approximation in each element (useful for non-uniform polynomial approximations)
mesh.pElem = refElem.p*ones(1, mesh.nOfElements);

%% Initialise mesh structures
% Connectivity matrix
mesh.indexT = zeros(mesh.nOfElements, 2);

% Nodes structure
mesh.X = zeros(mesh.nOfNodes, mesh.nsd);

%% Construct the connectivity matrix
indexIni = 1;
for iElem = 1:mesh.nOfElements
    indexEnd = indexIni + refElem.nOfNodes - 1;
    mesh.indexT(iElem, 1) = indexIni;
    mesh.indexT(iElem, 2) = indexEnd;
    mesh.X(indexIni:indexEnd, :) = X(T(iElem,:), :);
    indexIni = indexEnd + 1;
end
mesh.nOfNodes = size(mesh.X,1);