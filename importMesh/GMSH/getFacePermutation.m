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
% iPerm = getFacePermutation(vectorNodes, targetFirstNode, nsd, nOfVertices)
%
% Inputs:
%
% iPerm:                permutation of the nodes due to orientation
%
% Outputs:
%
% vectorNodes:          nodes of the face under analysis
% targetFirstNode:      first node of the face whose neighbour is sought
% nsd:                  number of spatial dimensions
% nOfFaceVertices:      number of face vertices
%

function iPerm = getFacePermutation(vectorNodes, targetFirstNode, nsd, nOfVertices)

if nsd==2
    vecNodes = [vectorNodes(1) vectorNodes(length(vectorNodes))];
else
    vecNodes = vectorNodes(1:nOfVertices);
end

iPerm = find(vecNodes==targetFirstNode);