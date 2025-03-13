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
% [jElem, jFace, jPerm] = findNeighbourElement(T, iElem, faceNodes, listPotentialNeigh, elemFaces, nsd, nOfElemFaces, nOfFaceVertices)
%
% Inputs:
%
% jElem:               ID of the neighbouring element
% jFace:               ID of the shared face
% jPerm:               permutation of the nodes due to orientation
%
% Outputs:
%
% T:                    structure containing the element connectivity
% iElem:                ID of the element under analysis
% facesNodes:           IDs of the nodes of the face under analysis
% listPotentialNeigh:   list of potential neighbours of the current element
% elemFaces:            structure containing the information on the faces
%                       of the reference element
% nsd:                  number of spatial dimensions
% nOfElemFaces:         number of element faces
% nOfFaceVertices:      number of face vertices
%
% See also: getFacePermutation
%

function [jElem, jFace, jPerm] = findNeighbourElement(T, iElem, faceNodes, listPotentialNeigh, elemFaces, nsd, nOfElemFaces, nOfFaceVertices)

sortFaceNodes = sort(faceNodes);

% Initialisation
jElem = 0;
jFace = 0;
jPerm = 0;

% Loop over potential neighbours
nOfPotentialNeigh = length(listPotentialNeigh(1, :));
for iPotentialNeigh = 1:nOfPotentialNeigh
    jIndexElem = listPotentialNeigh(1, iPotentialNeigh);
    if ~(jIndexElem == iElem)
        % Loop over faces of potential neighbour
        for jFace = 1:nOfElemFaces
            jFaceNodes = T(jIndexElem, elemFaces(jFace).nodes);
            % Compare set of nodes after sorting
            sortJFaceNode = sort(jFaceNodes);
            if sum(abs(sortFaceNodes-sortJFaceNode)) == 0
                jElem = jIndexElem;
                jPerm = getFacePermutation(jFaceNodes, faceNodes(1), nsd, nOfFaceVertices);
                return;
            end
        end
    end
end