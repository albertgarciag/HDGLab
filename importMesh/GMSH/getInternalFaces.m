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
% intFaces = getInternalFaces(T, refElem, refFace)
%
% Inputs:
%
% T:              structure containing the element connectivity
% refElem:        reference element 
% refFace:        reference face
%
% Outputs:
%
% intFaces:       structure containing the internal faces (seen from left 
%                 and right) and their connectivity
%
% See also: findNeighbourElement
%

function intFaces = getInternalFaces(T, refElem, refFace)

nOfElements = size(T, 1);
nOfNodes = max(max(T));
nOfElemNodes = refElem.nOfNodes;
nOfElemFaces = refElem.nOfFaces;

%% Count maximum number of elements sharing a node
vCount  = zeros(1, nOfNodes);
for kElem = 1:nOfElements
    vCount(T(kElem, :)) = vCount(T(kElem, :)) + 1;
end
nOfMaxElemsOneNode = max(vCount);

%% Store the list of elements containing each node
nodeToElem = zeros(nOfNodes, nOfMaxElemsOneNode);
indexNode = zeros(nOfNodes, 1);
for kElem = 1:nOfElements
    Te = T(kElem, :);
    for iLocalNode = 1:nOfElemNodes
        iGlobalNode = Te(iLocalNode);
        indexNode(iGlobalNode) = indexNode(iGlobalNode) + 1;
        nodeToElem(iGlobalNode, indexNode(iGlobalNode)) = kElem;
    end
end

%% Identify internal faces
nOfMaxFaces = nOfElements*nOfElemFaces;
intFaces = zeros(nOfMaxFaces, 5);
E = zeros(nOfElements, nOfElemFaces);

indexInt = 0;
for iElem = 1:nOfElements
    for iFace = 1:nOfElemFaces
        if E(iElem,iFace) == 0
            E(iElem,iFace) = 1;
            faceNodes = T(iElem, refElem.face(iFace).nodes);
            iNodeFirst = faceNodes(1);
            listPotentialNeigh = nodeToElem(iNodeFirst, nodeToElem(iNodeFirst, :)>0);
            [jElem, jFace, jPerm] = findNeighbourElement(T, iElem, faceNodes, listPotentialNeigh, refElem.face, refElem.nsd, nOfElemFaces, refFace.nOfVertices);
            if jElem > 0
                indexInt = indexInt + 1;
                E(jElem, jFace) = 1;
                intFaces(indexInt, :) = [iElem, iFace, jElem, jFace, jPerm];
            end
        end
    end
end
intFaces = intFaces(1:indexInt, :);