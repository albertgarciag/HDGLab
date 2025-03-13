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
% extFaces = getBoundaryFaces(T, faces, refElem, refFace)
%
% Inputs:
%
% T:              structure containing the element connectivity
% faces:          structure containing the information on the boundary faces
% refElem:        reference element 
% refFace:        reference face
%
% Outputs:
%
% extFaces:       structure containing the external faces, the ID of the 
%                 contour and the imposed boundary condition
%

function extFaces = getBoundaryFaces(T, faces, refElem, refFace)

nOfExtFaces = size(faces, 1);
nOfElemFaces = refElem.nOfFaces;
nOfFaceNodes = refFace.nOfNodes;

%% Initialise the structure for boundary faces
extFaces = zeros(nOfExtFaces,4);

%% Loop over external faces extracted by the mesh generator
for kFace = 1:nOfExtFaces
    iElem = faces(kFace,1);
    iNodes = sort(faces(kFace,2:nOfFaceNodes+1));
    iNurbs = faces(kFace,nOfFaceNodes+2);
    iBCfunc = faces(kFace,nOfFaceNodes+3);
    
    % Find the number of face by matching nodes
    for iFace = 1:nOfElemFaces
        nodesReferenceFace = refElem.face(iFace).nodes;
        nodesFace = sort(T(iElem, nodesReferenceFace));
        
        if all(ismember(nodesFace, iNodes))
            extFaces(kFace,:) = [iElem, iFace, iBCfunc, iNurbs];
            break;
        end
    end
end