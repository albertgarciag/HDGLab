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
% [X, T, faces, matElem] = scanMeshFileMSH(meshFile, pDegree, nsd)
%
% Inputs:
%
% meshFile:       name of the mesh file
% pDegree:        polynomial degree of approximation
% nsd:            number of spatial dimensions
%
% Outputs:
%
% X:              vector of the nodes
% T:              structure containing the element connectivity
% faces:          structure containing the information on the boundary faces
% matElem:        matrix storing the ID of the region to which each 
%                 element belongs
%
% See also: getElemFaceDataMSH, nodeRenumberingMSH
%

function [X, T, faces, matElem, nOfElemNodes] = scanMeshFileMSH(meshFile, pDegree, nsd)

%% Retrieve the geometrical information on mesh elements and faces
[elemID, faceID, lowerEntities, nOfElemNodes, nOfFaceNodes] = getElemFaceDataMSH(pDegree, nsd);

%% Open the mesh file
fileName = sprintf('%s.msh',meshFile);
fid = fopen(fileName, 'r');
for iLine=1:4
    neglectLine = fgets(fid);
end

%% Get nodes
formatPoint = '%d %f %f %f';

nOfNodes = fscanf(fid, '%d', 1);
vecNodes = ( fscanf(fid, formatPoint, [4,nOfNodes]) )';
X = vecNodes(:,2:nsd+1);

%% Read geometrical entities (from points to volumes)
for iLine=1:3
    neglectLine = fgets(fid);
end
nOfEntities = fscanf(fid, '%d', 1);
neglectLine = fgets(fid);

% Discard 1-node points
tline = fgetl(fid);
info = sscanf(tline,'%d')';
typeEntity = info(2)';

nOfDiscardedLines = 0;
while any(typeEntity==lowerEntities)
    nOfDiscardedLines = nOfDiscardedLines + 1;
    tline = fgetl(fid);
    info = sscanf(tline,'%d')';
    typeEntity = info(2);
    if ~any(typeEntity==lowerEntities)
        break;
    end
end
nOfAttributes = info(3);

%% Get faces
faceDetails = zeros(nOfEntities-nOfDiscardedLines,nOfAttributes);
faceVertices = zeros(nOfEntities-nOfDiscardedLines,nOfFaceNodes);

iFace = 1;
nOfFaceDetails = nOfAttributes + nOfFaceNodes;
faceDetails(iFace,:) = info(:,4:4+nOfAttributes-1);
faceVertices(iFace,:) = info(:,4+nOfAttributes:4+nOfFaceDetails-1);

while any(typeEntity==faceID)
    iFace = iFace+1;
    tline = fgetl(fid);
    info = sscanf(tline,'%d')';
    typeEntity = info(2);
    if ~any(typeEntity==faceID)
        break;
    end
    
    nOfFaceDetails = nOfAttributes + nOfFaceNodes;
    faceDetails(iFace,:) = info(:,4:4+nOfAttributes-1);
    faceVertices(iFace,:) = info(:,4+nOfAttributes:4+nOfFaceDetails-1);
end
nOfFaces = iFace - 1;
faceVertices = faceVertices(1:nOfFaces,:);
faceDetails = faceDetails(1:nOfFaces,:);

%% Get elements and region IDs
nOfElements = nOfEntities - nOfDiscardedLines - nOfFaces;
T = zeros(nOfElements,nOfElemNodes);
matElem = zeros(1,nOfElements);

iElem = 1;
nOfElemDetails = nOfAttributes + nOfElemNodes;
T(iElem,:) = info(:,4+nOfAttributes:4+nOfElemDetails-1);
matElem(iElem) = info(:,4);

while any(typeEntity==elemID) && iElem<nOfElements
    iElem = iElem+1;
    tline = fgetl(fid);
    info = sscanf(tline,'%d')';
    typeEntity = info(2);
    if ~any(typeEntity==elemID)
        break;
    end
    
    nOfElemDetails = nOfAttributes + nOfElemNodes;    
    T(iElem,:) = info(:,4+nOfAttributes:4+nOfElemDetails-1);
    matElem(iElem) = info(:,4);
end

%% Node permutation from GMSH to Matlab ordering
permNodes = nodeRenumberingMSH(pDegree, nsd);
T(:, 1:nOfElemNodes) = T(:, permNodes);

%% Gather information on boundary faces
faces = zeros(nOfFaces, nOfAttributes+nOfFaceNodes+1);
for iFace = 1:nOfFaces
    faceNodes = sort(faceVertices(iFace, :));
    
    for iElem = 1:nOfElements
        elemNodes = sort(T(iElem,1:nOfElemNodes));
        if ismember(faceNodes,elemNodes)
           break; 
        end
    end
    
    % Store
    faces(iFace,1) = iElem;
    faces(iFace,2:nOfFaceNodes+1) = faceVertices(iFace, :);
    faces(iFace,nOfFaceNodes+2) = faceDetails(iFace, 2); %NURBS ID
    faces(iFace,nOfFaceNodes+3) = faceDetails(iFace, 1); %BC type

end

fclose(fid);