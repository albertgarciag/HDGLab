%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2020 Ruben Sevilla, Matteo Giacomini and Antonio Huerta
% Contact:      r.sevilla@swansea.ac.uk
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
% [indexGlobalFace,nOfFaceDOF] = hdgElemToFaceIndex(indexTf,refFace,faceInfo,nOfFaces,nOfComponents)
%
% Inputs:
%
% indexTf:       global indices for the face nodal values
% refFace:       data structure for the reference face
% faceInfo:      data structure for the faces of an element
% nOfFaces:      number of element faces
% nOfComponents: number of components of solution field
%
% Outputs:
%
% indexGlobalFace:  global indices for all element face nodal values
% nOfFaceDOF:       number of degrees of freedom for all element faces
%
% See also: hdg_Poisson_GlobalSystem, hdg_Poisson_LocalProblem

function [indexGlobalFace,nOfFaceDOF] = hdgElemToFaceIndex(indexTf,refFace,faceInfo,nOfFaces,nOfComponents)

nOfNodesAllFaces = 0;
for iFace = 1:nOfFaces
    pHat = faceInfo.pHat(iFace);
    nOfNodesAllFaces = nOfNodesAllFaces + refFace(pHat,pHat).nOfNodes;
end

nOfFaceDOF = nOfNodesAllFaces*nOfComponents;
indexGlobalFace = zeros(1, nOfFaceDOF);
indexFaceIni = 1;
for iFace = 1:nOfFaces    
    pHat = faceInfo.pHat(iFace);
    nOfFaceNodes = nOfComponents*refFace(pHat,pHat).nOfNodes;
    indexFaceEnd = indexFaceIni + nOfFaceNodes - 1;
    indexFaceV = indexFaceIni:indexFaceEnd;
    
    iFaceGlobal = faceInfo.local2global(iFace);
    indexGlobalFace(indexFaceV) = indexTf(iFaceGlobal,1):indexTf(iFaceGlobal,2);
    
    indexFaceIni = indexFaceEnd + 1;
end