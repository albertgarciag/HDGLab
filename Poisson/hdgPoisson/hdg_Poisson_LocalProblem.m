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
% [u,q] = hdg_Poisson_LocalProblem(mesh,refElem,refFace,hdg,uHat,local,nOfComponents)
%
% Inputs:
%
% mesh:             data structure for the mesh
% refElem:          data structure for the reference element
% refFace:          data structure for the reference face
% hdg:              data structure for HDG
% uhat:             hybrid variable on mesh faces
% local:            local solver (from faces to elements)  
% nOfComponents:    number of  components of solution field
%
% Outputs:
%
% u:                primal variable on mesh elements 
% q:                mixed variable on mesh elements
%
% See also: hdg_Poisson_GlobalSystem, hdg_Poisson_ElementalMatrices
%           hdgElemToFaceIndex

function [u,q] = hdg_Poisson_LocalProblem(mesh,refElem,refFace,hdg,uHat,local,nOfComponents)

% Count the number of unknowns
nOfDOFU = 0;
for iElem = 1:mesh.nOfElements
    pElem = mesh.pElem(iElem);    
    nOfElementNodes = refElem(pElem).nOfNodes;
    nOfDOFU = nOfDOFU + nOfElementNodes*nOfComponents;
end

% Initialisation
u = zeros(nOfDOFU,1);
q = zeros(nOfDOFU*mesh.nsd,1);

indexUIni = 1;
indexQIni = 1;
for iElem = 1:mesh.nOfElements
    pElem = mesh.pElem(iElem);
    
    indexGlobalFace = hdgElemToFaceIndex(mesh.indexTf,refFace,hdg.faceInfo(iElem),refElem(pElem).nOfFaces,nOfComponents);
    
    nOfElementNodes = refElem(pElem).nOfNodes;
    nDOFsElemU = nOfElementNodes*nOfComponents;
    nDOFsElemQ = nDOFsElemU*mesh.nsd;
    
    indexUEnd = indexUIni + nDOFsElemU - 1;
    indexQEnd = indexQIni + nDOFsElemQ - 1;
    
    u(indexUIni:indexUEnd) = local(iElem).Zul*uHat(indexGlobalFace) + local(iElem).zuf;
    q(indexQIni:indexQEnd) = local(iElem).Zql*uHat(indexGlobalFace) + local(iElem).zqf;
    
    indexUIni = indexUEnd + 1;
    indexQIni = indexQEnd + 1;
end