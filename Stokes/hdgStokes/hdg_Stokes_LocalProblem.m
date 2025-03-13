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
% [u, p, L] = hdg_Stokes_LocalProblem(mesh, refElem, refFace, hdg, uHat, rho, local, nOfComponents)
%
% Inputs:
%
% mesh:             data structure for the mesh
% refElem:          data structure for the reference element
% refFace:          data structure for the reference face
% hdg:              data structure for HDG
% uhat:             hybrid variable on mesh faces
% rho:              constant value of the mean pressure in each element
% local:            local solver (from faces to elements)  
% nOfComponents:    number of  components of the velocity field
%
% Outputs:
%
% u:                primal velocity in mesh elements 
% p:                pressure in mesh elements
% L:                mixed variable in mesh elements
%

function [u, p, L] = hdg_Stokes_LocalProblem(mesh, refElem, refFace, hdg, uHat, rho, local, nOfComponents)

% Count the number of unknowns
nOfDOFP = 0;
for iElem = 1:mesh.nOfElements
    pElem = mesh.pElem(iElem);    
    nOfElementNodes = refElem(pElem).nOfNodes;
    nOfDOFP = nOfDOFP + nOfElementNodes;
end

% Initialisation
p = zeros(nOfDOFP,1);
u = zeros(nOfDOFP*nOfComponents,1);
L = zeros(nOfDOFP*nOfComponents*mesh.nsd,1);

indexPIni = 1;
indexUIni = 1;
indexLIni = 1;
for iElem = 1:mesh.nOfElements
    pElem = mesh.pElem(iElem);
    
    indexGlobalFace = hdgElemToFaceIndex(mesh.indexTf, refFace, hdg.faceInfo(iElem), refElem(pElem).nOfFaces, nOfComponents);
    
    nOfElementNodes = refElem(pElem).nOfNodes;
    nDOFsElemP = nOfElementNodes;
    nDOFsElemU = nOfElementNodes*nOfComponents;
    nDOFsElemL = nOfElementNodes*nOfComponents*mesh.nsd;
    
    indexPEnd = indexPIni + nDOFsElemP - 1;
    indexUEnd = indexUIni + nDOFsElemU - 1;
    indexLEnd = indexLIni + nDOFsElemL - 1;
    
    p(indexPIni:indexPEnd) = local(iElem).Zpl*uHat(indexGlobalFace) + local(iElem).zpr*rho(iElem) + local(iElem).zpf;
    u(indexUIni:indexUEnd) = local(iElem).Zul*uHat(indexGlobalFace) + local(iElem).zur*rho(iElem) + local(iElem).zuf;
    L(indexLIni:indexLEnd) = local(iElem).ZLl*uHat(indexGlobalFace) + local(iElem).zLr*rho(iElem) + local(iElem).zLf;
    
    indexPIni = indexPEnd + 1;
    indexUIni = indexUEnd + 1;
    indexLIni = indexLEnd + 1;
end