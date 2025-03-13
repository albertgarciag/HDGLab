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
% [mesh, hdg] = hdgPreprocess(mesh, hdg, refElem, refFace, ctt)
%
% Inputs:
%
% mesh:     data structure for the mesh
% refElem:  data structure for the reference element
% refFace:  data structure for the reference face
% hdg:      data structure for HDG
% ctt:      contains flags for boundary conditions and number of components
%           of solution field
%
% Outputs:
%
% mesh:     data structure for the mesh
%           adds field indexTf
% hdg:      data structure for HDG
%           adds fields faceInfo, nDOFglobal, vDOFtoSolve
%
% See also: getRefData, hdg_Poisson_Constants


function [mesh, hdg] = hdgPreprocess(mesh,refElem,refFace,hdg,ctt)

% Initialisation
hdg.faceInfo(mesh.nOfElements).local2global = [];
hdg.faceInfo(mesh.nOfElements).localNumFlux = [];
hdg.faceInfo(mesh.nOfElements).permutations = [];
hdg.faceInfo(mesh.nOfElements).pHat = [];
hdg.nDOFglobal = 0;

for iElem = 1:mesh.nOfElements
    pElem = mesh.pElem(iElem);
    nOfElementFaces = refElem(pElem).nOfFaces;
    hdg.faceInfo(iElem).local2global = zeros(1,nOfElementFaces);
    hdg.faceInfo(iElem).localNumFlux = zeros(1,nOfElementFaces);
    hdg.faceInfo(iElem).permutations = zeros(1,nOfElementFaces);
    hdg.faceInfo(iElem).pHat = zeros(1,nOfElementFaces);
end

% Internal faces ----------------------------------------------------------
indexFaceIni = 1;
for iFaceGlobal = 1:mesh.nOfIntFaces
    iElemL = mesh.intFaces(iFaceGlobal,1);
    iFaceL = mesh.intFaces(iFaceGlobal,2);
    iElemR = mesh.intFaces(iFaceGlobal,3);
    iFaceR = mesh.intFaces(iFaceGlobal,4);
    iPerm  = mesh.intFaces(iFaceGlobal,5);
    
    % From local to global numbers
    hdg.faceInfo(iElemL).local2global(iFaceL) = iFaceGlobal;
    hdg.faceInfo(iElemR).local2global(iFaceR) = iFaceGlobal;
    
    % Flag for the numerical flux
    hdg.faceInfo(iElemL).localNumFlux(iFaceL) = ctt.iBC_Interior;
    hdg.faceInfo(iElemR).localNumFlux(iFaceR) = ctt.iBC_Interior;
    
    % Permutation
    hdg.faceInfo(iElemR).permutations(iFaceR) = iPerm;
    
    % Number of degrees of freedom
    pL = mesh.pElem(iElemL);
    nOfFaceNodesL = refFace(pL,pL).nOfNodes;
    pR = mesh.pElem(iElemR);
    nOfFaceNodesR = refFace(pR,pR).nOfNodes;
    nOfFaceNodesMax = max([nOfFaceNodesL,nOfFaceNodesR]);
    nOfFaceUnknowns = nOfFaceNodesMax*ctt.nOfComponents;
    hdg.nDOFglobal = hdg.nDOFglobal + nOfFaceUnknowns;
    
    % p of uHat (max of each face shared)
    maxP = max([pL,pR]);
    hdg.faceInfo(iElemL).pHat(iFaceL) = maxP;
    hdg.faceInfo(iElemR).pHat(iFaceR) = maxP;
    
    % Global indexing of face unknowns
    indexFaceEnd = indexFaceIni + nOfFaceUnknowns - 1;
    mesh.indexTf(iFaceGlobal,:) = [indexFaceIni, indexFaceEnd];
    indexFaceIni = indexFaceEnd + 1;
end

% External faces ----------------------------------------------------------
for iFaceGlobalExt = 1:mesh.nOfExtFaces
    iElemL = mesh.extFaces(iFaceGlobalExt,1);
    iFaceL = mesh.extFaces(iFaceGlobalExt,2);
    iBC = mesh.extFaces(iFaceGlobalExt,3);
    
    % From local to global numbers
    iFaceGlobal = iFaceGlobalExt + mesh.nOfIntFaces;
    hdg.faceInfo(iElemL).local2global(iFaceL) = iFaceGlobal;
    
    % Flag for the numerical flux
    hdg.faceInfo(iElemL).localNumFlux(iFaceL) = iBC;
    
    % Number of degrees of freedom
    pL = mesh.pElem(iElemL);
    nOfFaceNodesL = refFace(pL,pL).nOfNodes;
    nOfFaceUnknowns = nOfFaceNodesL*ctt.nOfComponents;
    hdg.nDOFglobal = hdg.nDOFglobal + nOfFaceUnknowns;
    
    % p of uHat
    hdg.faceInfo(iElemL).pHat(iFaceL) = pL;
    
    % Global indexing of face unknowns
    indexFaceEnd = indexFaceIni + nOfFaceUnknowns - 1;
    mesh.indexTf(iFaceGlobal,:) = [indexFaceIni, indexFaceEnd];
    indexFaceIni = indexFaceEnd + 1;
end

% Degrees of freedom to solve for -----------------------------------------
if strcmp(hdg.problem,'Poisson') | strcmp(hdg.problem,'ConvDiff')
    isDOFtoSolve = ones(1,hdg.nDOFglobal);
elseif strcmp(hdg.problem,'Stokes')
    hdg.columnsGlobalSystem = hdg.nDOFglobal + mesh.nOfElements;
    % If incompressible flow with purely Dirichlet BC, add a line for the
    % restriction to solve the underdetermination of pressure
    %    if all(mesh.extFaces(:,3) == ctt.iBC_Dirichlet)
    if not(any(mesh.extFaces(:,3) == ctt.iBC_Neumann))
        hdg.pureDirichlet = 1;
        hdg.rowsGlobalSystem = hdg.columnsGlobalSystem + 1;
    else
        hdg.pureDirichlet = 0;
        hdg.rowsGlobalSystem = hdg.columnsGlobalSystem;
    end
    isDOFtoSolve = ones(1,hdg.rowsGlobalSystem);
end

for iFaceGlobalExt = 1:mesh.nOfExtFaces
    iFaceGlobal = iFaceGlobalExt + mesh.nOfIntFaces;
    indexGlobalFace = mesh.indexTf(iFaceGlobal,1):mesh.indexTf(iFaceGlobal,2);
    iBC = mesh.extFaces(iFaceGlobalExt,3);
    if iBC == ctt.iBC_Dirichlet
        isDOFtoSolve(indexGlobalFace) = 0;
    end
end

hdg.vDOFtoSolve = find(isDOFtoSolve==1);
