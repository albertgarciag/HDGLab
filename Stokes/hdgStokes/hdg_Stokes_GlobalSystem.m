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
% [uhat, rho, local] = hdg_Stokes_GlobalSystem(mesh, refElem, refFace, hdg, ctt, problemParams)
%
% Inputs:
%
% mesh:             data structure for the mesh
% refElem:          data structure for the reference element
% refFace:          data structure for the reference face
% hdg:              data structure for HDG
% ctt:              contains flags for boundary conditions and number of 
%                   components of solution field
% problemParams:    problem specific parameters
%
% Outputs:
%
% uhat:             hybrid variable on mesh faces
% rho:              constant value of the mean pressure in each element
% local:            local solver (from faces to elements)  
%
% See also: hdg_Stokes_ElementalMatrices, hdgElemToFaceIndex
%

function [uhat, rho, local] = hdg_Stokes_GlobalSystem(mesh, refElem, refFace, hdg, ctt,problemParams)

% Initialisation
mat.i = zeros(1, hdg.rowsGlobalSystem);
mat.j = zeros(1, hdg.columnsGlobalSystem);
mat.Kij = zeros(1, hdg.nDOFglobal);
f = zeros(hdg.rowsGlobalSystem, 1);

local(mesh.nOfElements).Zul = [];
local(mesh.nOfElements).Zpl = [];
local(mesh.nOfElements).ZLl = [];
local(mesh.nOfElements).zuf = [];
local(mesh.nOfElements).zpf = [];
local(mesh.nOfElements).zLf = [];
local(mesh.nOfElements).zur = [];
local(mesh.nOfElements).zpr = [];
local(mesh.nOfElements).zLr = [];

% Integral of the pressure over the domain (used only if purely Dirichlet BC)
intPressure = 0;

indexIni = 0;
indexIni2 = 0;
for iElem = 1:mesh.nOfElements
    Te = getElemConnectivity(mesh, iElem);
    Xe = mesh.X(Te,:);
    pElem = mesh.pElem(iElem);
    matElem = mesh.matElem(iElem);
    faceInfo = hdg.faceInfo(iElem);
    
    isAnySlipFace = any(faceInfo.localNumFlux == ctt.iBC_Slip) || any(faceInfo.localNumFlux == ctt.iBC_Axi);
    
    % Elemental matrices (from faces to elements)  
    if problemParams.axi==1
        [Zule, Zple, ZLle, zufe, zpfe, zLfe, zure, zpre, zLre, Ke, fe, intPressure] = ...
            hdg_Stokes_ElementalMatricesAxi(refElem, refFace, Xe, pElem, matElem, hdg.tau, faceInfo, ctt, problemParams, hdg.pureDirichlet, isAnySlipFace, intPressure);
    else
        [Zule, Zple, ZLle, zufe, zpfe, zLfe, zure, zpre, zLre, Ke, fe, intPressure] = ...
            hdg_Stokes_ElementalMatrices(refElem, refFace, Xe, pElem, matElem, hdg.tau, faceInfo, ctt, problemParams, hdg.pureDirichlet, isAnySlipFace, intPressure);
    end

    % Store the local solver
    local(iElem).Zul = Zule;
    local(iElem).Zpl = Zple;
    local(iElem).ZLl = ZLle;
    local(iElem).zuf = zufe;
    local(iElem).zpf = zpfe;
    local(iElem).zLf = zLfe;
    local(iElem).zur = zure;
    local(iElem).zpr = zpre;
    local(iElem).zLr = zLre;
    
    % Assembly
    [indexGlobalFace, ~] = hdgElemToFaceIndex(mesh.indexTf, refFace, hdg.faceInfo(iElem), refElem(pElem).nOfFaces, ctt.nOfComponents);
    if hdg.pureDirichlet
        iBlock = [indexGlobalFace, hdg.nDOFglobal+iElem, hdg.nDOFglobal+mesh.nOfElements+1];
        jBlock = [indexGlobalFace, hdg.nDOFglobal+iElem];
    else
        iBlock = [indexGlobalFace, hdg.nDOFglobal+iElem];
        jBlock = iBlock;
    end

    nI =numel(iBlock);
    nJ =numel(jBlock);
    currentIndex = indexIni + (1:nJ);
    currentIndex2 = indexIni2 + (1:nI*nJ);
    for i=1:nI
        mat.i(currentIndex) = iBlock(i);
        mat.j(currentIndex) = jBlock;
        currentIndex = currentIndex + nJ;
    end
    mat.Kij(currentIndex2) = reshape(Ke', 1, nI*nJ);
    
    % RHS
    f(iBlock) = f(iBlock) + fe;
    
    % Update indices
    indexIni = indexIni + nI*nJ;
    indexIni2 = indexIni2 + nI*nJ;
end

K = sparse(mat.i, mat.j, mat.Kij);
% If purely Dirichlet BC, add pressure constraint
if hdg.pureDirichlet
    f(hdg.rowsGlobalSystem) = -f(hdg.rowsGlobalSystem)+intPressure;
    K = [K, [K(end,:)'; 0]];
end

% Renumbering of the matrix
reOrd = symrcm(K(hdg.vDOFtoSolve, hdg.vDOFtoSolve));
% Solving the linear system
sol = zeros(hdg.columnsGlobalSystem, 1);
sol(hdg.vDOFtoSolve(reOrd)) = K(hdg.vDOFtoSolve(reOrd), hdg.vDOFtoSolve(reOrd))\f(hdg.vDOFtoSolve(reOrd)); 
% Solution
uhat = sol(1:hdg.nDOFglobal);
rho = sol(hdg.nDOFglobal+1:end);