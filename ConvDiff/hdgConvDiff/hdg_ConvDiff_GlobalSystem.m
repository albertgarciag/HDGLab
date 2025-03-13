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
% [uhat, local] = hdg_Poisson_GlobalSystem(mesh,refElem,refFace,hdg,ctt,problemParams)
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
% local:            local solver (from faces to elements)  
%
% See also: hdg_Poisson_ElementalMatrices, hdg_Poisson_LocalProblem
%           hdgElemToFaceIndex

function [uhat, local] = hdg_ConvDiff_GlobalSystem(mesh,refElem,refFace,hdg,ctt,problemParams)

% Initialisation
mat.i = zeros(1,hdg.nDOFglobal);
mat.j = zeros(1,hdg.nDOFglobal);
mat.Kij = zeros(1,hdg.nDOFglobal);
f = zeros(hdg.nDOFglobal,1);

local(mesh.nOfElements).Zql = [];
local(mesh.nOfElements).Zul = [];
local(mesh.nOfElements).zqf = [];
local(mesh.nOfElements).zuf = [];

indexIni = 0;
for iElem = 1:mesh.nOfElements
    Te = getElemConnectivity(mesh, iElem);
    Xe = mesh.X(Te,:);
    pElem = mesh.pElem(iElem);
    matElem = mesh.matElem(iElem);
    faceInfo = hdg.faceInfo(iElem);
    
    % Elemental matrices (from faces to elements)  
    [Zqle,Zule,zqfe,zufe,Ke,fe] = hdg_ConvDiff_ElementalMatrices(refElem,refFace,Xe,pElem,matElem,hdg.tau,faceInfo,ctt,problemParams);
    
    % Store the local solver
    local(iElem).Zql = Zqle;
    local(iElem).Zul = Zule;
    local(iElem).zqf = zqfe;
    local(iElem).zuf = zufe;
    
    % Assembly
    [indexGlobalFace, nOfFaceDOF] = hdgElemToFaceIndex(mesh.indexTf,refFace,hdg.faceInfo(iElem),refElem(pElem).nOfFaces,ctt.nOfComponents);
    nOfFaceDOF2 = nOfFaceDOF^2;
    currentIndex = indexIni + (1:nOfFaceDOF);
    currentIndex2 = indexIni + (1:nOfFaceDOF2);
    for i=1:nOfFaceDOF
        mat.i(currentIndex) = indexGlobalFace(i);
        mat.j(currentIndex) = indexGlobalFace;
        currentIndex = currentIndex + nOfFaceDOF;
    end
    mat.Kij(currentIndex2) = reshape(Ke',1,nOfFaceDOF2);
    
    % RHS
    f(indexGlobalFace) = f(indexGlobalFace) + fe;
    
    % Update indices
    indexIni = indexIni + nOfFaceDOF2;
end

K = sparse(mat.i,mat.j,mat.Kij);

% Solution
uhat = zeros(hdg.nDOFglobal*ctt.nOfComponents,1);
uhat(hdg.vDOFtoSolve) = K(hdg.vDOFtoSolve,hdg.vDOFtoSolve)\f(hdg.vDOFtoSolve);
