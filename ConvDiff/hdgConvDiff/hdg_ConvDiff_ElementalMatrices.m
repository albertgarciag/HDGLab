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
% [Zql,Zul,zqf,zuf,Ke,fe] = hdg_Poisson_ElementalMatrices(refElem,refFace,Xe,pElem,matElem,tau,faceInfo,ctt,problemParams)
%
% Inputs:
%
% refElem:          data structure for the reference element
% refFace:          data structure for the reference face
% Xe:               element nodal coordinates
% pElem:            degree of approximation for the primal variable
% matElem:          material index
% tau:              stabilisation parameter
% faceInfo:         data structure for the faces of an element
% ctt:              contains flags for boundary conditions and number of 
%                   components of solution field
% problemParams:    problem specific parameters
%
% Outputs:
%
% uhat:             solution on mesh faces (hybrid variable)
% local:            local solver (from faces to elements)  
%
% See also: hdg_Poisson_ElementalMatrices, hdg_Poisson_LocalProblem
%           hdgElemToFaceIndex
function [Zql,Zul,zqf,zuf,Ke,fe] = hdg_ConvDiff_ElementalMatrices(refElem,refFace,Xe,pElem,matElem,tau,tau_a,tau_d,faceInfo,ctt,problemParams)
kappa = problemParams.conductivity(matElem);
nOfFaces = refElem(pElem).nOfFaces;
tau_d = repmat(tau_d,1,nOfFaces);
tau_a = repmat(tau_a,1,nOfFaces);
tau = repmat(tau,1,nOfFaces);
%tau_d = repmat(tau_d,1,nOfFaces);
nsd = refElem(pElem).nsd;
nOfElementNodes = refElem(pElem).nOfNodes;
ndofU = nOfElementNodes;
ndofQ = ndofU*nsd;
ndofUHat = 0;
for iFace = 1:nOfFaces
    pHat = faceInfo.pHat(iFace);
    ndofUHat = ndofUHat + refFace(pHat,pHat).nOfNodes;
end
% Initialisation
fq = zeros(ndofQ,1);
fl = zeros(ndofUHat,1);
Auu = zeros(ndofU,ndofU);
Auq = zeros(ndofU,ndofQ);
Aul = zeros(ndofU,ndofUHat);
Alu = zeros(ndofUHat,ndofU);
Aqq = zeros(ndofQ,ndofQ);
Aql = zeros(ndofQ,ndofUHat);
All = zeros(ndofUHat,ndofUHat);
% Element computation -----------------------------------------------------
[N, dNx, wXY, Xg] = gaussElemCartesianInfo(refElem(pElem), Xe);
Nw = bsxfun(@times, N, wXY');
minusNN = -Nw*N';
for iNsd = 1:nsd
    vElem = iNsd:nsd:ndofQ;
    Auq(:,vElem)     = sqrt(kappa)*Nw*dNx(:,:,iNsd)';
    Aqq(vElem,vElem) = minusNN;
end
sourceTerm = convdiff_source(Xg,problemParams,matElem,nsd,problemParams.example);
fu = Nw*sourceTerm;
% Compute C = -dN * a * N' * w
a = convdiff_convection(Xg,problemParams,matElem,nsd,problemParams.example);
dNxa = zeros(size(dNx,[1 2]));
for iNsd=1:nsd
    dNxa = dNxa + bsxfun(@times,dNx(:,:,iNsd),a(:,iNsd)');
end
C = -dNxa*Nw';
% Face computation --------------------------------------------------------
indexFlip = zeros(1, ndofUHat);
indexFaceIni = 1;
nOfNodesPreviousFaces = 0;
for iFace = 1:nOfFaces
    iNumF = faceInfo.localNumFlux(iFace);
    iPerm = faceInfo.permutations(iFace);
    pHat  = faceInfo.pHat(iFace);
    
    nOfFaceNodesUhat = refFace(pHat,pHat).nOfNodes;
    indexFaceEnd = indexFaceIni + nOfFaceNodesUhat - 1;
    indexFaceV = indexFaceIni:indexFaceEnd;
    
    faceNodes = getNodesFace(iFace,0,refElem(pElem));
    Xf = Xe(faceNodes,:);
    
    faceNodesHDG = getNodesFaceHDG(iFace,iPerm,refElem(pHat));
    lastNodePreviousFaces = max(indexFlip);
    firstNodeNewFace = min(faceNodesHDG);
    indexFlip(indexFaceV) = faceNodesHDG + lastNodePreviousFaces + 1 - firstNodeNewFace;
    
    Nu = refFace(pElem,pHat).shapeFunctions;
    [N, NHat, n, ~, Xg, wXY] = gaussFaceCartesianInfo(Nu, refFace(pHat,pHat), Xf); 
    % calculo de conveccion en la frontera
    ahat = convdiff_convection(Xg,problemParams,matElem,nsd,problemParams.example);
    nahat = dot(n,ahat,2); 
    taunahat = tau(iFace) - nahat;
    
    if iNumF==ctt.iBC_Interior     
        Auu(faceNodes,faceNodes) = Auu(faceNodes,faceNodes) + N*bsxfun(@times, N', tau(iFace)*wXY);
        Aul(faceNodes,indexFaceV) = Aul(faceNodes,indexFaceV) + N*bsxfun(@times, NHat', taunahat.*wXY); % tau-n'a
        Alu(indexFaceV,faceNodes) = Alu(indexFaceV,faceNodes) + NHat*bsxfun(@times,N', tau(iFace)*wXY);
        All(indexFaceV,indexFaceV) = All(indexFaceV,indexFaceV) - NHat*bsxfun(@times, NHat', tau(iFace)*wXY);
        
        nNodesNsd = nsd*faceNodes;
        for iNsd = 1:nsd
            kNsd = nsd-iNsd;
            Aql(nNodesNsd-kNsd,indexFaceV) = Aql(nNodesNsd-kNsd,indexFaceV) + N*bsxfun(@times, NHat', sqrt(kappa)*wXY.*n(:,iNsd));
        end
    elseif iNumF==ctt.iBC_Dirichlet
        Auu(faceNodes,faceNodes) = Auu(faceNodes,faceNodes) + N*bsxfun(@times, N', tau(iFace)*wXY);
        uD = convdiff_Dirichlet(Xg,problemParams,matElem,nsd,problemParams.example);
        fu(faceNodes) = fu(faceNodes) + N*(uD.*wXY.*taunahat);  % tau-n'a
                
        nNodesNsd = nsd*faceNodes;
        for iNsd = 1:nsd
            kNsd = nsd-iNsd;
            fq(nNodesNsd-kNsd) = fq(nNodesNsd-kNsd) + N*(sqrt(kappa)*uD.*n(:,iNsd).*wXY);
        end        
    elseif iNumF==ctt.iBC_Neumann
        Auu(faceNodes,faceNodes) = Auu(faceNodes,faceNodes) + N*bsxfun(@times, N', tau(iFace)*wXY);
        Aul(faceNodes,indexFaceV) = Aul(faceNodes,indexFaceV) + N*bsxfun(@times, NHat', taunahat.*wXY); % tau-n'a
        % different matrices depending on the definition of global problem
        if problemParams.totalFluxNeumann==1
            Alu(indexFaceV,faceNodes) = Alu(indexFaceV,faceNodes) + NHat*bsxfun(@times,N', tau(iFace)*wXY); 
            All(indexFaceV,indexFaceV) = All(indexFaceV,indexFaceV) - NHat*bsxfun(@times, NHat', taunahat.*wXY); % tau-n'a
        elseif problemParams.totalFluxNeumann==0
            Alu(indexFaceV,faceNodes) = Alu(indexFaceV,faceNodes) + NHat*bsxfun(@times,N', tau_d(iFace)*wXY); 
            All(indexFaceV,indexFaceV) = All(indexFaceV,indexFaceV) - NHat*bsxfun(@times, NHat', tau_d(iFace)*wXY); 
        end
        t = convdiff_Neumann(Xg,n,problemParams,matElem,nsd,problemParams.example);
        fl(indexFaceV) = fl(indexFaceV) - NHat*(t.*wXY);
        
        nNodesNsd = nsd*faceNodes;
        for iNsd = 1:nsd
            kNsd = nsd-iNsd;
            Aql(nNodesNsd-kNsd,indexFaceV) = Aql(nNodesNsd-kNsd,indexFaceV) + N*bsxfun(@times, NHat', sqrt(kappa)*wXY.*n(:,iNsd));
        end
    end
    
    % Global indexing
    indexFaceIni = indexFaceEnd + 1;
    nOfNodesPreviousFaces = nOfNodesPreviousFaces + nOfFaceNodesUhat;
end
% Elemental mapping -------------------------------------------------------
Z = [Auu+C Auq; Auq' Aqq]\[Aul fu; Aql fq];
vU = 1:ndofU;
vQ = ndofU+1:ndofU+ndofQ;
Zul = Z(vU,1:ndofUHat);
zuf = Z(vU,ndofUHat+1);
Zql = Z(vQ,1:ndofUHat);
zqf = Z(vQ,ndofUHat+1);
% Flipping due to the different numbering of internal face nodes when
% seen from left or right element
Zql = Zql(:,indexFlip);
Zul = Zul(:,indexFlip);
Alq = Aql(:,indexFlip)';
Alu = Alu(indexFlip,:); %canvi (no estic segur del que estic fent)
All = All(indexFlip,indexFlip);
fl  = fl(indexFlip);
% Elemental matrices
Ke = Alq*Zql + Alu*Zul + All;
fe = fl - (Alq*zqf + Alu*zuf);