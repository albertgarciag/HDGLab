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
%     [Zul, Zpl, ZLl, zuf, zpf, zLf, zur, zpr, zLr, Ke, fe, intPressure] = 
%        hdg_Stokes_ElementalMatrices(refElem, refFace, Xe, pElem, matElem, tau, faceInfo, ctt, problemParams, isPureDirichlet, intPressure);
%
% Inputs:
%
% refElem:            data structure for the reference element
% refFace:            data structure for the reference face
% Xe:                 element nodal coordinates
% pElem:              degree of approximation for the primal variable
% matElem:            material index
% tau:                stabilisation parameter
% faceInfo:           data structure for the faces of an element
% ctt:                contains flags for boundary conditions and number of 
%                     components of solution field
% problemParams:      problem specific parameters
% isPureDirichlet:    boolean variable: 1 if Dirichlet BC on all the external
%                     boundaries, 0 otherwise
% isAnySlipFace:      boolean variable: 1 if at least one external boundary 
%                     has a slip condition, 0 otherwise
% intPressure:        cumulative value of the integral of the pressure on the domain
%
% Outputs:
%
% Zul:                elemental contribution of uhat to the computation of u
% Zpl:                elemental contribution of uhat to the computation of p
% ZLl:                elemental contribution of uhat to the computation of L
% zuf:                elemental contribution of the source term to the
%                     computation of u 
% zpf:                elemental contribution of the source term to the
%                     computation of p 
% zLf:                elemental contribution of the source term to the 
%                     computation of L 
% zur:                elemental contribution of rho to the computation of u
% zpr:                elemental contribution of rho to the computation of p
% zLr:                elemental contribution of rho to the computation of L 
% Ke:                 elemental contribution to the matrix of the global problem
% fe:                 elemental contribution to the RHS of the global problem
% intPressure:        cumulative value of the integral of the pressure on the domain
%

function [Zul, Zpl, ZLl, zuf, zpf, zLf, zur, zpr, zLr, Ke, fe, intPressure] = ...
    hdg_Stokes_ElementalMatrices(refElem, refFace, Xe, pElem, matElem, tau, faceInfo, ctt, problemParams, isPureDirichlet, isAnySlipFace, intPressure)

nu = problemParams.viscosity(matElem);

nOfFaces = refElem(pElem).nOfFaces;
tau = repmat(tau*nu/problemParams.charLength,1,nOfFaces);

nsd = refElem(pElem).nsd;
nsd2 = nsd*nsd;

nOfElementNodes = refElem(pElem).nOfNodes;
ndofP = nOfElementNodes;
ndofU = nOfElementNodes*ctt.nOfComponents;
ndofL = nOfElementNodes*ctt.nOfComponents*nsd;

nOfNodesAllFaces = 0;
for iFace = 1:nOfFaces
    pHat = faceInfo.pHat(iFace);
    nOfNodesAllFaces = nOfNodesAllFaces + refFace(pHat,pHat).nOfNodes;
end
ndofUhat = nOfNodesAllFaces*ctt.nOfComponents;

% Initialisation of the elemental contributions
% Local problem: Mixed variable equation
ALL = zeros(ndofL,ndofL);
ALu = zeros(ndofL,ndofU);
ALl = zeros(ndofL,ndofUhat);
fL  = zeros(ndofL,1);
% Local problem: Momentum equation
Auu = zeros(ndofU,ndofU);
Aul = zeros(ndofU,ndofUhat);
fu  = zeros(ndofU,1);
% Local problem: Continuity equation
Apu = zeros(ndofP,ndofU);
Apl = zeros(ndofP,ndofUhat);
fp  = zeros(ndofP,1);
% Local problem: Solvability condition for uniqueness of pressure
Arp = zeros(ndofP,1);
% Global problem: Transmission condition equation
All = zeros(ndofUhat,ndofUhat);
fl  = zeros(ndofUhat,1);
% Global problem: Incompressibility constraint
Arl = zeros(1, ndofUhat);
fn = 0;
% Additional restriction for purely Dirichlet BC
ArlExtra = zeros(ndofP,1);
% Additional contributions if any of the faces features a slip BC
if isAnySlipFace
    AlL = zeros(ndofUhat,ndofL);
    Alu = zeros(ndofUhat,ndofU);
    Alp = zeros(ndofUhat,ndofP);
end

% Element computation -----------------------------------------------------
[N, dNx, wXY, Xg] = gaussElemCartesianInfo(refElem(pElem), Xe);
areaElem = sum(wXY);
Nw = bsxfun(@times, N, wXY');
NwN = Nw*N';

sourceTerm = stokes_source(Xg,problemParams,matElem,nsd,problemParams.example);
analyticalPressure = stokes_analyticalPressure(Xg,problemParams,matElem,nsd,problemParams.example);

for iNsd = 1:nsd2
    vElem = iNsd:nsd2:ndofL;
    ALL(vElem,vElem) = -NwN;
end

for iNsd = 1:nsd
    wElem = iNsd:nsd:ndofU;
    Apu(:,wElem) = dNx(:,:,iNsd)*Nw';
    fu(wElem) = Nw*sourceTerm(:,iNsd);
    for jNsd = 1:nsd
       zElem = (iNsd-1)*nsd+jNsd:nsd2:ndofL;
       ALu(zElem,wElem) = sqrt(nu)*dNx(:,:,jNsd)*Nw';
    end
end

Arp = N*wXY/areaElem;

% If purely Dirichlet BC
if isPureDirichlet
    ArlExtra = sum(Nw,2);
    intPressure = intPressure + sum(Nw*analyticalPressure,1);
end

% Face computation --------------------------------------------------------
indexFlip = zeros(1, ndofUhat);
indexFaceIni = 1;
nOfNodesPreviousFaces = 0;
areaFaces = 0;
for iFace = 1:nOfFaces
    iNumF = faceInfo.localNumFlux(iFace);
    iPerm = faceInfo.permutations(iFace);
    pHat  = faceInfo.pHat(iFace);
    
    nOfFaceNodesUhat = refFace(pHat,pHat).nOfNodes;
    ndofFaceUhat = nOfFaceNodesUhat*ctt.nOfComponents;
    indexFaceEnd = indexFaceIni + ndofFaceUhat - 1;
    indexFaceV = indexFaceIni:indexFaceEnd; 
    
    faceNodes = getNodesFace(iFace,0,refElem(pElem));
    Xf = Xe(faceNodes,:);
    
    faceNodesHDG = getNodesFaceHDG(iFace,iPerm,refElem(pHat));
    lastNodePreviousFaces = max(indexFlip)/nsd;
    firstNodeNewFace = min(faceNodesHDG);
    nodeVec = zeros(1,2*length(faceNodesHDG));
    for iNode = 1:nsd
        nodeVec(iNode:nsd:length(faceNodesHDG)*nsd-(nsd-iNode)) = nsd*faceNodesHDG - (nsd-iNode);
    end
    nodeVec = nodeVec + (lastNodePreviousFaces+1-firstNodeNewFace)*nsd;
    indexFlip(indexFaceV) = nodeVec;
    
    Nu = refFace(pElem,pHat).shapeFunctions;
    [N, NHat, n, t, Xg, wXY] = gaussFaceCartesianInfo(Nu, refFace(pHat,pHat), Xf);
    areaFaces = areaFaces + sum(wXY);
    
    nNodesNsd = faceNodes*nsd;
    nNodesNsd2 = faceNodes*nsd2;
    if iNumF==ctt.iBC_Interior     
%        Arp(faceNodes) = Arp(faceNodes) + N*wXY;

        for iNsd = 1:nsd
            kNsd = nsd-iNsd;
            Auu(nNodesNsd-kNsd,nNodesNsd-kNsd) = Auu(nNodesNsd-kNsd,nNodesNsd-kNsd) + N*bsxfun(@times, N', tau(iFace)*wXY);
            Aul(nNodesNsd-kNsd,indexFaceV(iNsd:nsd:ndofFaceUhat)) = Aul(nNodesNsd-kNsd,indexFaceV(iNsd:nsd:ndofFaceUhat)) + N*bsxfun(@times, NHat', tau(iFace)*wXY);
            Apl(faceNodes,indexFaceV(iNsd:nsd:ndofFaceUhat)) = Apl(faceNodes,indexFaceV(iNsd:nsd:ndofFaceUhat)) + N*bsxfun(@times, NHat', wXY.*n(:,iNsd));
            All(indexFaceV(iNsd:nsd:ndofFaceUhat),indexFaceV(iNsd:nsd:ndofFaceUhat)) = All(indexFaceV(iNsd:nsd:ndofFaceUhat),indexFaceV(iNsd:nsd:ndofFaceUhat)) - NHat*bsxfun(@times, NHat', tau(iFace)*wXY);
            Arl(1,indexFaceV(iNsd:nsd:ndofFaceUhat)) = Arl(1,indexFaceV(iNsd:nsd:ndofFaceUhat)) + (NHat*(wXY.*n(:,iNsd)))';
            
            if isAnySlipFace
                Alp(indexFaceV(iNsd:nsd:ndofFaceUhat),faceNodes) = Alp(indexFaceV(iNsd:nsd:ndofFaceUhat),faceNodes) + NHat*bsxfun(@times, N', wXY.*n(:,iNsd));
                Alu(indexFaceV(iNsd:nsd:ndofFaceUhat),nNodesNsd-kNsd) = Alu(indexFaceV(iNsd:nsd:ndofFaceUhat),nNodesNsd-kNsd) + NHat*bsxfun(@times, N', tau(iFace)*wXY);
            end
            
            for jNsd = 1:nsd
                mNsd = nsd2-(iNsd-1)*nsd-jNsd;
                ALl(nNodesNsd2-mNsd,indexFaceV(iNsd:nsd:ndofFaceUhat)) = ALl(nNodesNsd2-mNsd,indexFaceV(iNsd:nsd:ndofFaceUhat)) + N*bsxfun(@times, NHat', sqrt(nu)*wXY.*n(:,jNsd));
                
                if isAnySlipFace
                   AlL(indexFaceV(iNsd:nsd:ndofFaceUhat),nNodesNsd2-mNsd) = AlL(indexFaceV(iNsd:nsd:ndofFaceUhat),nNodesNsd2-mNsd) + NHat*bsxfun(@times, N', sqrt(nu)*wXY.*n(:,jNsd)); 
                end
            end
        end
    elseif iNumF==ctt.iBC_Dirichlet
        uD = stokes_dirichlet(Xg,problemParams,matElem,nsd,problemParams.example);
        
%        Arp(faceNodes) = Arp(faceNodes) + N*wXY;
        fp(faceNodes) = fp(faceNodes) + N*(sum(uD.*n,2).*wXY);
        fn = fn - sum(sum(uD.*n,2).*wXY,1);
        
        for iNsd = 1:nsd
            kNsd = nsd-iNsd;
            Auu(nNodesNsd-kNsd,nNodesNsd-kNsd) = Auu(nNodesNsd-kNsd,nNodesNsd-kNsd) + N*bsxfun(@times, N', tau(iFace)*wXY);
            fu(nNodesNsd-kNsd) = fu(nNodesNsd-kNsd) + N*(uD(:,iNsd)*tau(iFace).*wXY);
            
            for jNsd = 1:nsd
                mNsd = nsd2-(iNsd-1)*nsd-jNsd;
                fL(nNodesNsd2-mNsd) = fL(nNodesNsd2-mNsd) + N*(sqrt(nu)*uD(:,iNsd).*n(:,jNsd).*wXY);
            end
        end
    elseif iNumF==ctt.iBC_Neumann
        neumannTension = stokes_neumann(Xg,n,problemParams,matElem,nsd,problemParams.example);
        
%        Arp(faceNodes) = Arp(faceNodes) + N*wXY;

        for iNsd = 1:nsd
            kNsd = nsd-iNsd;
            Auu(nNodesNsd-kNsd,nNodesNsd-kNsd) = Auu(nNodesNsd-kNsd,nNodesNsd-kNsd) + N*bsxfun(@times, N', tau(iFace)*wXY);
            Aul(nNodesNsd-kNsd,indexFaceV(iNsd:nsd:ndofFaceUhat)) = Aul(nNodesNsd-kNsd,indexFaceV(iNsd:nsd:ndofFaceUhat)) + N*bsxfun(@times, NHat', tau(iFace)*wXY);
            Apl(faceNodes,indexFaceV(iNsd:nsd:ndofFaceUhat)) = Apl(faceNodes,indexFaceV(iNsd:nsd:ndofFaceUhat)) + N*bsxfun(@times, NHat', wXY.*n(:,iNsd));
            All(indexFaceV(iNsd:nsd:ndofFaceUhat),indexFaceV(iNsd:nsd:ndofFaceUhat)) = All(indexFaceV(iNsd:nsd:ndofFaceUhat),indexFaceV(iNsd:nsd:ndofFaceUhat)) - NHat*bsxfun(@times, NHat', tau(iFace)*wXY);
            Arl(1,indexFaceV(iNsd:nsd:ndofFaceUhat)) = Arl(1,indexFaceV(iNsd:nsd:ndofFaceUhat)) + (NHat*(wXY.*n(:,iNsd)))';
            fl(indexFaceV(iNsd:nsd:ndofFaceUhat)) = fl(indexFaceV(iNsd:nsd:ndofFaceUhat)) - NHat*(neumannTension(:,iNsd).*wXY);
            
            if isAnySlipFace
                Alp(indexFaceV(iNsd:nsd:ndofFaceUhat),faceNodes) = Alp(indexFaceV(iNsd:nsd:ndofFaceUhat),faceNodes) + NHat*bsxfun(@times, N', wXY.*n(:,iNsd));
                Alu(indexFaceV(iNsd:nsd:ndofFaceUhat),nNodesNsd-kNsd) = Alu(indexFaceV(iNsd:nsd:ndofFaceUhat),nNodesNsd-kNsd) + NHat*bsxfun(@times, N', tau(iFace)*wXY);
            end
            
            for jNsd = 1:nsd
                mNsd = nsd2-(iNsd-1)*nsd-jNsd;
                ALl(nNodesNsd2-mNsd,indexFaceV(iNsd:nsd:ndofFaceUhat)) = ALl(nNodesNsd2-mNsd,indexFaceV(iNsd:nsd:ndofFaceUhat)) + N*bsxfun(@times, NHat', sqrt(nu)*wXY.*n(:,jNsd));
                
                if isAnySlipFace
                   AlL(indexFaceV(iNsd:nsd:ndofFaceUhat),nNodesNsd2-mNsd) = AlL(indexFaceV(iNsd:nsd:ndofFaceUhat),nNodesNsd2-mNsd) + NHat*bsxfun(@times, N', sqrt(nu)*wXY.*n(:,jNsd)); 
                end
            end
        end 
    elseif iNumF==ctt.iBC_Slip        
        slipForce = stokes_slip(Xg,n,t,problemParams,matElem,nsd,problemParams.example);
        D = [n, problemParams.betaSlip*t];
        E = [problemParams.alphaSlip*n, t];
                
%        Arp(faceNodes) = Arp(faceNodes) + N*wXY;

        for iNsd = 1:nsd
            kNsd = nsd-iNsd;
            Auu(nNodesNsd-kNsd,nNodesNsd-kNsd) = Auu(nNodesNsd-kNsd,nNodesNsd-kNsd) + N*bsxfun(@times, N', tau(iFace)*wXY);
            Aul(nNodesNsd-kNsd,indexFaceV(iNsd:nsd:ndofFaceUhat)) = Aul(nNodesNsd-kNsd,indexFaceV(iNsd:nsd:ndofFaceUhat)) + N*bsxfun(@times, NHat', tau(iFace)*wXY);
            Apl(faceNodes,indexFaceV(iNsd:nsd:ndofFaceUhat)) = Apl(faceNodes,indexFaceV(iNsd:nsd:ndofFaceUhat)) + N*bsxfun(@times, NHat', wXY.*n(:,iNsd));
            Arl(1,indexFaceV(iNsd:nsd:ndofFaceUhat)) = Arl(1,indexFaceV(iNsd:nsd:ndofFaceUhat)) + (NHat*(wXY.*n(:,iNsd)))';
            fl(indexFaceV(iNsd:nsd:ndofFaceUhat)) = fl(indexFaceV(iNsd:nsd:ndofFaceUhat)) + NHat*(slipForce(:,iNsd).*wXY);
            
            for jNsd = 1:nsd
                mNsd = nsd2-(iNsd-1)*nsd-jNsd;
                ALl(nNodesNsd2-mNsd,indexFaceV(iNsd:nsd:ndofFaceUhat)) = ALl(nNodesNsd2-mNsd,indexFaceV(iNsd:nsd:ndofFaceUhat)) + N*bsxfun(@times, NHat', sqrt(nu)*wXY.*n(:,jNsd));
                
                pNsd = nsd-jNsd;
                Alu(indexFaceV(iNsd:nsd:ndofFaceUhat),nNodesNsd-pNsd) = Alu(indexFaceV(iNsd:nsd:ndofFaceUhat),nNodesNsd-pNsd) - NHat*bsxfun(@times, N', tau(iFace)*wXY.*E(:,nsd*(iNsd-1)+jNsd));
                All(indexFaceV(iNsd:nsd:ndofFaceUhat),indexFaceV(jNsd:nsd:ndofFaceUhat)) = All(indexFaceV(iNsd:nsd:ndofFaceUhat),indexFaceV(jNsd:nsd:ndofFaceUhat)) + NHat*bsxfun(@times, NHat', wXY.*(D(:,nsd*(iNsd-1)+jNsd)+tau(iFace)*E(:,nsd*(iNsd-1)+jNsd)));
                Alp(indexFaceV(iNsd:nsd:ndofFaceUhat),faceNodes) = Alp(indexFaceV(iNsd:nsd:ndofFaceUhat),faceNodes) - NHat*bsxfun(@times, N', wXY.*n(:,jNsd).*E(:,nsd*(iNsd-1)+jNsd));
                
                for lNsd = 1:nsd
                    qNsd = nsd2-(jNsd-1)*nsd-lNsd;
                    AlL(indexFaceV(iNsd:nsd:ndofFaceUhat),nNodesNsd2-qNsd) = AlL(indexFaceV(iNsd:nsd:ndofFaceUhat),nNodesNsd2-qNsd) - NHat*bsxfun(@times, N', sqrt(nu)*wXY.*n(:,lNsd).*E(:,nsd*(iNsd-1)+jNsd));
                end
            end
        end 
    end
    
    % Global indexing
    indexFaceIni = indexFaceEnd + 1;
    nOfNodesPreviousFaces = nOfNodesPreviousFaces + nOfFaceNodesUhat*ctt.nOfComponents;
end

% To ensure symmetry of the global matrix
%Arp = Arp/areaFaces;

% Elemental mapping -------------------------------------------------------
Z = [ALL                 ALu             zeros(ndofL,ndofP)  zeros(ndofL,1);
     ALu'                Auu             Apu'                zeros(ndofU,1);
     zeros(ndofP,ndofL)  Apu             zeros(ndofP)        Arp;
     zeros(1,ndofL)      zeros(1,ndofU)  Arp'                0]\[ALl                fL  zeros(ndofL,1); 
                                                                 Aul                fu  zeros(ndofU,1); 
                                                                 Apl                fp  zeros(ndofP,1); 
                                                                 zeros(1,ndofUhat)  0   1];


% Contribution of the local element
vL = 1:ndofL;
vU = ndofL+1:ndofL+ndofU;
vP = ndofL+ndofU+1:ndofL+ndofU+ndofP;
ZLl = Z(vL,1:ndofUhat);
Zul = Z(vU,1:ndofUhat);
Zpl = Z(vP,1:ndofUhat);
zLf = Z(vL,ndofUhat+1);
zuf = Z(vU,ndofUhat+1);
zpf = Z(vP,ndofUhat+1);
zLr = Z(vL,ndofUhat+2);
zur = Z(vU,ndofUhat+2);
zpr = Z(vP,ndofUhat+2);

% Flipping due to the different numbering of internal face nodes when
% seen from left or right element
ZLl = ZLl(:,indexFlip);
Zul = Zul(:,indexFlip);
Zpl = Zpl(:,indexFlip);
if isAnySlipFace
    AlL = AlL(indexFlip,:);
    Alu = Alu(indexFlip,:);
    Alp = Alp(indexFlip,:);
else
    AlL = ALl(:,indexFlip)';
    Alu = Aul(:,indexFlip)';
    Alp = Apl(:,indexFlip)';
end
All = All(indexFlip,indexFlip);    
fl  = fl(indexFlip);
Arl = Arl(:,indexFlip);

% Elemental matrices
if isPureDirichlet
    Ke = [All + (AlL*ZLl + Alu*Zul + Alp*Zpl),  AlL*zLr + Alu*zur + Alp*zpr;
          Arl,                                  0;
          ArlExtra'*Zpl                         ArlExtra'*zpr];
    fe = [fl - (AlL*zLf + Alu*zuf + Alp*zpf);
          fn; 
          ArlExtra'*zpf];
else
    Ke = [All + (AlL*ZLl + Alu*Zul + Alp*Zpl),  AlL*zLr + Alu*zur + Alp*zpr;
          Arl,                                  0];
    fe = [fl - (AlL*zLf + Alu*zuf + Alp*zpf);
          fn];
end