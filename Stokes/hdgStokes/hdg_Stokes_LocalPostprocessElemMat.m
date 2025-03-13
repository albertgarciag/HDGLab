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
% [Ke, BLe, intUStar, intU] = hdg_Stokes_LocalPostprocessElemMat(refElemStar, Xe, ug, Lg, nu, problemParams)
%
% Inputs:
%
% refElemStar:  element shape functions
% Xe:           element nodal coordinates
% ug:           primal velocity at Gauss points of higher order element 
% Lg:           mixed variable at gauss points of higher order element 
% nu:           viscosity
%
% Outputs:
%
% Ke:           elemental stiffness matrix
% BLe:          integral of the right hand side
% intUStar:     integral of shape functions of higher order element 
% intU:         element integral of the primal velocity
%

function [Ke, BLe, intUStar, intU] = hdg_Stokes_LocalPostprocessElemMat(refElemStar, Xe, ug, Lg, nu, problemParams)

nsd = refElemStar.nsd;
nOfNodes = refElemStar.nOfNodes;
nsd2 = nsd^2;
nsdNOfNodes = nOfNodes*nsd;
nsd2NOfNodes = nOfNodes*nsd2;

% Initialisation
Ke = zeros(nsdNOfNodes);
Be = zeros(nsdNOfNodes,nsd2NOfNodes);
intUStar = zeros(nsdNOfNodes,nsd);

% Element contribution
[N, dNx, wXY, ~] = gaussElemCartesianInfo(refElemStar, Xe);
if problemParams.axi==1
    axiCoefElem = 2*pi*N'*Xe(:,2);
    dNw = bsxfun(@times, dNx, (axiCoefElem.*wXY)');
    NwT = N*(axiCoefElem.*wXY);
else
    dNw = bsxfun(@times, dNx, wXY);
    NwT = N*wXY;
end
for iNsd=1:nsd
    vElem = iNsd:nsd:nsdNOfNodes;
    %intUStar(vElem,iNsd) = N*wXY;
    intUStar(vElem,iNsd) = NwT;

    for jNsd = 1:nsd
       wElem = (iNsd-1)*nsd+jNsd:nsd2:nsd2NOfNodes;
       Be(vElem,wElem) = -sqrt(nu)*dNw(:,:,jNsd)*N';
       Ke(vElem,vElem) = Ke(vElem,vElem) + sqrt(nu)*dNx(:,:,jNsd)*dNw(:,:,jNsd)';
    end
end
    
BLe = Be*Lg;
intU = intUStar'*ug;  
