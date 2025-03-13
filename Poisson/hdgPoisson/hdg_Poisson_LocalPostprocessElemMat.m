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
% [Ke,Bq,intUStar,intU] = hdg_Poisson_LocalPostprocessElemMat(refElemStar,Xe,ug,qg,kappa)
%
% Inputs:
%
% refElemStar:  element shape functions
% Xe:           element nodal coordinates
% ug:           solution on gauss points of higher order element 
% qg:           gradient on gauss points of higher order element 
% kappa:        heat conductivity
%
% Outputs:
%
% Ke:           elemental stiffness matrix
% Bq:           integral of the right hand side
% intUStar:     integral of shape functions of higher order element 
% intU:         element integral of the primal variable
%
% See also: hdg_Poisson_LocalPostprocess, hdg_Poisson_LocalPostprocessInterpolation

function [Ke,Bqe,intUStar,intU] = hdg_Poisson_LocalPostprocessElemMat(refElemStar,Xe,ug,qg,kappa)

nsd = refElemStar.nsd;
nOfNodes = refElemStar.nOfNodes;
nsdNOfNodes = nsd*nOfNodes;

% Initialisation
Ke = zeros(nOfNodes);
Be = zeros(nOfNodes,nsdNOfNodes);

[N, dNx, wXY, ~] = gaussElemCartesianInfo(refElemStar, Xe);
dNw = bsxfun(@times, dNx, wXY');
for iNsd=1:nsd
    vElem = iNsd:nsd:nsdNOfNodes;
    Ke = Ke + kappa*dNx(:,:,iNsd)*dNw(:,:,iNsd)';
    Be(:,vElem) = - sqrt(kappa)*dNw(:,:,iNsd)*N';
end
intUStar = N*wXY;
    
Bqe = Be*qg;
intU = intUStar'*ug;  