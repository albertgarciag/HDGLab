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
% [N, dNx, wXY, Xg] = gaussElemCartesianInfo(refElem, Xe)
%
% Inputs:
%
% refElem:  data structure for the reference element
% Xe:       element nodal coordinates
%
% Outputs:
%
% N:        shape functions 
% dNx:      derivatives with Cartesian coordinates
% wXY:      element Gauss weights with Cartesian coordinates
% Xg:       element Gauss points with Cartesian coordinates
%
% See also: hdg_Poisson_ElementalMatrices, getRefData
%           isoparametricInvJ2D, isoparametricInvJ3D

function [N, dNx, wXY, Xg] = gaussElemCartesianInfo(refElem, Xe)

N = refElem.shapeFunctions(:,:,1)';
dNRef = zeros(refElem.nOfNodes,refElem.nOfGauss,refElem.nsd);

for iNsd =1:refElem.nsd
    dNRef(:,:,iNsd) = refElem.shapeFunctions(:,:,iNsd+1)';
end

% Derivatives with Cartesian coordinates
dNx = zeros(refElem.nOfNodes,refElem.nOfGauss,refElem.nsd);
if refElem.nsd==2
    [detJ, invJ] = isoparametricInvJ2D(Xe, dNRef(:,:,1)', dNRef(:,:,2)');
    dNx(:,:,1) = bsxfun(@times, dNRef(:,:,1), invJ.a11) + bsxfun(@times, dNRef(:,:,2), invJ.a12);
    dNx(:,:,2) = bsxfun(@times, dNRef(:,:,1), invJ.a21) + bsxfun(@times, dNRef(:,:,2), invJ.a22);
elseif refElem.nsd==3
    [detJ, invJ] = isoparametricInvJ3D(Xe, dNRef(:,:,1)', dNRef(:,:,2)', dNRef(:,:,3)');
    dNx(:,:,1) = bsxfun(@times, dNRef(:,:,1), invJ.a11) + bsxfun(@times, dNRef(:,:,2), invJ.a12) + bsxfun(@times, dNRef(:,:,3), invJ.a13);
    dNx(:,:,2) = bsxfun(@times, dNRef(:,:,1), invJ.a21) + bsxfun(@times, dNRef(:,:,2), invJ.a22) + bsxfun(@times, dNRef(:,:,3), invJ.a23);
    dNx(:,:,3) = bsxfun(@times, dNRef(:,:,1), invJ.a31) + bsxfun(@times, dNRef(:,:,2), invJ.a32) + bsxfun(@times, dNRef(:,:,3), invJ.a33);
else
    err('gaussElemCartesianInfo: wrong nsd')
end

% Gauss points (for source evaluation)
Xg = N'*Xe;

% Weights in Cartesian coordinates
wXY = refElem.gaussWeights.*detJ;