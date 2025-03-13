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
% [detJ, invJ] = isoparametricInvJ3D(Xe, Nxi, Neta, Nzet)
%
% Inputs:
%
% Xe:            element nodal coordinates
% Nxi,Neta,Nzet: shape function derivatives in the reference element
%
% Outputs:
%
% detJ:          determinant of the Jacobian of the isoparametric mapping
% invJ:          inverse of the Jacobian of the isoparametric mapping
%
% See also: hdg_Poisson_ElementalMatrices, getRefData

function [detJ, invJ] = isoparametricInvJ3D(Xe, Nxi, Neta, Nzet)

J1 = Nxi*Xe;
J2 = Neta*Xe;
J3 = Nzet*Xe;

detJ = J1(:,1).*(J2(:,2).*J3(:,3) - J3(:,2).*J2(:,3)) + ...
       J1(:,2).*(J2(:,3).*J3(:,1) - J3(:,3).*J2(:,1)) + ...
       J1(:,3).*(J2(:,1).*J3(:,2) - J3(:,1).*J2(:,2));
   
invJ.a11 = ( ( J2(:,2).*J3(:,3) - J2(:,3).*J3(:,2) )./detJ )';
invJ.a12 = ( ( J1(:,3).*J3(:,2) - J1(:,2).*J3(:,3) )./detJ )';
invJ.a13 = ( ( J1(:,2).*J2(:,3) - J1(:,3).*J2(:,2) )./detJ )';
invJ.a21 = ( ( J2(:,3).*J3(:,1) - J2(:,1).*J3(:,3) )./detJ )';
invJ.a22 = ( ( J1(:,1).*J3(:,3) - J1(:,3).*J3(:,1) )./detJ )';
invJ.a23 = ( ( J1(:,3).*J2(:,1) - J1(:,1).*J2(:,3) )./detJ )';
invJ.a31 = ( ( J2(:,1).*J3(:,2) - J2(:,2).*J3(:,1) )./detJ )';
invJ.a32 = ( ( J1(:,2).*J3(:,1) - J1(:,1).*J3(:,2) )./detJ )';
invJ.a33 = ( ( J1(:,1).*J2(:,2) - J1(:,2).*J2(:,1) )./detJ )';