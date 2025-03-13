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
% [N, NHat, n, t, Xg, wXY] = gaussFaceCartesianInfo(Nu, refFace, Xf)
%
% Inputs:
%
% Nu:       face shape functions and derivaties for the primal variable
% refElem:  data structure for the reference face
% Xf:       face nodal coordinates
%
% Outputs:
%
% N:        face shape functions for the primal variable
% NHat:     face shape functions for the hybrid variable
% n:        outward unit normal vector
% t:        tangential unit vector(s)
% Xg:       face Gauss points with Cartesian coordinates
% wXY:      face Gauss weights with Cartesian coordinates
%
% See also: hdg_Poisson_ElementalMatrices, getRefData


function [N, NHat, n, t, Xg, wXY] = gaussFaceCartesianInfo(Nu, refFace, Xf)

% Shape functions for the element solution
N = Nu(:,:,1)';

% Shape functions for the face solution
NHat = refFace.shapeFunctions(:,:,1)';

% Outward unit normal
if refFace.nsd==1
    t = Nu(:,:,2)*Xf;     
    n = [t(:,2), -t(:,1)];
    
    normN = sqrt(n(:,1).^2 + n(:,2).^2);
    n = [n(:,1)./normN, n(:,2)./normN];
    
    normT = sqrt(t(:,1).^2 + t(:,2).^2);
    t = [t(:,1)./normT, t(:,2)./normT];
elseif refFace.nsd==2
    t1 = Nu(:,:,2)*Xf; 
    t2 = Nu(:,:,3)*Xf; 
    n = cross3DColumn(t1,t2);
   
    normN = sqrt(n(:,1).^2 + n(:,2).^2 + n(:,3).^2);
    n = [n(:,1)./normN, n(:,2)./normN, n(:,3)./normN];
    
    normT1 = sqrt(t1(:,1).^2 + t1(:,2).^2 + t1(:,3).^2);
    normT2 = sqrt(t2(:,1).^2 + t2(:,2).^2 + t2(:,3).^2);    
    t = [t1(:,1)./normT1, t1(:,2)./normT1, t1(:,3)./normT1, t2(:,1)./normT2, t2(:,2)./normT2, t2(:,3)./normT2];
else
    err('gaussFaceCartesianInfo: wrong nsd')
end

% Gauss points (for traction evaluation)
Xg = N'*Xf;

% Weight
wXY = normN.*refFace.gaussWeights;