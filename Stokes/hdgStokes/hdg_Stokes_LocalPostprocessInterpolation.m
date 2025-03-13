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
% [ug, Lg] = hdg_Stokes_LocalPostprocessInterpolation(shapeFunctions, u, L, nOfNodes, Te, nsd)
%
% Inputs:
%
% shapeFunctions: element shape functions
% u:              primal velocity in mesh elements
% L:              mixed variable in mesh elements
% nOfNodes:       number of element nodes
% Te:             element connectivity
% nsd:            number of spatial dimensions
%
% Outputs:
%
% ug:      velocity at Gauss points of higher order element 
% Lg:      mixed variable at Gauss points of higher order element 
%

function [ug, Lg] = hdg_Stokes_LocalPostprocessInterpolation(shapeFunctions, u, L, nOfNodes, Te, nsd)

nsd2 = nsd^2;
nsdNOfNodes = nOfNodes*nsd;
nsd2NOfNodes = nOfNodes*nsd2;

ug = zeros(nsdNOfNodes,1);
for iNsd = 1:nsd
    uIndex = nsd*Te-nsd+iNsd;
    ug(iNsd:nsd:nsdNOfNodes) = shapeFunctions*u(uIndex);
end

Lg = zeros(nsd2NOfNodes,1);
for iNsd = 1:nsd2
    LIndex = nsd2*Te-nsd2+iNsd;
    Lg(iNsd:nsd2:nsd2NOfNodes) = shapeFunctions*L(LIndex);
end