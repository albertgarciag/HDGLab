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
% [ug,qg] = hdg_Poisson_LocalPostprocessInterpolation(shapeFunctions,u,q,nOfNodes,Te,nsd)
%
% Inputs:
%
% shapeFunctions: element shape functions
% u:              solution on mesh elements (primal variable)
% q:              gradient of the solution on mesh elements (mixed variable)
% nOfNodes:       number of element nodes
% Te:             element connectivity
% nsd:            number of spatial dimensions
%
% Outputs:
%
% ug:      solution on gauss points of higher order element 
% qg:      gradient of the solution on gauss points of higher order element 
%
% See also: hdg_Poisson_LocalPostprocess, hdg_Poisson_LocalPostprocessElemMat

function [ug,qg] = hdg_Poisson_LocalPostprocessInterpolation(shapeFunctions,u,q,nOfNodes,Te,nsd)

nsdNOfPoints = nsd*nOfNodes;

ug = shapeFunctions*u(Te);

qg = zeros(nsd*nOfNodes,1);
for iNsd=1:nsd
    qIndex = nsd*Te-nsd+iNsd;
    qg(iNsd:nsd:nsdNOfPoints) = shapeFunctions*q(qIndex);
end