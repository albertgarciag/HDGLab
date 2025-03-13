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
% uNew = hdgReshape(u,nOfComponents)
%
% Inputs:
%
% u:             vector of nodal values, 
%                dimension [nOfNodes x nOfComponents, 1]
% nOfComponents: number of  components of solution field
%
% Outputs:
%
% uNew:          matrix of nodal values
%                dimension [nOfNodes, nOfComponents]
%
% See also: hdg_Poisson_GlobalSystem, hdg_Poisson_LocalProblem

function uNew = hdgReshape(u,nOfComponents)

nOfUnknowns = numel(u);

uNew = zeros(nOfUnknowns/nOfComponents, nOfComponents);

for iComp = 1:nOfComponents
    uNew(:,iComp) = u(iComp:nOfComponents:nOfUnknowns);
end