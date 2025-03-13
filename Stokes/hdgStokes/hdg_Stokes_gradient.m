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
% gradU = hdg_Stokes_gradient(mesh, L, problemParams)
%
% Inputs:
%
% mesh:             data structure for the mesh
% L:                mixed variable in mesh elements
% problemParams:    problem specific parameters
%
% Outputs:
%
% gradU:            gradient in mesh elements 
%

function gradU = hdg_Stokes_gradient(mesh, L, problemParams)

gradU = zeros(size(L));
for iElem=1:mesh.nOfElements
    Te = getElemConnectivity(mesh, iElem);
    iMat = mesh.matElem(iElem);
    gradU(Te,:) = -L(Te,:)/sqrt(problemParams.viscosity(iMat));
end