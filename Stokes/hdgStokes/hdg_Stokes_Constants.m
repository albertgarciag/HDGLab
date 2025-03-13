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
% ctt = hdg_Stokes_Constants(nsd)
%
% Input: 
%
% nsd:   number of spatial dimensions of the problem
%
% Output: 
%
% ctt:   contains flags for boundary conditions and number of components of
%        the velocity field
%
% The flag corresponds to the third column of the field mesh.extFaces
% It can be changed to match your mesh specification
%

function ctt = hdg_Stokes_Constants(nsd)

ctt.iBC_Interior  = 0;
ctt.iBC_Dirichlet = 1;
ctt.iBC_Neumann   = 2;
ctt.iBC_Slip      = 3;
ctt.iBC_Axi       = 4;

ctt.nOfComponents  = nsd;