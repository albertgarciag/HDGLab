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
% order = nodeRenumberingMSH(pDegree, nsd)
%
% Inputs:
%
% pDegree:     polynomial degree of approximation
% nsd:         number of spatial dimensions
%
% Outputs:
%
% order:      renumbering of the elemental nodes to match Matlab order
%

function order = nodeRenumberingMSH(pDegree, nsd)

if nsd == 2 %%%%%%%%%%%%%%%%%%%%% TRI
   switch pDegree
       case 1
           order = [1, 2, 3];
       case 2
           order = [1, 2, 3, 4, 6, 5];
       case 3
           order = [1, 2, 3, 4, 5, 9, 10, 6, 8, 7];
       case 4
           order = [1, 2, 3, 4, 5, 6, 12, 13, 14, 7, 11, 15, 8, 10, 9];
       case 5
           order = [1, 2, 3, 4, 5, 6, 7, 15, 16, 19, 17, 8, 14, 21, 20, 9, 13, 18, 10, 12, 11];
   end
elseif nsd==3 %%%%%%%%%%%%%%%%%%% TET
   switch pDegree
       case 1
           order = [1, 4, 2, 3];
       case 2
           order = [1, 4, 2, 3, 8, 5, 10, 7, 9, 6];
       case 3
            order = [1, 4, 2, 3, 12, 11, 5, 18, 15, 6, 16, 10, 19, 13, 17, 20, 7, 9, 14, 8];
       case 4
           order = [1, 4, 2, 3, 16, 15, 14, 5, 26, 28, 20, 6, 27, 21, 7, 22, 13, 29, 30, 17, 23, 35, 32, 25, 33, 8, 12, 31, 18, 24, 34, 9, 11, 19, 10];
       case 5
           order = [1, 4, 2, 3, 20, 19, 18, 17, 5, 35, 40, 37, 25, 6, 38, 39, 26, 7, 36, 27, 8, 28, 16, 41, 44, 42, 21, 29, 53, 56, 47, 34, 54, 50, 31, 48, 9, 15, 46, 45, 22, 32, 55, 52, 33, 51, 10, 14, 43, 23, 30, 49, 11, 13, 24, 12];
   end
end