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
% [refElem, refFace] = getRefData(pElem, nsd, optionNodes)
%
% Inputs:
%
% pElem:         a vector of dimension nOfElements with the degree of 
%                approximation to be used in each element
% nsd:           number of spatial dimensions
% optionNodes:   flag for Fekete (=1) or uniform (=0) nodal distributions
%
% Outputs:
%
% refElem:       data structure for the reference element
% refFace:       data structure for the reference face
%

function [refElem, refFace] = getRefData(pElem, nsd, optionNodes)

if nsd==2
    data.elem = load(sprintf('refElemTRI_opt%d', optionNodes));
    data.face = load(sprintf('refFaceTRI_opt%d', optionNodes));
elseif nsd==3
    data.elem = load(sprintf('refElemTET_opt%d', optionNodes));
    data.face = load(sprintf('refFaceTET_opt%d', optionNodes));
else
    err('getIsoRefData: wrong nsd')
end

pMax = max(pElem);

refElem = data.elem.refElem(1:pMax+1);
refFace = data.face.refFace(1:pMax+1,1:pMax+1);