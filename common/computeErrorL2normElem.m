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
% [eL2, uL2] = computeErrorL2normElem(refElem,Xe,ue,problemParams,matElem,fAnalytical)
%
% Inputs:
%
% refElem:          data structure for the reference element
% Xe:               element nodal coordinates
% ue:               primal variable on element
% problemParams:    problem specific parameters
% matElem:          material index
% fAnalytical:      function handle to analytical solution
%
% Outputs:
%
% uErr:             data structure for the solution error
%
% See also: computeErrorL2norm, gaussElemCartesianInfo

function [eL2, uL2] = computeErrorL2normElem(refElem,Xe,ue,fAnalytical,problemParams,matElem)

[N, ~, wXY, Xg] = gaussElemCartesianInfo(refElem, Xe);
uG = N'*ue;

uaG  = fAnalytical(Xg,problemParams,matElem,refElem.nsd,problemParams.example);

eL2 = wXY'*(uG - uaG).^2;
uL2 = wXY'*uaG.^2;