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
% uStar = hdg_Stokes_LocalPostprocess(mesh, refElem, u, L, problemParams, nOfComponents)
%
% Inputs:
%
% mesh:            data structure for the mesh
% refElem:         data structure for the reference element
% u:               primal velocity in mesh elements 
% L:               mixed variable in mesh elements
% problemParams:   problem specific parameters
% nOfComponents:   number of  components of solution field
%
% Outputs:
%
% uStar:           superconvergent postprocessed velocity in mesh elements 
%
% See also: hdg_Stokes_LocalPostprocessInterpolation,
%           hdg_Stokes_LocalPostprocessElemMat
%

function uStar = hdg_Stokes_LocalPostprocess(mesh, refElem, u, L, problemParams, nOfComponents)

% Initialisation
nOfDOFUStar = 0;
for iElem = 1:mesh.nOfElements
    pElem = mesh.pElem(iElem);
    pElemStar = pElem + 1;
    nOfElementNodesStar = refElem(pElemStar).nOfNodes;
    nOfDOFUStar = nOfDOFUStar + nOfElementNodesStar*nOfComponents;
end
uStar = zeros(nOfDOFUStar,1);

% Computation
indexIni = 1;
for iElem = 1:mesh.nOfElements
    Te = getElemConnectivity(mesh, iElem);
    Xe = mesh.X(Te,:);
    pElem = mesh.pElem(iElem);
    pElemStar = pElem + 1;
    
    nOfDOFElemUStar = refElem(pElemStar).nOfNodes*nOfComponents;
    
    NXeStar = refElem(pElem).shapeFunctionsNodesPPp1;
    XeStar = NXeStar*Xe;
    
    iMat = mesh.matElem(iElem);
    nu = problemParams.viscosity(iMat);
    
    [ug, Lg] = hdg_Stokes_LocalPostprocessInterpolation(NXeStar, u, L, refElem(pElemStar).nOfNodes, Te, mesh.nsd);
    [Ke, BLe, intUStar, intU] = hdg_Stokes_LocalPostprocessElemMat(refElem(pElemStar), XeStar, ug, Lg, nu, problemParams);
        
    indexEnd = indexIni + nOfDOFElemUStar - 1;
   
    % Introduce the constraint via a Lagrange multiplier
%     K = [Ke intUStar; intUStar' zeros(mesh.nsd)];
%     f = [BLe; intU];
%     sol = K\f;
    sol = [Ke intUStar; 
           intUStar' zeros(mesh.nsd)]\[BLe; 
                                       intU];
    uStar(indexIni:indexEnd) = sol(1:nOfDOFElemUStar);
    
    % Indexing
    indexIni = indexEnd + 1;
end
