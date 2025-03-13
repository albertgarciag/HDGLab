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
% uStar = hdg_Poisson_LocalPostprocess(mesh,refElem,u,q,problemParams,nOfComponents)
%
% Inputs:
%
% mesh:          data structure for the mesh
% refElem:       data structure for the reference element
% u:             primal variable on mesh elements 
% q:             mixed variable on mesh elements
% problemParams: problem specific parameters
% nOfComponents: number of  components of solution field
%
% Outputs:
%
% uStar:         super convergent solution on mesh elements 
%
% See also: hdg_Poisson_GlobalSystem, hdg_Poisson_LocalProblem
%           hdg_Poisson_LocalPostprocessInterpolation
%           hdg_Poisson_LocalPostprocessElemMat

function uStar = hdg_Poisson_LocalPostprocess(mesh,refElem,u,q,problemParams,nOfComponents)

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
    
    NXeStar = refElem(pElem).shapeFunctionsNodesPPp1;
    XeStar = NXeStar*Xe;
    
    iMat = mesh.matElem(iElem);
    kappa = problemParams.conductivity(iMat);
    
    [ug,qg] = hdg_Poisson_LocalPostprocessInterpolation(NXeStar,u,q,refElem(pElemStar).nOfNodes,Te,mesh.nsd);
    [Ke,Bqe,intUStar,intU] = hdg_Poisson_LocalPostprocessElemMat(refElem(pElemStar),XeStar,ug,qg,kappa);
    
    % Constraint with Lagrange multipliers
    K = [Ke intUStar; intUStar' 0];
    f = [Bqe; intU];
    
    nOfDOFUStar = refElem(pElemStar).nOfNodes*nOfComponents;
    indexEnd = indexIni + nOfDOFUStar - 1;
    
    % Elemental solution
    sol = K\f;
    uStar(indexIni:indexEnd) = sol(1:end-1);
    
    % Indexing
    indexIni = indexEnd + 1;
end