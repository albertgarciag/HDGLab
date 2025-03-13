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
% uErr = computeErrorL2norm(mesh,refElem,u,problemParams,nOfComponents,fAnalytical,uStar)
%
% Inputs:
%
% mesh:             data structure for the mesh
% refElem:          data structure for the reference element
% u:                primal variable on mesh elements
% problemParams:    problem specific parameters
% nOfComponents:    number of components of solution field
% fAnalytical:      function handle to analytical solution
% uStar:            super convergent solution on mesh elements 
%
% Outputs:
%
% uErr:             data structure for the solution error
%
% See also: computeErrorL2normElem, hdg_Poisson_LocalPostprocess
%           hdg_Poisson_LocalProblem

function uErr = computeErrorL2norm(mesh,refElem,u,problemParams,nOfComponents,fAnalytical,uStar)

TOLrel = 1e-6;

% Compute error of u ------------------------------------------------------
uErr.L2 = zeros(1,nOfComponents);
normUaL2 = zeros(1,nOfComponents);
for iElem = 1:mesh.nOfElements
    Te = getElemConnectivity(mesh, iElem);
    Xe = mesh.X(Te,:);
    pElem = mesh.pElem(iElem);
    ue = u(Te,:);
    
    [eL2, uL2] = computeErrorL2normElem(refElem(pElem),Xe,ue,fAnalytical,problemParams,mesh.matElem(iElem));
    
    uErr.L2 = uErr.L2 + eL2;
    normUaL2 = normUaL2 + uL2;
end

uErr.L2 = sqrt(uErr.L2);
normUaL2 = sqrt(normUaL2);
if normUaL2>TOLrel
    uErr.L2 = uErr.L2./normUaL2;
end

% Compute error of uStar (if provided) ------------------------------------
% In addition the L2 norm of the difference between u and uStar is computed
% This is useful to provide an error indicator that can drive an automatic 
% adaptivity process

if ~isempty(uStar)
    uErr.L2Star = zeros(1,nOfComponents);
    uErr.L2vsStarElem = zeros(mesh.nOfElements,nOfComponents);
    normUrefL2Elem = zeros(mesh.nOfElements,nOfComponents);
    normUaL2 = zeros(1,nOfComponents);
    
    indexIni = 1;
    for iElem = 1:mesh.nOfElements
        Te = getElemConnectivity(mesh, iElem);
        Xe = mesh.X(Te,:);
        pElem = mesh.pElem(iElem);
        pElemStar = pElem + 1;
        
        % Nodal distribution for p+1
        XeStar = refElem(pElem).shapeFunctionsNodesPPp1*Xe;
        NAtGaussStar = refElem(pElem).shapeFunctionsGaussPPp1;
        
        nOfDOFUStar = refElem(pElemStar).nOfNodes;
        indexEnd = indexIni + nOfDOFUStar - 1;
        
        ue = u(Te,:);
        ueStar = uStar(indexIni:indexEnd,:);
        
        [eL2, uL2] = computeErrorL2normElem(refElem(pElemStar),XeStar,ueStar,fAnalytical,problemParams,mesh.matElem(iElem));
        uErr.L2Star = uErr.L2Star + eL2;
        normUaL2 = normUaL2 + uL2;
        
        [dL2, uRefL2] = computeDiffL2normElem(NAtGaussStar,ue,refElem(pElemStar),XeStar,ueStar);
        uErr.L2vsStarElem(iElem,:) = dL2;
        normUrefL2Elem(iElem,:) = uRefL2;
        
        indexIni = indexEnd + 1;
    end
    
    uErr.L2Star = sqrt(uErr.L2Star);
    normUaL2 = sqrt(normUaL2);
    if normUaL2>TOLrel
        uErr.L2Star = uErr.L2Star./normUaL2;
    end
    
    uErr.L2vsStarElem = sqrt(uErr.L2vsStarElem);
    normUrefL2Elem = sqrt(normUrefL2Elem);
    if all(normUrefL2Elem>TOLrel)
        uErr.L2vsStarElem = uErr.L2vsStarElem./normUrefL2Elem;
    end
end


