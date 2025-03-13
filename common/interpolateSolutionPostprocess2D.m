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
% visual = interpolateSolutionPostprocess2D(mesh, refElem, postproc, visual, colorField, fieldType)
%
% Inputs:
%
% mesh:             data structure for the mesh
% refElem:          data structure for the reference element
% postproc:          data structure for postprocess
% visual:         data structure for visualisation
% colorField:       elemental or nodal field values
% fieldType:        type of field: 'nodal' or 'elemental'
%
% Outputs:
%
% visual:         data structure for visualisation
%                   adds field U
%
% See also: buildSubmeshPostprocess2D

function visual = interpolateSolutionPostprocess2D(mesh, refElem, postproc, visual, colorField, fieldType)

if isempty(colorField) && strcmp(fieldType,'nodal')
    colorField = zeros(mesh.nOfNodes,1);
elseif isempty(colorField) && strcmp(fieldType,'elemental')
    colorField = zeros(mesh.nOfElements,1);
end

nOfComponents = size(colorField,2);

% Initialisation
nOfNodesSubMesh = size(visual.X,1);
visual.U = zeros(nOfNodesSubMesh, nOfComponents);

% Populate
iNodeIni = 1;
iNodePlotIni = 1;
for iElem = 1:mesh.nOfElements
    % Element info
    p = mesh.pElem(iElem);
    Te = mesh.indexT(iElem,1):mesh.indexT(iElem,2);
    
    % Submesh for plotting and mesh nodes
    iNodeEnd = iNodeIni + postproc.nOfNodesPlot - 1;
    iNodePlotEnd = iNodePlotIni + refElem(p).nOfNodes - 1;
    if strcmp(fieldType,'nodal')
        N = postproc.elem(p).N(:,:,1);
        interpSol = N*colorField(Te,:);
    elseif strcmp(fieldType,'elemental')
        interpSol = colorField(iElem)*ones(postproc.nOfNodesPlot,nOfComponents);
    end
    
    visual.U(iNodeIni:iNodeEnd,:) = interpSol;
    
    iNodeIni = iNodeEnd + 1;
    iNodePlotIni = iNodePlotEnd + 1;
end