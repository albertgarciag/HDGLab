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
% visual = interpolateSolutionPostprocess3D(mesh,refFace,postproc,visual,plotOpts,colorField,fieldType)
%
% Inputs:
%
% mesh:             data structure for the mesh
% refFace:          data structure for the reference face
% postproc:          data structure for postprocess
% visual:         data structure for visualisation
% plotOpts:         user options for visualisation
% colorField:       elemental or nodal field values
% fieldType:        type of field: 'nodal' or 'elemental'
%
% Outputs:
%
% visual:         data structure for visualisation
%
% See also: buildSubmeshPostprocess3D

function visual = interpolateSolutionPostprocess3D(mesh,refFace,postproc,visual,plotOpts,colorField,fieldType)

plotFaces = plotOpts.plotFacesCond;
nOfPlotFaces = size(plotFaces,1);

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
for j = 1:nOfPlotFaces
    iElem = plotFaces(j,1);
    iFace = plotFaces(j,2);
    p = mesh.pElem(iElem);
    
    Te = mesh.indexT(iElem,1):mesh.indexT(iElem,2);  
    
    iNodeEnd = iNodeIni + postproc.Face.nNodesPlot - 1;
    iNodePlotEnd = iNodePlotIni + refFace(p,p).nOfNodes - 1;    
    
    if strcmp(fieldType,'nodal')
        interpSol = postproc.Elem(p).face(iFace).N*colorField(Te,:);
    elseif strcmp(fieldType,'elemental')
        interpSol = repmat( colorField(iElem,:), postproc.Face.nNodesPlot, 1);
    end
        
    visual.U(iNodeIni:iNodeEnd,:) = interpSol;
    
    iNodeIni = iNodeIni + postproc.Face.nNodesPlot;
    iNodePlotIni = iNodePlotEnd + 1;  
end
