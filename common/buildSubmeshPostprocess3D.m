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
% visual = buildSubmeshPostprocess3D(mesh,refElem,refFace,postproc,plotOpts)
%
% Inputs:
%
% mesh:             data structure for the mesh
% refElem:          data structure for the reference element
% refFace:          data structure for the reference face
% postproc:          data structure for postprocess
% plotOpts:         user options for visualisation
%
% Outputs:
%
% visual:         data structure for visualisation
%
% See also: interpolateSolutionPostprocess3D

function visual = buildSubmeshPostprocess3D(mesh,refElem,refFace,postproc,plotOpts)

plotFaces = plotOpts.plotFacesCond;
nOfPlotFaces = size(plotFaces,1);

% Initialisation
nOfNodesSubMesh = nOfPlotFaces*postproc.Face.nNodesPlot;
nOfElementsSubMesh = nOfPlotFaces*postproc.Face.nElemPlot;
nOfEdgeNodesPlot = length(postproc.Face.edgeNodesPlot);
visual.edges.X = zeros(nOfEdgeNodesPlot,nOfPlotFaces,3);

visual.X = zeros(nOfNodesSubMesh, 3);
visual.T = zeros(nOfElementsSubMesh, 3);

nOfNodesPlot = 0;
for j = 1:nOfPlotFaces
    iElem = plotFaces(j,1);
    p = mesh.pElem(iElem);
    nOfNodesPlot = nOfNodesPlot + refFace(p,p).nOfNodes;
end
visual.Xnodes = zeros(nOfNodesPlot,mesh.nsd);

% Populate
iNodeIni = 1;
iElemIni = 1;
iNodePlotIni = 1;
kEdge = 1;
for j = 1:nOfPlotFaces
    iElem = plotFaces(j,1);
    iFace = plotFaces(j,2);
    p = mesh.pElem(iElem);
    
    Te = mesh.indexT(iElem,1):mesh.indexT(iElem,2);
    Xe = mesh.X(Te,:);
    
    % Mesh nodes
    iNodePlotEnd = iNodePlotIni + refFace(p,p).nOfNodes - 1;
    visual.Xnodes(iNodePlotIni:iNodePlotEnd,:) = Xe(refElem(p).face(iFace).nodes,:);
    iNodePlotIni = iNodePlotEnd + 1;
    
    interpNodesElem = postproc.Elem(p).face(iFace).N*Xe;
    
    iNodeEnd = iNodeIni + postproc.Face.nNodesPlot - 1;
    iElemEnd = iElemIni + postproc.Face.nElemPlot - 1;
    
    visual.X(iNodeIni:iNodeEnd,:) = interpNodesElem;
    visual.T(iElemIni:iElemEnd,:) = postproc.Face.connecPlot + (iNodeIni-1);
    
    iNodeIni = iNodeIni + postproc.Face.nNodesPlot;
    iElemIni = iElemIni + postproc.Face.nElemPlot;
    
    visual.edges.X(:,kEdge,:) = interpNodesElem(postproc.Face.edgeNodesPlot,:);
    kEdge = kEdge + 1;
end