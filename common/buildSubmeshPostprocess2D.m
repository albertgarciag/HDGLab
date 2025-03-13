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
% visual = buildSubmeshPostprocess2D(mesh, refElem, postproc)
%
% Inputs:
%
% mesh:             data structure for the mesh
% refElem:          data structure for the reference element
% postproc:          data structure for postprocess
%
% Outputs:
%
% visual:         data structure for visualisation
%
% See also: interpolateSolutionPostprocess2D

function visual = buildSubmeshPostprocess2D(mesh, refElem, postproc)

% Triangles 
nElementEdges = 3;

% Initialisation
nOfNodesSubMesh = 0;
nOfElementsSubMesh = 0;
nOfNodesSubMesh = nOfNodesSubMesh + postproc.nOfNodesPlot*mesh.nOfElements;
nOfElementsSubMesh = nOfElementsSubMesh + postproc.nSubElemsOnePlot*mesh.nOfElements;

visual.X = zeros(nOfNodesSubMesh, 3);
visual.T = zeros(nOfElementsSubMesh, 3);

nOfEdgeNodes = size(postproc.edgeNodesSplit,1);
visual.edges.X = zeros(nOfEdgeNodes,mesh.nOfElements*nElementEdges,3);

visual.nEdges = zeros(1,mesh.nOfElements);

nOfNodesPlot = size(mesh.X,1);
visual.Xnodes = zeros(nOfNodesPlot,3);

% Populate
iNodeIni = 1;
iElemIni = 1;
iNodePlotIni = 1;
indexEdge = 1;
for iElem = 1:mesh.nOfElements
    % Element info
    p = mesh.pElem(iElem);
    Te = mesh.indexT(iElem,1):mesh.indexT(iElem,2);
    Xe = mesh.X(Te,:);
    
    % Mesh nodes
    iNodePlotEnd = iNodePlotIni + refElem(p).nOfNodes - 1;    
    visual.Xnodes(iNodePlotIni:iNodePlotEnd,1:2) = Xe;
    iNodePlotIni = iNodePlotEnd + 1;
    
    % Submesh for plotting
    iNodeEnd = iNodeIni + postproc.nOfNodesPlot - 1;
    iElemEnd = iElemIni + postproc.nSubElemsOnePlot - 1;
    
    interpNodesElem = postproc.elem(p).N(:,:,1)*Xe;
    connecElem = postproc.connecNodesPlot + (iNodeIni-1);
    
    visual.X(iNodeIni:iNodeEnd,1:2) = interpNodesElem;
    visual.T(iElemIni:iElemEnd,:) = connecElem;
    
    iNodeIni = iNodeEnd + 1;
    iElemIni = iElemEnd + 1;
    
    % Edge info
    edgeNodesPlot = postproc.edgeNodesSplit;
    nOfEdges = postproc.nOfEdges;
    for iEdge = 1:nOfEdges
        visual.edges.X(:,indexEdge,1:2) = interpNodesElem(edgeNodesPlot(:,iEdge),1:2);
        indexEdge = indexEdge + 1;
    end
end

visual.edges.X = visual.edges.X(:,1:indexEdge-1,:);