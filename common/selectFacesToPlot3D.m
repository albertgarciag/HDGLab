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
% plotFacesCond = selectFacesToPlot3D(mesh,postproc,plotFaces,conditionPlot)
%
% Inputs:
%
% mesh:             data structure for the mesh
% postproc:          data structure for postprocess
% plotFaces:        list of faces
% conditionPlot:    a condition to select a subset of faces
%
% Outputs:
%
% plotFacesCond:    subset of faces that satisfy the condition
%
% See also: buildSubmeshPostprocess3D, interpolateSolutionPostprocess3D

function plotFacesCond = selectFacesToPlot3D(mesh,postproc,plotFaces,conditionPlot)

if nargin<4
    conditionPlot = 'x<Inf';
end

% Count number of faces of each type
nOfInputFaces = size(plotFaces,1);
vectorPlotFaces = zeros(1,nOfInputFaces);

for j = 1:nOfInputFaces
    iElem = plotFaces(j,1);
    iFace = plotFaces(j,2);    
    Te = mesh.indexT(iElem,1):mesh.indexT(iElem,2);    
    faceVertices = postproc.faceVertices(iFace,:);
    Tfv = Te(faceVertices);
    Xfv = mesh.X(Tfv,:);
    x = Xfv(:,1);
    y = Xfv(:,2);
    z = Xfv(:,3);
    if all(eval(conditionPlot))
        vectorPlotFaces(j) = 1;
    end
end

plotFacesCond = plotFaces(vectorPlotFaces==1,:);