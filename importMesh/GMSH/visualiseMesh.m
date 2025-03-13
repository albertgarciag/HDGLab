%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Visualise the imported mesh with the node location
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

disp('Visualising...')
plotOpts.resolution = 2;            % 1: low resolution, 2: high resolution
plotOpts.fieldsWithMesh = 1;        % 0: false, 1: true
plotOpts.fieldsWithNodes = 1;       % 0: false, 1: true   
plotOpts.alphaFace = 1;             % Transparency (0 to 1)
plotOpts.componentsU = 1;           % 0: false, 1: true   
plotOpts.moduleU = 1;               % 0: false, 1: true   
if ~exist('iFig','var')
    iFig = 0;
end

[refElemPlot, refFacePlot] = getRefData(pDegree, nsd, optionNodes);

if mesh.nsd==2
    load(sprintf('postprocessTRI_res%d', plotOpts.resolution))
    
    visual = buildSubmeshPostprocess2D(mesh, refElemPlot, postproc);
    
    % Plot mesh
    plotOpts.componentsToPlot = 1;
    visual = interpolateSolutionPostprocess2D(mesh, refElemPlot, postproc, visual, [], 'elemental');
    iFig = postprocessField2D(visual, plotOpts, iFig);
elseif mesh.nsd==3
    load(sprintf('postprocessTET_res%d', plotOpts.resolution))
    
    plotOpts.plotFaces = mesh.extFaces(:, 1:2);
    conditionPlot = 'x<100';
    
    plotOpts.plotFacesCond = selectFacesToPlot3D(mesh, postproc, plotOpts.plotFaces, conditionPlot);
    visual = buildSubmeshPostprocess3D(mesh, refElemPlot, refFacePlot, postproc, plotOpts);
    
    % Plot mesh
    plotOpts.componentsToPlot = 1;
    visual = interpolateSolutionPostprocess3D(mesh, refFacePlot, postproc, visual, plotOpts, [], 'nodal');
    iFig = postprocessField3D(visual, plotOpts, iFig);
end