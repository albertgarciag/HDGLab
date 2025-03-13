%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Postprocessing for the Poisson equation using isoparametric elements   
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
disp('Postprocess...')
plotOpts.resolution = 2;            % 1: low resolution, 2: high resolution
plotOpts.fieldsWithMesh = 1;        % 0: false, 1: true
plotOpts.fieldsWithNodes = 0;       % 0: false, 1: true
plotOpts.componentsToPlot = 1;      % vector of components to plot
plotOpts.alphaFace = 1;             % Transparency (0 to 1)
if ~exist('iFig','var')
    iFig = 0;
end

if mesh.nsd==2
    load(sprintf('postprocessTRI_opt%d_res%d', mesh.optionNodes, plotOpts.resolution))
elseif mesh.nsd==3
    load(sprintf('postprocessTET_opt%d_res%d', mesh.optionNodes, plotOpts.resolution))
end

meshStar = buildIsoparametricHOmesh(mesh, refElem);

if mesh.nsd==2
    visual = buildSubmeshPostprocess2D(mesh, refElem, postproc);
    visualStar = buildSubmeshPostprocess2D(meshStar,refElem,postproc);
    
    % Plot mesh
    visual = interpolateSolutionPostprocess2D(mesh,refElem,postproc,visual,[],'elemental');
    iFig = postprocessField2D(visual,plotOpts,iFig);
    
    % Plot u
    visual = interpolateSolutionPostprocess2D(mesh,refElem,postproc,visual,u,'nodal');
    iFig = postprocessField2D(visual,plotOpts,iFig);
    
    % Plot uStar
    visualStar = interpolateSolutionPostprocess2D(meshStar,refElem,postproc,visualStar,uStar,'nodal');
    iFig = postprocessField2D(visualStar,plotOpts,iFig);
    
elseif mesh.nsd==3
    plotOpts.plotFaces = mesh.extFaces(:,1:2);
     conditionPlot = 'x<100';
    
    plotOpts.plotFacesCond = selectFacesToPlot3D(mesh,postproc,plotOpts.plotFaces,conditionPlot);
    visual = buildSubmeshPostprocess3D(mesh,refElem,refFace,postproc,plotOpts);
    visualStar = buildSubmeshPostprocess3D(meshStar,refElem,refFace,postproc,plotOpts);
    
    % Plot mesh
    visual = interpolateSolutionPostprocess3D(mesh,refFace,postproc,visual,plotOpts,[], 'nodal');
    iFig = postprocessField3D(visual,plotOpts,iFig);
    
    % Plot u
    visual = interpolateSolutionPostprocess3D(mesh,refFace,postproc,visual,plotOpts,u, 'nodal');
    iFig = postprocessField3D(visual,plotOpts,iFig);
    
    % Plot uStar
    visualStar = interpolateSolutionPostprocess3D(meshStar,refFace,postproc,visualStar,plotOpts,uStar, 'nodal');
    iFig = postprocessField3D(visualStar,plotOpts,iFig);
end