%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Postprocessing for the Stokes equation using isoparametric elements   
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

disp('Postprocess...')
plotOpts.resolution = 2;            % 1: low resolution, 2: high resolution
plotOpts.fieldsWithMesh = 1;        % 0: false, 1: true
plotOpts.fieldsWithNodes = 0;       % 0: false, 1: true   
plotOpts.alphaFace = 1;             % Transparency (0 to 1)
plotOpts.componentsU = 0;           % 0: false, 1: true   
plotOpts.moduleU = 1;               % 0: false, 1: true   
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
    visualStar = buildSubmeshPostprocess2D(meshStar, refElem, postproc);
    
    % Plot mesh
    plotOpts.componentsToPlot = 1;
    visual = interpolateSolutionPostprocess2D(mesh, refElem, postproc, visual, [], 'elemental');
    iFig = postprocessField2D(visual, plotOpts, iFig);
        
    visual = interpolateSolutionPostprocess2D(mesh, refElem, postproc, visual, u, 'nodal');
    visualStar = interpolateSolutionPostprocess2D(meshStar, refElem, postproc, visualStar, uStar, 'nodal');
    
    if plotOpts.componentsU
        plotOpts.componentsToPlot = 1:ctt.nOfComponents;   % vector of components to plot
        % Plot components of u    
        iFig = postprocessField2D(visual, plotOpts, iFig);
        % Plot components of uStar
        iFig = postprocessField2D(visualStar, plotOpts, iFig);
    end

    if plotOpts.moduleU
        modU = sqrt(sum(visual.U.^2,2));
        modUStar = sqrt(sum(visualStar.U.^2,2));
        
        plotOpts.componentsToPlot = 1;
        % Plot module of u
        visual.U = modU;
        iFig = postprocessField2D(visual, plotOpts, iFig);        
        % Plot module of uStar
        visualStar.U = modUStar;
        iFig = postprocessField2D(visualStar, plotOpts, iFig);
    end
    
    % Plot p
    plotOpts.componentsToPlot = 1;
    visual = interpolateSolutionPostprocess2D(mesh, refElem, postproc, visual, p, 'nodal');
    iFig = postprocessField2D(visual, plotOpts, iFig);
    
elseif mesh.nsd==3
    plotOpts.plotFaces = mesh.extFaces(:, 1:2);                      
    conditionPlot = 'x<100';
    
    plotOpts.plotFacesCond = selectFacesToPlot3D(mesh, postproc, plotOpts.plotFaces, conditionPlot);
    visual = buildSubmeshPostprocess3D(mesh, refElem, refFace, postproc, plotOpts);
    visualStar = buildSubmeshPostprocess3D(meshStar, refElem, refFace, postproc, plotOpts);
    
    % Plot mesh
    plotOpts.componentsToPlot = 1;
    visual = interpolateSolutionPostprocess3D(mesh, refFace, postproc, visual, plotOpts, [], 'nodal');
    iFig = postprocessField3D(visual, plotOpts, iFig);
    
    visual = interpolateSolutionPostprocess3D(mesh, refFace, postproc, visual, plotOpts, u, 'nodal');
    visualStar = interpolateSolutionPostprocess3D(meshStar, refFace, postproc, visualStar, plotOpts, uStar, 'nodal');
            
    if plotOpts.componentsU
        plotOpts.componentsToPlot = 1:ctt.nOfComponents;   % vector of components to plot
        % Plot components of u
        iFig = postprocessField3D(visual, plotOpts, iFig);        
        % Plot components of uStar
        iFig = postprocessField3D(visualStar, plotOpts, iFig);
    end

    if plotOpts.moduleU
        modU = sqrt(sum(visual.U.^2,2));
        modUStar = sqrt(sum(visualStar.U.^2,2));
        
        plotOpts.componentsToPlot = 1;
        % Plot module of u
        visual.U = modU;
        iFig = postprocessField3D(visual, plotOpts, iFig);
        % Plot module of uStar
        visualStar.U = modUStar;
        iFig = postprocessField3D(visualStar, plotOpts, iFig);
    end  
    
    % Plot p
    plotOpts.componentsToPlot = 1;
    visual = interpolateSolutionPostprocess3D(mesh, refFace, postproc, visual, plotOpts, p, 'nodal');
    iFig = postprocessField3D(visual, plotOpts, iFig);    
end