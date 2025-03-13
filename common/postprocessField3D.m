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
% iFig = postprocessField3D(visual,plotOpts,iFig)
%
% Inputs:
%
% visual:         data structure for visualisation
% plotOpts:         user options for visualisation
% iFig:             figure number
%
% Outputs:
%
% iFig:             figure number plus one
%
% See also: buildSubmeshPostprocess2D, interpolateSolutionPostprocess2D
%           trisurf, plot3, colormap, colorbar

function iFig = postprocessField3D(visual,plotOpts,iFig)

nOfComponents = numel(plotOpts.componentsToPlot);

for kComp = 1:nOfComponents
    iComp = plotOpts.componentsToPlot(kComp);
    
    iFig = iFig + 1;
    figure(iFig)
    
    trisurf(visual.T, visual.X(:,1), visual.X(:,2), visual.X(:,3), visual.U(:,iComp), ...
        'EdgeColor','none','FaceColor','interp','FaceLighting','flat','FaceAlpha',plotOpts.alphaFace)
    hold on
    
    if plotOpts.fieldsWithMesh
        plot3(visual.edges.X(:,:,1), visual.edges.X(:,:,2), visual.edges.X(:,:,3), 'k')
    end
    if plotOpts.fieldsWithNodes
        plot3(visual.Xnodes(:,1), visual.Xnodes(:,2), visual.Xnodes(:,3), 'ko', 'MarkerSize', 2, 'MarkerFaceColor', 'k')
    end
    view(3)
    axis equal
    axis off
    colormap jet
    hcb=colorbar;
    set(hcb,'FontSize', 18,'FontName','mwa_cmr10')
end