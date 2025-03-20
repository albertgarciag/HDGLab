clearvars
close all

addpath ./export_fig

addpath ./resConvDiff/convergence/
% addpath ../resStokes(convergenceGMSH

%% Setup
domain = 'square';
% domain = 'cube';
% rescale = 1;
% domain = 'ring';
rescale = 5;
nsd = 2;
exID = 3; % 3 only example implemented

isExportFig = 1;

% hMeshes = [4 4 4 4 3];
hMeshes = [4 4 4 4 4];
% hMeshes = [3 3 3 3 3];
pRefinements = [1 2 3 4 5];

col = {'b', 'r', 'g', 'm', 'k'};
sty = {'-', '--'};
sym = {'o', 's', 'd', '^', '*'};

%% PLOT u, u^*
figure(2), hold on

for iDeg = 1:numel(pRefinements)
    P = pRefinements(iDeg);
    nOfMeshes = hMeshes(iDeg);
    
    h = zeros(nOfMeshes,1);
    errU = zeros(nOfMeshes,1);
    errUstar = zeros(nOfMeshes,1);

    for iMesh = 1:nOfMeshes
        H = iMesh;
        
        solutionFile = sprintf('%sH%dP%d_ex%d.mat', domain, H, P, exID);
        load(solutionFile,'mesh','uErr')
        
        h(iMesh) = computeElementDiameter(mesh)/rescale;
        errU(iMesh) = sqrt(sum((uErr.L2).^2,2));
        errUstar(iMesh) = sqrt(sum((uErr.L2Star).^2,2));
    end
    plot(log10(h), log10(errU), sprintf('%s%s',col{iDeg},sty{1},sym{iDeg}), 'LineWidth', 2, 'MarkerSize', 8)
    plot(log10(h), log10(errUstar), sprintf('%s%s',col{iDeg},sty{2},sym{iDeg}), 'LineWidth', 2, 'MarkerSize', 8)        
    %
    xx = [log10(h(end)) log10(h(end-1))];
    yy = [log10(errU(end)) log10(errU(end-1))];
    addTriangle(xx,yy,col{iDeg});
    %
    if iDeg==numel(pRefinements)
        xx1 = [log10(h(end)) log10(h(end-1))];
        yy1 = [log10(errUstar(end)) log10(errUstar(end-1))];
        addTriangle(xx1,yy1,col{iDeg});        
    end
end

box on
grid on
xlabel('log$_{10}$($h$)','Interpreter','latex','FontName','cmr12')
ylabel('log$_{10}(||E||_{L^2(\Omega)})$','Interpreter','latex','FontName','cmr12')

if nsd == 2
    switch exID
        case 1
            % not implemented
            set(gca,'XLim',[-1.58 -0.3],'XTick',[-1.5:.2:-0.3],'YLim',[-12 0],'YTick',[-12:2:0],'FontName','cmr12')
        case 2
            % not implemented
            set(gca,'XLim',[-1.56 -0.3],'XTick',[-1.5:.2:-0.3],'YLim',[-10 2],'YTick',[-10:2:2],'FontName','cmr12')
        case 3
            set(gca,'XLim',[-2 -0.5],'XTick',[-1.9:.3:-0.8],'YLim',[-12 0],'YTick',[-12:2:0],'FontName','cmr12')
    end
else
    % not implemented
    set(gca,'XLim',[-0.85 0],'XTick',[-0.8:.2:0],'YLim',[-12 0],'YTick',[-12:2:0],'FontName','cmr12')
end

set(gca,'FontSize',22,'FontName','cmr12','TickLabelInterpreter','latex');
leg = {sprintf('{\\boldmath{$u$}}${, {\\texttt{p}}{=}%d}$', pRefinements(1)), sprintf('{\\boldmath{$u_{\\star}$}}${, {\\texttt{p}}{=}%d}$', pRefinements(1)), ...
       sprintf('{\\boldmath{$u$}}${, {\\texttt{p}}{=}%d}$', pRefinements(2)), sprintf('{\\boldmath{$u_{\\star}$}}${, {\\texttt{p}}{=}%d}$', pRefinements(2)), ...
       sprintf('{\\boldmath{$u$}}${, {\\texttt{p}}{=}%d}$', pRefinements(3)), sprintf('{\\boldmath{$u_{\\star}$}}${, {\\texttt{p}}{=}%d}$', pRefinements(3)), ...
       sprintf('{\\boldmath{$u$}}${, {\\texttt{p}}{=}%d}$', pRefinements(4)), sprintf('{\\boldmath{$u_{\\star}$}}${, {\\texttt{p}}{=}%d}$', pRefinements(4)), ...
       sprintf('{\\boldmath{$u$}}${, {\\texttt{p}}{=}%d}$', pRefinements(5)), sprintf('{\\boldmath{$u_{\\star}$}}${, {\\texttt{p}}{=}%d}$', pRefinements(5))
       };
hLeg = legend(leg, 'Location', 'northeast', 'FontName','cmr12');
set(hLeg,'Interpreter','latex');

%% Export
set(gcf,'Color','w');
if isExportFig
    if nsd == 2
        switch exID
            case 1
                export_fig convdiffD_hConvUUstar_ex1 -pdf -transparent -r800
            case 2
                export_fig convdiff_hConvUUstar_ex2 -pdf -transparent -r800
            case 3
                export_fig convdiff2D_hConvUUstar_ex3 -pdf -transparent -r800
        end
    else
        export_fig stokes3D_hConvUUstar_ex1 -pdf -transparent -r800
    end
end
hold off