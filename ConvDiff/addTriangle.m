function m = addTriangle(x,y,color)
m = abs(diff(y(1:2))/diff(x(1:2)));

xlims=get(gca,'XLim');
xlimsize = xlims(2)-xlims(1);
latexFont = 'mwa_cmr10';
finalPoint = [x(1),y(1)];
tangentalDistance = xlimsize*0.1;

xWidth = diff(x(1:2));
Tri = nan(3,2);
lineVector = [1,m];
lineVector = lineVector./sqrt(sum(lineVector.^2));
lineVector = lineVector*(tangentalDistance/2);

normalVector = [m,-1];
normalVector = normalVector./sqrt(sum(normalVector.^2));
normalVector = normalVector*(tangentalDistance/2);

Tri(1,:) = finalPoint + normalVector + lineVector;
Tri(2,:) = Tri(1,:) + [xWidth/2,0];
Tri(3,:) = [Tri(2,1), m*Tri(2,1) + (Tri(1,2)-m*Tri(1,1))];

YAxisTextPos = [Tri(2,1),diff(Tri(2:3,2))*0.5+Tri(2,2)];
XAxisTextPos = [diff(Tri(1:2,1))*0.5+Tri(1,1),Tri(1,2)];

YTxt = text(YAxisTextPos(1),YAxisTextPos(2),sprintf('%2.1f',m));
XTxt = text(XAxisTextPos(1),XAxisTextPos(2),'1');

% move X text down by its own width * 0.6
XTxtExtent = get(XTxt,'Extent');
XTxtHeight = XTxtExtent(4);
set(XTxt,'Position',XAxisTextPos-[0,XTxtHeight*0.9]);

% move Y text across by its own width * 0.25 and down height * 0.5
YTxtExtent = get(YTxt,'Extent');
YTxtWidth = YTxtExtent(3);
set(YTxt,'VerticalAlignment','middle');
set(YTxt,'Position',YAxisTextPos+[YTxtWidth,0]*0.25);

triFontSize = 14;
set(YTxt,'Fontname','cmr12','Color',color,'FontSize',triFontSize,'Interpreter','latex');
set(XTxt,'Fontname','cmr12','Color',color,'FontSize',triFontSize,'Interpreter','latex');

hold on;
hGradPlot = plot(Tri([1:3,1],1),Tri([1:3,1],2),'lineWidth',1,'Color',color);
set(get(get(hGradPlot,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

end