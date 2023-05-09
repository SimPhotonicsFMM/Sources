function PlotCourbe(x,y,LabelX,LabelY,Title,PropLine,Legend)

% PlotCourbe
%   Linear plot of 1D curves as function "plot" with some properties
%
% Syntax
%   PlotCourbe(x,y,LabelX,LabelY,Title,PropLine,Legend);
%
% Description
%   x : Abscissa array
%   y : Ordinate array
%   LabelX : label in x
%   LabelY : label in y
%   Title  : Title of the figure
%   PropLine : curve properties (symbols, color)
%   Legend : Legend

% Date of the latest version : 10 February 2023
% Author : Mondher Besbes (LCF / CNRS / IOGS)

plot(x,y,PropLine,'LineWidth',2)   

axis tight
xlabel(LabelX,'FontWeight','bold','FontSize',12)
ylabel(LabelY,'FontWeight','bold','FontSize',12)
title(Title,'FontWeight','bold','FontSize',12)
%grid on
box on
set(gca,'LineWidth',2,'FontWeight','bold','FontSize',12)

if nargin == 7
    legend1 = legend(Legend);
    set(legend1,...
        'FontSize',11,...
        'LineWidth',1.5,...
        'FontWeight','bold');
end

set(gcf,'Color','w')

return
