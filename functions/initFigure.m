function [fig,ax] = initFigure(fig)
%INITFIGURE Initialize figures
if nargin<1 ; fig = gcf ; end

fig = clf(fig,'reset') ;
    fig.NumberTitle = 'off' ;
    fig.Name = 'MASK CREATION: select ROI, substract background and apply threshold' ;
    fig.MenuBar = 'none' ; 
    fig.ToolBar = 'none' ;
ax = axes('Units','normalized','Outerposition',[0 0 1 1]) ;
    hold on ;
    axis(ax,'equal')
    axis(ax,'tight')
    axis(ax,'ij') ;
    
end

