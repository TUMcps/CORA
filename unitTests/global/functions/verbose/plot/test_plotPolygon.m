function res = test_plotPolygon
% test_plotPolygon - unit test function for plotting of 2D/3D polygons
%
% Syntax:
%    res = test_plotPolygon
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger, Tobias Ladner
% Written:       30-April-2023
% Last update:   05-May-2023 (TL, check if plotted correctly)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

try
    % test if plotting works ----------------------------------------------

    figure;

    % 2D polygon
    V = [1 0; 1 2; 0 3; -2 2; -3 0; 0 -1; 1 0]';
    plotPolygon(V);
    plotPolygon(V,'r');
    plotPolygon(V,'FaceColor','r');

    resvec(end+1) = true;

    % point cloud with convex hull
    V = randn(2,100);
    plotPolygon(V,'convHull',true);

    % plot at position
    V = randn(2,100);
    plotPolygon(V,'XPos',2);
    plotPolygon(V,'YPos',2);
    plotPolygon(V,'ZPos',2);

    resvec(end+1) = true;

    % 3D polygon
    V = [1 0 1; 1 2 -1; 0 3 0; -2 2 1; -3 0 2; 0 -1 0; 1 0 1]';
    plotPolygon(V);
    plotPolygon(V,'r');
    plotPolygon(V,'FaceColor','r');

    resvec(end+1) = true;
    close;

    % test if plotted correctly -------------------------------------------

    figure;
    ax = gca();

    colorOrder = ax.ColorOrder;

    V = [1 3 5 4 1; -2 4 5 4 -2];
    plotPolygon(V);

    % test points
    resvec(end+1) = isequal(V, [ax.Children(1).XData;ax.Children(1).YData]);
    % test color
    resvec(end+1) = isequal(colorOrder(1,:), ax.Children(1).Color);

    % replot
    plotPolygon(V);
    % first color as hold is off
    resvec(end+1) = isequal(colorOrder(1,:), ax.Children(1).Color);
    
    hold on;
    ax = gca();
    % now choose next colors (axes are in reverse order)
    resvec(end+1) = isequal(colorOrder(1,:), ax.Children(1).Color);
    plotPolygon(V);
    resvec(end+1) = isequal(colorOrder(2,:), ax.Children(1).Color);
    plotPolygon(V);
    resvec(end+1) = isequal(colorOrder(3,:), ax.Children(1).Color);
    % check whether first plot is still there
    resvec(end+1) = isequal(colorOrder(1,:), ax.Children(3).Color);

    close

    figure; ax = gca();
    
    % test convHull
    plotPolygon(V, 'ConvHull',true);
    resvec(end+1) = isequal([1 5 3 1; -2 5 4 -2], [ax.Children(1).XData;ax.Children(1).YData]);

    % test XPos,YPos,ZPos
    plotPolygon(V(1,:),'XPos',2);
    resvec(end+1) = isequal([2 2 2 2 2; V(1,:)], [ax.Children(1).XData;ax.Children(1).YData]);
    plotPolygon(V(1,:),'YPos',3);
    resvec(end+1) = isequal([V(1,:); 3 3 3 3 3], [ax.Children(1).XData;ax.Children(1).YData]);
    plotPolygon(V,'ZPos',4);
    resvec(end+1) = isequal([V; 4 4 4 4 4], [ax.Children(1).XData;ax.Children(1).YData;ax.Children(1).ZData]);
    % no y given
    plotPolygon(V(1,:),'ZPos',1);
    resvec(end+1) = isequal([V(1,:); 0 0 0 0 0; 1 1 1 1 1], [ax.Children(1).XData;ax.Children(1).YData;ax.Children(1).ZData]);

    close;

catch ME
    close
    resvec(end+1) = false;
end

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
