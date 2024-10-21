function res = test_plotPoints
% test_plotPoints - unit test function for plotting points in 2D/3D
%
% Syntax:
%    res = test_plotPoints
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

% Authors:       Tobias Ladner
% Written:       08-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

try
    % test if plotting works ----------------------------------------------

    figure;

    % 2D
    V = [1 0; 1 2; 0 3; -2 2; -3 0; 0 -1; 1 0]';
    plotPoints(V);
    plotPoints(V,1:2,'.r');
    plotPoints(V,1:2,'FaceColor','r');

    % point cloud with convex hull
    V = randn(2,100);
    plotPoints(V,1:2,'ConvHull',true);

    % plot at position
    V = randn(2,100);
    plotPoints(V,1:2,'XPos',2);
    plotPoints(V,1:2,'YPos',2);
    plotPoints(V,1:2,'ZPos',2);

    % 3D polygon
    V = [1 0 1; 1 2 -1; 0 3 0; -2 2 1; -3 0 2; 0 -1 0; 1 0 1]';
    plotPoints(V,1:2);
    plotPoints(V,1:2,'.r');
    plotPoints(V,1:2,'FaceColor','r');

    close;

    % test if plotted correctly -------------------------------------------

    figure;
    ax = gca();

    colorOrder = ax.ColorOrder;

    V = [1 3 5 4 1; -2 4 5 4 -2];
    plotPoints(V);

    % test points
    assert(isequal(V, [ax.Children(1).XData;ax.Children(1).YData]));
    % test color
    assert(isequal(colorOrder(1,:), ax.Children(1).Color));

    % replot
    plotPoints(V);
    % first color as hold is off
    assert(isequal(colorOrder(1,:), ax.Children(1).Color));
    
    hold on;
    ax = gca();
    % now choose next colors (axes are in reverse order)
    assert(isequal(colorOrder(1,:), ax.Children(1).Color));
    plotPoints(V);
    assert(isequal(colorOrder(2,:), ax.Children(1).Color));
    plotPoints(V);
    assert(isequal(colorOrder(3,:), ax.Children(1).Color));
    % check whether first plot is still there
    assert(isequal(colorOrder(1,:), ax.Children(3).Color));

    close

    figure; ax = gca();
    
    % test convHull
    plotPoints(V,1:2,'ConvHull',true);
    assert(isequal([1 5 3 1; -2 5 4 -2], [ax.Children(1).XData;ax.Children(1).YData]));

    % test XPos,YPos,ZPos
    plotPoints(V,1,'XPos',2);
    assert(isequal([2 2 2 2 2; V(1,:)], [ax.Children(1).XData;ax.Children(1).YData]));
    plotPoints(V,1,'YPos',3);
    assert(isequal([V(1,:); 3 3 3 3 3], [ax.Children(1).XData;ax.Children(1).YData]));
    plotPoints(V,1:2,'ZPos',4);
    assert(isequal([V; 4 4 4 4 4], [ax.Children(1).XData;ax.Children(1).YData;ax.Children(1).ZData]));
    % no y given
    plotPoints(V,1,'ZPos',1);
    assert(isequal([V(1,:); 0 0 0 0 0; 1 1 1 1 1], [ax.Children(1).XData;ax.Children(1).YData;ax.Children(1).ZData]));

    close;

catch ME
    close
    rethrow(ME)
end

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
