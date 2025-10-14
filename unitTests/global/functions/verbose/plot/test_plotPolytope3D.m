function res = test_plotPolytope3D
% test_plotPolytope3D - unit test function for plotting of 3D polytopes
%
% Syntax:
%    res = test_plotPolytope3D
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
% Written:       18-April-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

I = interval([-1; -2; 1], [3; 1; 4]);
V = vertices(I);
tol = 1e-12;

try
    % test if plotting works ----------------------------------------------

    figure;

    % 3D polygon via interval and test

    plotPolytope3D(V);
    plotPolytope3D(V, 'r');
    plotPolytope3D(V, 'FaceColor', 'r');

    close

    % test if plotted correctly -------------------------------------------

    figure;
    ax = gca();

    colorOrder = ax.ColorOrder;

    % plot vertices
    plotPolytope3D(V);

    % check if all facets were plotted in 1 handle
    children = allchild(ax);
    assert(numel(children) == 1);

    % check correct plotting of facets 
    V_facets = [ ...
     -1, 3, 3, -1, -1, NaN, -1, 3, 3, -1, -1, NaN, -1, 3, 3, -1, -1, NaN, -1, 3, 3, -1, -1, NaN, 3, 3, 3, 3, 3, NaN, -1, -1, -1, -1, -1, NaN ; ...
     -2, -2, 1, 1, -2, NaN, -2, -2, 1, 1, -2, NaN, -2, -2, -2, -2, -2, NaN, 1, 1, 1, 1, 1, NaN, -2, 1, 1, -2, -2, NaN, -2, 1, 1, -2, -2, NaN ; ...
     1, 1, 1, 1, 1, NaN, 4, 4, 4, 4, 4, NaN, 1, 1, 4, 4, 1, NaN, 1, 1, 4, 4, 1, NaN, 1, 1, 4, 4, 1, NaN, 1, 1, 4, 4, 1, NaN ; ...
     ];
     % check points
    assert(all(withinTol(V_facets, [ax.Children(1).XData;ax.Children(1).YData;ax.Children(1).ZData]) | isnan(V_facets),"all"));
    
    % test color
    assert(isequal(colorOrder(1, :), ax.Children(1).Color));

    % replot and check if previous plot was deleted and color order is ok
    plotPolytope3D(V);
    children = allchild(ax);
    assert(numel(children) == 1);
    assert(isequal(colorOrder(1, :), ax.Children(1).Color));

    % check specified color
    plotPolytope3D(V, 'r');
    % test color
    if CORA_PLOT_FILLED
        assert(isequal([1, 0, 0], ax.Children(1).EdgeColor));
        assert(isequal([1, 0, 0], ax.Children(1).FaceColor));
    else
        assert(isequal([1, 0, 0], ax.Children(1).Color));
    end

    % check face color
    plotPolytope3D(V, 'FaceColor','r','FaceAlpha',0.1);
    V_facets = { [ -1 3 3 -1 ; -2 -2 1 1 ; 1 1 1 1 ] [ -1 3 3 -1 ; -2 -2 1 1 ; 4 4 4 4 ] [ -1 3 3 -1 ; -2 -2 -2 -2 ; 1 1 4 4 ] [ -1 3 3 -1 ; 1 1 1 1 ; 1 1 4 4 ] [ 3 3 3 3 ; -2 1 1 -2 ; 1 1 4 4 ] [ -1 -1 -1 -1 ; -2 1 1 -2 ; 1 1 4 4 ] };
    assert(isequal([1, 0, 0], ax.Children(1).FaceColor));
    % check faces
    children = allchild(ax);
    V_facets_plotted = arrayfun(@(facet) [facet.XData'; facet.YData'; facet.ZData'], children(end:-1:1), 'UniformOutput', false)';
    % when plotting with face color, it's not important to close the regions
    assert(all(arrayfun(@(i) compareMatrices(V_facets{i},V_facets_plotted{i},tol,'subset'), 1:numel(V_facets))))
    close;

catch ME
    close
    rethrow(ME)
end

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
