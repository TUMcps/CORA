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

    % check if all facets were plotted
    children = allchild(ax);
    assert(numel(children) == 6 + 1);% +1 due to nan

    % check correct plotting of facets 
    i = 1; % facet 1
    V_facet = [ ...
        -1.000, -1.000, -1.000, -1.000, -1.000; ...
        -2.000, 1.000, 1.000, -2.000, -2.000; ...
        1.000, 1.000, 4.000, 4.000, 1.000; ...
        ];
    assert(compareMatrices(V_facet, [children(i).XData; children(i).YData; children(i).ZData], 1e-4, 'equal', true));
    i = 2; % facet 2
    V_facet = [ ...
        3.000, 3.000, 3.000, 3.000, 3.000; ...
        -2.000, 1.000, 1.000, -2.000, -2.000; ...
        1.000, 1.000, 4.000, 4.000, 1.000; ...
        ];
    assert(compareMatrices(V_facet, [children(i).XData; children(i).YData; children(i).ZData], 1e-4, 'equal', true));
    i = 3; % facet 3
    V_facet = [ ...
        -1.000, 3.000, 3.000, -1.000, -1.000; ...
        1.000, 1.000, 1.000, 1.000, 1.000; ...
        1.000, 1.000, 4.000, 4.000, 1.000; ...
        ];
    assert(compareMatrices(V_facet, [children(i).XData; children(i).YData; children(i).ZData], 1e-4, 'equal', true));
    i = 4; % facet 4
    V_facet = [ ...
        -1.000, 3.000, 3.000, -1.000, -1.000; ...
        -2.000, -2.000, -2.000, -2.000, -2.000; ...
        1.000, 1.000, 4.000, 4.000, 1.000; ...
        ];
    assert(compareMatrices(V_facet, [children(i).XData; children(i).YData; children(i).ZData], 1e-4, 'equal', true));
    i = 5; % facet 4
    V_facet = [ ...
        -1.000, 3.000, 3.000, -1.000, -1.000; ...
        -2.000, -2.000, 1.000, 1.000, -2.000; ...
        4.000, 4.000, 4.000, 4.000, 4.000; ...
        ];
    assert(compareMatrices(V_facet, [children(i).XData; children(i).YData; children(i).ZData], 1e-4, 'equal', true));
    i = 6; % facet 6
    V_facet = [ ...
        -1.000, 3.000, 3.000, -1.000, -1.000; ...
        -2.000, -2.000, 1.000, 1.000, -2.000; ...
        1.000, 1.000, 1.000, 1.000, 1.000; ...
        ];
    assert(compareMatrices(V_facet, [children(i).XData; children(i).YData; children(i).ZData], 1e-4, 'equal', true));
    i = 7; % last child should be nan, used to delete previous plots as hold off
    assert(all(isnan([children(i).XData; children(i).YData; children(i).ZData])));

    % test color
    assert(isequal(colorOrder(1, :), ax.Children(1).Color));

    % replot and check if previous plot was deleted and color order is ok
    plotPolytope3D(V);
    children = allchild(ax);
    assert(numel(children) == 6 + 1);% +1 due to nan
    assert(isequal(colorOrder(1, :), ax.Children(1).Color));

    % check specified color
    plotPolytope3D(V, 'r');
    assert(isequal([1, 0, 0], ax.Children(1).Color));

    close;

catch ME
    close
    rethrow(ME)
end

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
