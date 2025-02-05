function res = test_zonotope_plot
% test_zonotope_plot - unit test function of plot; this function aims
%    to go through many variations of input arguments
%    note: only run-time errors checked, go through manually to check for bugs
%
% Syntax:
%    res = test_zonotope_plot
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
% See also: -

% Authors:       Mark Wetzlinger
% Written:       04-August-2020
% Last update:   09-May-2023 (TL, added plotted point checks)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% instantiate zonotope
Z = zonotope([1;-1;2],[2 -1 3; 0 1 -1; -1 4 2]);

try
    % try all variations in plotting
    figure;
    
    % one argument: object
    plot(Z);
    
    % two arguments: object, dimensions
    plot(Z,1);
    plot(Z,[1,2]);
    plot(Z,[2,3]);
    
    % three arguments: object, dimensions, linespec
    plot(Z,[1,2],'r+');
    
    % three arguments: object, dimensions, NVpairs
    plot(Z,[1,2],'LineWidth',2);
    plot(Z,[1,2],'Color',[.6 .6 .6],'LineWidth',2);
    plot(Z,[1,2],'EdgeColor','k','FaceColor',[.8 .8 .8]);
    
    % four arguments: object, dimensions, linespec, NVpairs
    plot(Z,[1,2],'r','LineWidth',2);
    plot(Z,[1,2],'r','LineWidth',2,'EdgeColor',[.6 .6 .6]);

    % check if plotted correctly
    ax = gca();
    colorOrder = ax.ColorOrder;
    
    % plot set
    plot(Z,[1,2]);
    V = [ ...
     3, -3, -5, -1, 5, 7, 3 ; ...
     -3, -1, 1, 1, -1, -3, -3 ; ...
    ];
    % check points
    assert(compareMatrices(V, [ax.Children(1).XData;ax.Children(1).YData],1e-4,'equal',true));
    % test color
    assert(isequal(colorOrder(1,:), ax.Children(1).Color));

    % check FaceAlpha
    plot(Z, 1:2, 'r', 'FaceAlpha', 0.1)
    assert(isequal([1 0 0], ax.Children(1).FaceColor))
    assert(isequal([1 0 0], ax.Children(1).EdgeColor))
    assert(isequal(0.1, ax.Children(1).FaceAlpha))
    assert(isequal(1, ax.Children(1).EdgeAlpha))
    plot(Z, 1:2, 'r', 'FaceColor', [0 1 0 0.1])

    % check EdgeAlpha
    plot(Z, 1:2, 'r', 'EdgeAlpha', 0.1)
    plot(Z, 1:2, 'EdgeColor','r', 'EdgeAlpha', 0.1,'FaceColor','b')
    plot(Z, 1:2, 'EdgeColor','r', 'EdgeAlpha', 0.1,'FaceColor','b','FaceAlpha',0.2)
    plot(Z, 1:3, 'r', 'EdgeAlpha', 0.1)
    plot(Z, 1:3, 'EdgeColor','r', 'EdgeAlpha', 0.1,'FaceColor','b')
    plot(Z, 1:3, 'EdgeColor','r', 'EdgeAlpha', 0.1,'FaceColor','b','FaceAlpha',0.2)
    plot(Z, 1:2, 'EdgeColor', [0 1 0 0.1])
    
    % close figure
    close;
catch ME
    close;
    rethrow(ME)
end

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
