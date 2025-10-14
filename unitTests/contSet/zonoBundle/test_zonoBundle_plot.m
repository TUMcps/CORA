function res = test_zonoBundle_plot
% test_zonoBundle_plot - unit test function of plot; this function aims
%    to go through many variations of input arguments
%    note: only run-time errors checked, go through manually to check for bugs
%
% Syntax:
%    res = test_zonoBundle_plot
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
% Written:       25-May-2022
% Last update:   09-May-2023 (TL, added plotted point checks)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% instantiate zonotope bundle
Z1 = zonotope([1;-1;2],[2 -1 3; 0 1 -1; -1 4 2]);
Z2 = Z1 + [1;0;0];
zB = zonoBundle({Z1,Z2});

try
    % try all variations in plotting
    figure;
    
    % one argument: object
    plot(zB);

    % two arguments: object, dimensions
    plot(zB,1);
    plot(zB,[2,3]);

    % three arguments: object, dimensions, linespec
    plot(zB,[1,2],'r+');
    
    % three arguments: object, dimensions, NVpairs
    plot(zB,[1,2],'LineWidth',2);
    plot(zB,[1,2],'Color',[.6 .6 .6],'LineWidth',2);
    plot(zB,[1,2],'EdgeColor','k','FaceColor',[.8 .8 .8]);
    
    % four arguments: object, dimensions, linespec, NVpairs
    plot(zB,[1,2],'r','LineWidth',2);
    plot(zB,[1,2],'r','LineWidth',2,'EdgeColor',[.6 .6 .6]);

    % check if plotted correctly
    ax = gca();
    colorOrder = ax.ColorOrder;
    
    % plot set
    plot(zB,[1,2]);
    V = [ ...
     -1, 5, 7, 4, -2, -4, -1 ; ...
     1, -1, -3, -3, -1, 1, 1 ; ...
    ];
    % check points
    assert(compareMatrices(V, readVerticesFromFigure(ax.Children(1)),1e-4,'equal',true));
    % test color
    if CORA_PLOT_FILLED
        assert(isequal(colorOrder(1,:), ax.Children(1).EdgeColor));
        assert(isequal(colorOrder(1,:), ax.Children(1).FaceColor));
    else
        assert(isequal(colorOrder(1,:), ax.Children(1).Color));
    end

    % check barely intersecting sets
    Z{1} = zonotope([1;1],[1 0; 0 1]);
    Z{2} = zonotope([2;2],[1;0]);
    zB = zonoBundle(Z);
    plot(zB)

    V = [
        2 1 2
        2 2 2
    ];
    assert(compareMatrices(V, readVerticesFromFigure(ax.Children(1)),1e-4,'equal',true));
        
    % close figure
    close;
catch ME
    close;
    rethrow(ME)
end

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
