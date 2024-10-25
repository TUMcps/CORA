function res = test_polytope_plot
% test_polytope_plot - unit test function of plot
%
% Syntax:
%    res = test_polytope_plot
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

% Authors:       Mark Wetzlinger
% Written:       30-November-2022
% Last update:   09-May-2023 (TL, added plotted point checks)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% instantiate polytope (via conversion from zonotope)
Z = zonotope([1;-1;2],[2 -1 3; 0 1 -1; -1 4 2]);
P = polytope(Z);

try
    % try all variations in plotting
    figure;
    
    % one argument: object
    plot(P);
    
    % two arguments: object, dimensions
    plot(P,1);
    plot(P,[1,2]);
    plot(P,[2,3]);
   
    % three arguments: object, dimensions, linespec
    plot(P,[1,2],'r+');
    
    % three arguments: object, dimensions, NVpairs
    plot(P,[1,2],'LineWidth',2);
    plot(P,[1,2],'Color',[.6 .6 .6],'LineWidth',2);
    plot(P,[1,2],'EdgeColor','k','FaceColor',[.8 .8 .8]);
    
    % four arguments: object, dimensions, linespec, NVpairs
    plot(P,[1,2],'r','LineWidth',2);
    plot(P,[1,2],'r','LineWidth',2,'EdgeColor',[.6 .6 .6]);

    % check if plotted correctly
    ax = gca();
    colorOrder = ax.ColorOrder;
    
    % plot set
    plot(P,[1,2]);
    V = [ ...
     -5, -1, 5, 7, 3, -3, -5 ; ...
     1, 1, -1, -3, -3, -1, 1 ; ...
    ];
    % check points
    assert(compareMatrices(V, [ax.Children(1).XData;ax.Children(1).YData],1e-4,'equal',true));
    % test color
    assert(isequal(colorOrder(1,:), ax.Children(1).Color));

    % check polytope with redundant vertices
    V = [ 1.000 4.000 4.000 1.000 1.000 4.000 4.000 4.000 4.000 7.000 7.000 4.000 4.000 7.000 7.000 1.000 ; 3.000 3.000 6.000 3.000 3.000 3.000 6.000 0.000 0.000 0.000 3.000 2.000 2.000 2.000 5.000 3.000 ];
    P = polytope(V);
    plot(P)

    V_true = [ ...
     1, 4, 7, 7, 4, 1 ; ...
     3, 6, 5, 0, 0, 3 ; ...
    ];
    % check points
    assert(compareMatrices(V_true, [ax.Children(1).XData;ax.Children(1).YData],1e-4,'equal',true));
    assert(all(contains(P,V)))
    
    % close figure
    close;
catch ME
    close;
    rethrow(ME)
end

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
