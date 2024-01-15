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

resvec = [];

% instantiate polytope (via conversion from zonotope)
Z = zonotope([1;-1;2],[2 -1 3; 0 1 -1; -1 4 2]);
P = polytope(Z);

try
    % try all variations in plotting
    figure;
    
    % one argument: object
    plot(P);
    resvec(end+1) = true;
    
    % two arguments: object, dimensions
    plot(P,1);
    plot(P,[1,2]);
    plot(P,[2,3]);
    resvec(end+1) = true;
   
    % three arguments: object, dimensions, linespec
    plot(P,[1,2],'r+');
    resvec(end+1) = true;
    
    % three arguments: object, dimensions, NVpairs
    plot(P,[1,2],'LineWidth',2);
    plot(P,[1,2],'Color',[.6 .6 .6],'LineWidth',2);
    plot(P,[1,2],'EdgeColor','k','FaceColor',[.8 .8 .8]);
    resvec(end+1) = true;
    
    % four arguments: object, dimensions, linespec, NVpairs
    plot(P,[1,2],'r','LineWidth',2);
    plot(P,[1,2],'r','LineWidth',2,'EdgeColor',[.6 .6 .6]);
    resvec(end+1) = true;

    % check if plotted correctly
    ax = gca();
    colorOrder = ax.ColorOrder;
    
    % plot set
    plot(P,[1,2]);
    V = [-5 -3 3 7 5 -1 -5; 1 -1 -3 -3 -1 1 1];
    % check points
    resvec(end+1) = compareMatrices(V, [ax.Children(1).XData;ax.Children(1).YData],1e-4,'equal',true);
    % test color
    resvec(end+1) = isequal(colorOrder(1,:), ax.Children(1).Color);
    
    % close figure
    close;
catch ME
    close;
    resvec(end+1) = false;
end

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
