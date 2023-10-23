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

resvec = [];

% instantiate zonotope bundle
Z1 = zonotope([1;-1;2],[2 -1 3; 0 1 -1; -1 4 2]);
Z2 = Z1 + [1;0;0];
zB = zonoBundle({Z1,Z2});

try
    % try all variations in plotting
    figure;
    
    % one argument: object
    plot(zB);
    resvec(end+1) = true;

    % two arguments: object, dimensions
    plot(zB,1);
    plot(zB,[2,3]);
    resvec(end+1) = true;

    % three arguments: object, dimensions, linespec
    plot(zB,[1,2],'r+');
    resvec(end+1) = true;
    
    % three arguments: object, dimensions, NVpairs
    plot(zB,[1,2],'LineWidth',2);
    plot(zB,[1,2],'Color',[.6 .6 .6],'LineWidth',2);
    plot(zB,[1,2],'EdgeColor','k','FaceColor',[.8 .8 .8]);
    resvec(end+1) = true;
    
    % four arguments: object, dimensions, linespec, NVpairs
    plot(zB,[1,2],'r','LineWidth',2);
    plot(zB,[1,2],'r','LineWidth',2,'EdgeColor',[.6 .6 .6]);
    resvec(end+1) = true;

    % check if plotted correctly
    ax = gca();
    colorOrder = ax.ColorOrder;
    
    % plot set
    plot(zB,[1,2]);
    V = [7 7 5 -1 -4 -4 -2 4 7
        -3 -3 -1 1 1 1 -1 -3 -3];
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
