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

resvec = [];

% instantiate zonotope
Z = zonotope([1;-1;2],[2 -1 3; 0 1 -1; -1 4 2]);

try
    % try all variations in plotting
    figure;
    
    % one argument: object
    plot(Z);
    resvec(end+1) = true;
    
    % two arguments: object, dimensions
    plot(Z,1);
    plot(Z,[1,2]);
    plot(Z,[2,3]);
    resvec(end+1) = true;
    
    % three arguments: object, dimensions, linespec
    plot(Z,[1,2],'r+');
    resvec(end+1) = true;
    
    % three arguments: object, dimensions, NVpairs
    plot(Z,[1,2],'LineWidth',2);
    plot(Z,[1,2],'Color',[.6 .6 .6],'LineWidth',2);
    plot(Z,[1,2],'EdgeColor','k','FaceColor',[.8 .8 .8]);
    resvec(end+1) = true;
    
    % four arguments: object, dimensions, linespec, NVpairs
    plot(Z,[1,2],'r','LineWidth',2);
    plot(Z,[1,2],'r','LineWidth',2,'EdgeColor',[.6 .6 .6]);
    resvec(end+1) = true;

    % check if plotted correctly
    ax = gca();
    colorOrder = ax.ColorOrder;
    
    % plot set
    plot(Z,[1,2]);
    V = [3 7 5 -1 -5 -3 3
        -3 -3 -1 1 1 -1 -3];
    % check points
    resvec(end+1) = compareMatrices(V, [ax.Children(1).XData;ax.Children(1).YData],1e-4,'equal',true);
    % test color
    resvec(end+1) = isequal(colorOrder(1,:), ax.Children(1).Color);

    % check FaceAlpha
    plot(Z, 1:2, 'r', 'FaceAlpha', 0.1)
    resvec(end+1) = isequal([1 0 0], ax.Children(1).FaceColor) ...
        && isequal([1 0 0], ax.Children(1).EdgeColor) ...
        && isequal(0.1, ax.Children(1).FaceAlpha) ...
        && isequal(1, ax.Children(1).EdgeAlpha);
    
    % close figure
    close;
catch ME
    close;
    resvec(end+1) = false;
end

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
