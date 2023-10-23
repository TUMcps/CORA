function res = test_halfspace_plot
% test_halfspace_plot - unit test function of plot;
%    this function aims to go through many variations of input arguments
%    note: only run-time errors checked, manual bug check necessary
%
% Syntax:
%    res = test_halfspace_plot
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

% normal vector and offset
a = [-0.5; -1; 0.1];
b = -1;

% instantiate constrained hyperplane
hs = halfspace(a,b);

try
    % try all variations in plotting
    figure;
    xlim([-2,3]);
    ylim([-2,2]);
    
    % one argument: object
    plot(hs);
    resvec(end+1) = true;
    
    % two arguments: object, dimensions
    plot(hs,1);
    plot(hs,[1,2]);
    resvec(end+1) = true;

    % check if plotted correctly
    ax = gca();
    colorOrder = ax.ColorOrder;
    
    % plot set
    plot(hs,[1,2]);
    V = [-2 2; 3 -0.5; 3 2; -2 2];
    % check points
    resvec(end+1) = compareMatrices(V, [ax.Children(1).Vertices],1e-4,'equal',true);
    % test color
    resvec(end+1) = isequal(colorOrder(1,:), ax.Children(1).FaceColor);

    % check edge color
    hs = halfspace([1;1],1);
    plot(hs,[1,2],'EdgeColor',[1 0 0]);
    ax = gca();
    resvec(end+1) = isequal(ax.Children(1).EdgeColor, [1 0 0]);

    % but default edge color should be none
    hs = halfspace([1;1],1);
    plot(hs);
    ax = gca();
    resvec(end+1) = isequal(ax.Children(1).EdgeColor, 'none');
  
    % close figure
    close;
catch ME
    close;
    resvec(end+1) = false;
end

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
