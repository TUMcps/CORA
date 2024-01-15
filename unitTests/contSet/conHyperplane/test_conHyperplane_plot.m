function res = test_conHyperplane_plot
% test_conHyperplane_plot - unit test function of plot;
%    this function aims to go through many variations of input arguments
%    note: only run-time errors checked, manual bug check necessary
%
% Syntax:
%    res = test_conHyperplane_plot
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

% instantiate constrained hyperplane
a = [1, 1]; b = 2;
C = [1, 0];
d = 2.5;
hyp = conHyperplane(a, b, C, d);

try
    % try all variations in plotting
    figure;
    xlim([-2,2]);
    ylim([-2,2]);
    
    % one argument: object
    plot(hyp);
    resvec(end+1) = true;
    
    % two arguments: object, dimensions
    plot(hyp,1);
    plot(hyp,[1,2]);
    resvec(end+1) = true;
    
    % three arguments: object, dimensions, NVpairs
    plot(hyp,[1,2],'LineWidth',2);
    plot(hyp,[1,2],'Color',[.6 .6 .6],'LineWidth',2);
    plot(hyp,[1,2],'FaceColor',[.6 .6 .6]);
    resvec(end+1) = true;

    % check if plotted correctly
    ax = gca();
    colorOrder = ax.ColorOrder;
    
    % plot set
    plot(hyp,[1,2]);
    V = [2.5 1; -0.5 1];
    % check points
    resvec(end+1) = compareMatrices(V, [ax.Children(1).XData;ax.Children(1).YData],1e-14,'equal',true);
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
