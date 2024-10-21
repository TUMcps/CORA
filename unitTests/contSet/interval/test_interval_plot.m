function res = test_interval_plot
% test_interval_plot - unit test function of plot; this function aims
%    to go through many variations of input arguments
%    note: only run-time errors checked
%
% Syntax:
%    res = test_interval_plot
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
% Last update:   08-May-2023 (TL, added plotted point checks)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

I = interval([1;1;2],[3;4;7]);

try
    % try all variations in plotting
    figure;
    
    % one argument: object
    plot(I);
    
    % two arguments: object, dimensions
    plot(I,1);
    plot(I,[1,2]);
    plot(I,[2,3]);
    
    % three arguments: object, dimensions, linespec
    plot(I,[1,2],'r+');
    
    % three arguments: object, dimensions, NVpairs
    plot(I,[1,2],'LineWidth',2);
    plot(I,[1,2],'Color',[.6 .6 .6],'LineWidth',2);
    plot(I,[1,2],'EdgeColor','k','FaceColor',[.8 .8 .8]);
    
    % four arguments: object, dimensions, linespec, NVpairs
    plot(I,[1,2],'r','LineWidth',2);
    plot(I,[1,2],'r','LineWidth',2,'EdgeColor',[.6 .6 .6]);

    % plot 3d
    plot(I,[1,2,3]);

    close;

    % check if plotted correctly
    figure; hold on;
    ax = gca();
    colorOrder = ax.ColorOrder;
    
    % plot first set
    plot(I,[1,2]);
    V = [1 1 3 3 1; 1 4 4 1 1];
    % check points
    assert(compareMatrices(V, [ax.Children(1).XData;ax.Children(1).YData],1e-4,'equal',true));
    % test color
    assert(isequal(colorOrder(1,:), ax.Children(1).Color));

    % plot second set
    plot(I,[1,3]);
    V = [1 1 3 3 1; 2 7 7 2 2];
    % check points
    assert(compareMatrices(V, [ax.Children(1).XData;ax.Children(1).YData],1e-4,'equal',true));
    % test color
    assert(isequal(colorOrder(2,:), ax.Children(1).Color));

    % plot 3d set
    plot(I,[1,2,3]);
    V = [ ...
     1.000, 3.000, 3.000, 1.000, 1.000 ; ...
     1.000, 1.000, 4.000, 4.000, 1.000 ; ...
     2.000, 2.000, 2.000, 2.000, 2.000 ; ...
    ];
    % check points (only first facet)
    assert(compareMatrices(V, [ax.Children(1).XData;ax.Children(1).YData;ax.Children(1).ZData],1e-4,'equal',true));
    % test color
    assert(isequal(colorOrder(3,:), ax.Children(1).Color));
    
    % close figure
    close;

    % test interval with infinity bounds
    figure;

    % set bounds
    xlim([1,2]);
    ylim([-2,3]);
    ax = gca();

    % plot interval with all inf bounds
    I = interval([-Inf;-Inf],[Inf;Inf]);
    plot(I);

    % check points
    V = [1 1 2 2 1; -2 3 3 -2 -2];
    assert(compareMatrices(V, [ax.Children(1).XData';ax.Children(1).YData'],1e-4,'equal',true));

    % plot interval with some inf bounds
    I = interval([1.5;-Inf],[Inf;2]);
    plot(I);

    % check points
    V = [1.5 1.5 2 2 1.5; -2 2 2 -2 -2];
    assert(compareMatrices(V, [ax.Children(1).XData';ax.Children(1).YData'],1e-4,'equal',true));

    % plot interval outside of xlim
    I = interval([-Inf;4],[2;Inf]);
    plot(I);

    % check points
    V = [1 2 1; 4 4 4];
    assert(compareMatrices(V, [ax.Children(1).XData';ax.Children(1).YData'],1e-4,'equal',true));

    % check single point
    p = [1.5;1];
    I = interval(p);
    plot(I);
    assert(compareMatrices(p, [ax.Children(1).XData;ax.Children(1).YData],1e-4,'equal',true));

    % close figure
    close;
    

catch ME
    close;
    rethrow(ME)
end

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
