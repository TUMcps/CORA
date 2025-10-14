function res = test_emptySet_plot
% test_emptySet_plot - unit test function of plot
%
% Syntax:
%    res = test_emptySet_plot
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

% Authors:       Tobias Ladner
% Written:       03-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

O = emptySet(3);

try
    % try all variations in plotting
    figure;
    
    % one argument: object
    plot(O);
    
    % two arguments: object, dimensions
    plot(O,1);
    plot(O,[1,2]);
    plot(O,[2,3]);
    
    % three arguments: object, dimensions, linespec
    plot(O,[1,2],'r+');
    
    % three arguments: object, dimensions, NVpairs
    plot(O,[1,2],'LineWidth',2);
    plot(O,[1,2],'Color',[.6 .6 .6],'LineWidth',2);
    plot(O,[1,2],'EdgeColor','k','FaceColor',[.8 .8 .8]);
    
    % four arguments: object, dimensions, linespec, NVpairs
    plot(O,[1,2],'r','LineWidth',2);
    plot(O,[1,2],'r','LineWidth',2,'EdgeColor',[.6 .6 .6]);

    close;

    % check if plotted correctly
    figure; hold on;
    ax = gca();
    colorOrder = ax.ColorOrder;
    
    % plot first set
    plot(O,[1,2]);
    V = zeros(2,0);
    % check points
    assert(all(isnan(readVerticesFromFigure(ax.Children(1)))));
    % test color
    if CORA_PLOT_FILLED
        assert(isequal(colorOrder(1,:), ax.Children(1).EdgeColor));
        assert(isequal(colorOrder(1,:), ax.Children(1).FaceColor));
    else
        assert(isequal(colorOrder(1,:), ax.Children(1).Color));
    end

    % plot second set
    plot(O,[1,3]);
    V = zeros(2,0);
    % check points
    assert(all(isnan(readVerticesFromFigure(ax.Children(1)))));
    
    % close figure
    close;

    % test interval with infinity bounds
    figure;

    % set bounds
    xlim([1,2]);
    ylim([-2,3]);
    ax = gca();

    % plot interval with all inf bounds
    V = zeros(2,0);
    plot(O);

    % check points
    V = zeros(2,0);
    assert(all(isnan(readVerticesFromFigure(ax.Children(1)))));

    % close figure
    close;

catch ME
    close;
    rethrow(ME)
end

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
