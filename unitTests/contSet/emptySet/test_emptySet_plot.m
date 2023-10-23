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

resvec = [];

O = emptySet(3);

try
    % try all variations in plotting
    figure;
    
    % one argument: object
    plot(O);
    resvec(end+1) = true;
    
    % two arguments: object, dimensions
    plot(O,1);
    plot(O,[1,2]);
    plot(O,[2,3]);
    resvec(end+1) = true;
    
    % three arguments: object, dimensions, linespec
    plot(O,[1,2],'r+');
    resvec(end+1) = true;
    
    % three arguments: object, dimensions, NVpairs
    plot(O,[1,2],'LineWidth',2);
    plot(O,[1,2],'Color',[.6 .6 .6],'LineWidth',2);
    plot(O,[1,2],'EdgeColor','k','FaceColor',[.8 .8 .8]);
    resvec(end+1) = true;
    
    % four arguments: object, dimensions, linespec, NVpairs
    plot(O,[1,2],'r','LineWidth',2);
    plot(O,[1,2],'r','LineWidth',2,'EdgeColor',[.6 .6 .6]);
    resvec(end+1) = true;

    close;

    % check if plotted correctly
    figure; hold on;
    ax = gca();
    colorOrder = ax.ColorOrder;
    
    % plot first set
    plot(O,[1,2]);
    V = zeros(2,0);
    % check points
    resvec(end+1) = all(isnan([ax.Children(1).XData;ax.Children(1).YData]));
    % test color
    resvec(end+1) = isequal(colorOrder(1,:), ax.Children(1).Color);

    % plot second set
    plot(O,[1,3]);
    V = zeros(2,0);
    % check points
    resvec(end+1) = all(isnan([ax.Children(1).XData;ax.Children(1).YData]));
    % test color
    resvec(end+1) = isequal(colorOrder(2,:), ax.Children(1).Color);
    
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
    resvec(end+1) = all(isnan([ax.Children(1).XData;ax.Children(1).YData]));

    % close figure
    close;

catch ME
    close;
    resvec(end+1) = false;
end

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
