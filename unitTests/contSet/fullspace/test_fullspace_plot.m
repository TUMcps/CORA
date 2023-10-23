function res = test_fullspace_plot
% test_fullspace_plot - unit test function of plot
%
% Syntax:
%    res = test_fullspace_plot
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

fs = fullspace(3);

try
    % try all variations in plotting
    figure;
    
    % one argument: object
    plot(fs);
    resvec(end+1) = true;
    
    % two arguments: object, dimensions
    plot(fs,1);
    plot(fs,[1,2]);
    plot(fs,[2,3]);
    resvec(end+1) = true;
    
    % three arguments: object, dimensions, linespec
    plot(fs,[1,2],'r+');
    resvec(end+1) = true;
    
    % three arguments: object, dimensions, NVpairs
    plot(fs,[1,2],'LineWidth',2);
    plot(fs,[1,2],'Color',[.6 .6 .6],'LineWidth',2);
    plot(fs,[1,2],'EdgeColor','k','FaceColor',[.8 .8 .8]);
    resvec(end+1) = true;
    
    % four arguments: object, dimensions, linespec, NVpairs
    plot(fs,[1,2],'r','LineWidth',2);
    plot(fs,[1,2],'r','LineWidth',2,'EdgeColor',[.6 .6 .6]);
    resvec(end+1) = true;

    close;

    % check if plotted correctly
    figure; hold on;
    ax = gca();
    colorOrder = ax.ColorOrder;
    
    % plot first set
    plot(fs,[1,2]);
    V = [0 1 1 0 0; 0 0 1 1 0];
    % check points
    resvec(end+1) = compareMatrices(V, [reshape(ax.Children(1).XData,1,[]);reshape(ax.Children(1).YData,1,[])],1e-4,'equal',true);
    % test color
    resvec(end+1) = isequal(colorOrder(1,:), ax.Children(1).FaceColor);

    % plot second set
    plot(fs,[1,3]);
     V = [0 1 1 0 0; 0 0 1 1 0];
    % check points
    resvec(end+1) = compareMatrices(V, [reshape(ax.Children(1).XData,1,[]);reshape(ax.Children(1).YData,1,[])],1e-4,'equal',true);
    % test color
    resvec(end+1) = isequal(colorOrder(2,:), ax.Children(1).FaceColor);
    
    % close figure
    close;

    % test interval with infinity bounds
    figure;

    % set bounds
    xlim([1,2]);
    ylim([-2,3]);
    ax = gca();

    % plot interval with all inf bounds
    fs = interval([-Inf;-Inf],[Inf;Inf]);
    plot(fs);

    % check points
    V = [1 2 2 1 1; -2 -2 3 3 -2];
    resvec(end+1) = compareMatrices(V, [ax.Children(1).XData;ax.Children(1).YData],1e-4,'equal',true);

    % close figure
    close;
    

catch ME
    close;
    resvec(end+1) = false;
end

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
