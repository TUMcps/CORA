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

% Author:       Mark Wetzlinger
% Written:      04-August-2020
% Last update:  08-May-2023 (TL: added plotted point checks)
% Last revision:---

%------------- BEGIN CODE --------------

resvec = [];

I = interval([1;1;2],[3;4;7]);

try
    % try all variations in plotting
    figure;
    
    % one argument: object
    plot(I);
    resvec(end+1) = true;
    
    % two arguments: object, dimensions
    plot(I,1);
    plot(I,[1,2]);
    plot(I,[2,3]);
    resvec(end+1) = true;
    
    % three arguments: object, dimensions, linespec
    plot(I,[1,2],'r+');
    resvec(end+1) = true;
    
    % three arguments: object, dimensions, NVpairs
    plot(I,[1,2],'LineWidth',2);
    plot(I,[1,2],'Color',[.6 .6 .6],'LineWidth',2);
    plot(I,[1,2],'EdgeColor','k','FaceColor',[.8 .8 .8]);
    resvec(end+1) = true;
    
    % four arguments: object, dimensions, linespec, NVpairs
    plot(I,[1,2],'r','LineWidth',2);
    plot(I,[1,2],'r','LineWidth',2,'EdgeColor',[.6 .6 .6]);
    resvec(end+1) = true;

    close;

    % check if plotted correctly
    figure; hold on;
    ax = gca();
    colorOrder = ax.ColorOrder;
    
    % plot first set
    plot(I,[1,2]);
    V = [1 3 3 1 1; 1 1 4 4 1];
    % check points
    resvec(end+1) = compareMatrices(V, [ax.Children(1).XData;ax.Children(1).YData],1e-4,'equal',true);
    % test color
    resvec(end+1) = isequal(colorOrder(1,:), ax.Children(1).Color);

    % plot second set
    plot(I,[1,3]);
    V = [1 3 3 1 1; 2 2 7 7 2];
    % check points
    resvec(end+1) = compareMatrices(V, [ax.Children(1).XData;ax.Children(1).YData],1e-4,'equal',true);
    % test color
    resvec(end+1) = isequal(colorOrder(2,:), ax.Children(1).Color);
    
    % close figure
    close;
catch ME
    close;
    resvec(end+1) = false;
end

% gather results
res = all(resvec);

%------------- END OF CODE --------------
