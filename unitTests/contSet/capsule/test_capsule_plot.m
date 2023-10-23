function res = test_capsule_plot
% test_capsule_plot - unit test function of plot; this function aims
%    to go through many variations of input arguments
%    note: only run-time errors checked, go through manually to check for bugs
%
% Syntax:
%    res = test_capsule_plot
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

resvec = [];

% instantiate capsule
C = capsule([1;1;2],[0;1;-0.5],0.5);

try
    % try all variations in plotting
    figure;
    
    % one argument: object
    plot(C);
    
    % two arguments: object, dimensions
    plot(C,1);
    plot(C,[1,2]);
    plot(C,[2,3]);
    resvec(end+1) = true;
    
    % three arguments: object, dimensions, linespec
    plot(C,[1,2],'r+');
    
    % three arguments: object, dimensions, NVpairs
    plot(C,[1,2],'LineWidth',2);
    plot(C,[1,2],'Color',[.6 .6 .6],'LineWidth',2);
    plot(C,[1,2],'FaceColor',[.6 .6 .6]);
    
    % four arguments: object, dimensions, linespec, NVpairs
    plot(C,[1,2],'FaceColor','r','LineWidth',2);
    plot(C,[1,2],'FaceColor','r','LineWidth',2,'EdgeColor',[.6 .6 .6]);
    resvec(end+1) = true;
   
    close;

    % check if plotted correctly
    figure; hold on;
    ax = gca();
    colorOrder = ax.ColorOrder;
    
    % plot set
    plot(C,[1,2]);
    V = [1.5 1.38682 1.01731 0.660629 0.5 0.5 0.588577 1.00472 1.46239 1.5 1.5; 
        2 2.31682 2.4997 2.36719 2 0 -0.284132 -0.499978 -0.190251 0 2];
    % check points
    resvec(end+1) = compareMatrices(V, [ax.Children(1).XData;ax.Children(1).YData],1e-5,'subset',true);
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
