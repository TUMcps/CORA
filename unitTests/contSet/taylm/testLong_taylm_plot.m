function res = testLong_taylm_plot
% testLong_taylm_plot - unit test of plot function
%
% Syntax:
%    res = testLong_taylm_plot
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: taylm, polyZonotope
% Subfunctions: none
% MAT-files required: none

% Authors:       Tobias Ladner
% Written:       10-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

try

    figure;
    
    % 1d case
    tx = taylm(interval(-2,2),4,'x');
    plot(tx,1);
    resvec(end+1) = true;
    
    % 2d case
    tay = [tx; sin(tx)];
    plot(tay,1:2);
    resvec(end+1) = true;
    
    plot(tay,[2,1]);
    resvec(end+1) = true;
    
    % three arguments: object, dimensions, linespec
    plot(tay,[1,2],'r+');
    resvec(end+1) = true;
    
    % three arguments: object, dimensions, NVpairs
    plot(tay,[1,2],'LineWidth',2);
    plot(tay,[1,2],'Color',[.6 .6 .6],'LineWidth',2);
    plot(tay,[1,2],'EdgeColor','k','FaceColor',[.8 .8 .8]);
    resvec(end+1) = true;
    
    % four arguments: object, dimensions, linespec, NVpairs
    plot(tay,[1,2],'r','LineWidth',2);
    plot(tay,[1,2],'r','LineWidth',2,'EdgeColor',[.6 .6 .6]);
    resvec(end+1) = true;

    % check if plotted correctly
    ax = gca();
    colorOrder = ax.ColorOrder;

    % plot set
    plot(tay,[1,2]);
    V = [ ...
        2.000 1.9844 0.5781 -0.1406 -0.7656 -1.4805 -1.6875 -0.9180 -0.0469 1.2891 1.8984 2.0000; ...
        0.400 0.4154 0.2793 -0.4068 -0.9575 -1.2063 -0.6199 -0.5224  0.2198 1.1987 1.0248 0.4000; ...
        ];
    % check points
    resvec(end+1) = compareMatrices(V, [ax.Children(1).XData;ax.Children(1).YData],1e-4,'subset',true);
    % test color
    resvec(end+1) = isequal(colorOrder(1,:), ax.Children(1).Color);
    
    % close figure
    close

catch ME
    resvec(end+1) = false;

    % close figure
    close
end

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
