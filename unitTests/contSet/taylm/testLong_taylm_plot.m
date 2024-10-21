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

try

    figure;
    
    % 1d case
    tx = taylm(interval(-2,2),4,'x');
    plot(tx,1);
    
    % 2d case
    tay = [tx; sin(tx)];
    plot(tay,1:2);
    
    plot(tay,[2,1]);
    
    % three arguments: object, dimensions, linespec
    plot(tay,[1,2],'r+');
    
    % three arguments: object, dimensions, NVpairs
    plot(tay,[1,2],'LineWidth',2);
    plot(tay,[1,2],'Color',[.6 .6 .6],'LineWidth',2);
    plot(tay,[1,2],'EdgeColor','k','FaceColor',[.8 .8 .8]);
    
    % four arguments: object, dimensions, linespec, NVpairs
    plot(tay,[1,2],'r','LineWidth',2);
    plot(tay,[1,2],'r','LineWidth',2,'EdgeColor',[.6 .6 .6]);

    % check if plotted correctly
    ax = gca();
    colorOrder = ax.ColorOrder;

    % plot set
    plot(tay,[1,2]);
    V = [ ...
     2.00000, 1.84375, 1.68750, 1.53125, 1.37500, 1.21875, 1.06250, 0.90625, 0.75000, 0.59375, 0.43750, 0.28125, 0.12500, -0.03125, -0.18750, -0.34375, -0.50000, -0.65625, -0.81250, -0.96875, -1.12500, -1.28125, -1.43750, -1.59375, -1.75000, -1.90625, -1.93750, -1.78125, -1.62500, -1.46875, -1.31250, -1.15625, -1.00000, -0.84375, -0.68750, -0.53125, -0.37500, -0.21875, -0.06250, 0.09375, 0.25000, 0.40625, 0.56250, 0.71875, 0.87500, 1.03125, 1.18750, 1.34375, 1.50000, 1.65625, 1.81250, 1.96875, 2.00000 ; ...
     0.40000, 0.53247, 0.61993, 0.66619, 0.67507, 0.65037, 0.59592, 0.51553, 0.41302, 0.29220, 0.15688, 0.01088, -0.14199, -0.29791, -0.45307, -0.60366, -0.74585, -0.87583, -0.98980, -1.08392, -1.15440, -1.19741, -1.20913, -1.18577, -1.12349, -1.01849, -0.45863, -0.57264, -0.64316, -0.67401, -0.66900, -0.63195, -0.56667, -0.47697, -0.36667, -0.23959, -0.09954, 0.04966, 0.20421, 0.36028, 0.51407, 0.66175, 0.79952, 0.92355, 1.03004, 1.11516, 1.17511, 1.20606, 1.20421, 1.16574, 1.08683, 0.96367, 0.40000 ; ...
    ];
    % check points
    assert(compareMatrices(V, [ax.Children(1).XData;ax.Children(1).YData],1e-4,'subset',true));
    % test color
    assert(isequal(colorOrder(1,:), ax.Children(1).Color));
    
    % close figure
    close

catch ME
    close
    rethrow(ME)
end

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
