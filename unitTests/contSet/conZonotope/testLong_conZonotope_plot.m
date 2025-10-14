function res = testLong_conZonotope_plot
% testLong_conZonotope_plot - unit test function of plot;
%    this function aims to go through many variations of input arguments
%    note: only run-time errors checked, manual bug check necessary
%
% Syntax:
%    res = testLong_conZonotope_plot
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
% Last update:   09-May-2023 (TL, added plotted point checks)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% construct zonotope
Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
A = [1 1 1];
b = 1;
cZ = conZonotope(Z,A,b);

try
    % try all variations in plotting
    figure;
    
    % one argument: object
    plot(cZ);
    
    % two arguments: object, dimensions
    plot(cZ,1);
    plot(cZ,[1,2]);
    
    % three arguments: object, dimensions, linespec
    plot(cZ,[1,2],'r+');
    
    % three arguments: object, dimensions, NVpairs
    plot(cZ,[1,2],'LineWidth',2);
    plot(cZ,[1,2],'Color',[.6 .6 .6],'LineWidth',2);
    
    % three arguments: object, dimensions, NVpair 'Splits'
    plot(cZ,[1,2],'Splits',4);
    plot(cZ,[1,2],'Splits',4,'LineWidth',2);
    plot(cZ,[1,2],'Splits',4,'EdgeColor','k','FaceColor',[.8 .8 .8]);

    % three arguments: object, dimensions, NVpair 'Template'
    plot(cZ,[1,2],'Template',16);
    plot(cZ,[1,2],'Template',16,'LineWidth',2);
    plot(cZ,[1,2],'Template',16,'EdgeColor','k','FaceColor',[.8 .8 .8]);
    
    % four arguments: object, dimensions, linespec, NVpairs
    plot(cZ,[1,2],'r','LineWidth',2);
    plot(cZ,[1,2],'r','LineWidth',2,'EdgeColor',[.6 .6 .6]);

    % check if plotted correctly
    ax = gca();
    colorOrder = ax.ColorOrder;
    
    % plot set
    plot(cZ,[1,2]);
    V = [3.5 -2.5 -0.5 3.5
        -0.5 -1.5 2.5 -0.5];
    % check points
    assert(compareMatrices(V, readVerticesFromFigure(ax.Children(1)),1e-4,'equal',true));
    % test color
    if CORA_PLOT_FILLED
        assert(isequal(colorOrder(1,:), ax.Children(1).EdgeColor));
        assert(isequal(colorOrder(1,:), ax.Children(1).FaceColor));
    else
        assert(isequal(colorOrder(1,:), ax.Children(1).Color));
    end

    
    % close figure
    close;
catch ME
    close;
    rethrow(ME)
end

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
