function res = test_conPolyZono_plot
% test_conPolyZono_plot - unit test function of plot;
%    this function aims to go through many variations of input arguments
%    note: only run-time errors checked, manual bug check necessary
%
% Syntax:
%    res = test_conPolyZono_plot
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

% construct constrained polynomial zonotope
c = [0;0];
G = [1 0 1 -1; 0 1 1 1];
E = [1 0 1 2; 0 1 1 0; 0 0 1 1];
A = [1 -0.5 0.5];
b = 0.5;
EC = [0 1 2; 1 0 0; 0 1 0];
cPZ = conPolyZono(c,G,E,A,b,EC);

try
    % try all variations in plotting
    figure;
    
    % one argument: object
    plot(cPZ);
    resvec(end+1) = true;
    
    % two arguments: object, dimensions
    plot(cPZ,1);
    plot(cPZ,[1,2]);
    resvec(end+1) = true;
    
    % three arguments: object, dimensions, NVpairs
    plot(cPZ,[1,2],'LineWidth',2);
    plot(cPZ,[1,2],'Color',[.6 .6 .6],'LineWidth',2);
    plot(cPZ,[1,2],'FaceColor',[.6 .6 .6]);
    resvec(end+1) = true;

    % check if plotted correctly
    ax = gca();
    colorOrder = ax.ColorOrder;
    
    % plot set
    plot(cPZ,[1,2]);
    V = [0.504898 0.497589 0.508057 0.645584 0.802383 1.58653 2.49087 1.01437 0.444302 -0.869492 -1.5 -0.886902 -0.00170898 0.150778 0.461334 0.505696 0.504898;
        1.69337 1.71927 2.00171 1.55672 0.393692 -0.592482 -1.00864 -0.397009 -0.0164938 -0.0579681 1 0.347473 0.72876 0.643607 0.983759 1.69184 1.69337];
    % check points
    resvec(end+1) = compareMatrices(V, [ax.Children(1).XData;ax.Children(1).YData],1e-4,'subset',true);
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
