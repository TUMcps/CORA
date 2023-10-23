function res = testLong_polyZonotope_plot
% testLong_polyZonotope_plot - unit test function of plot;
%    this function aims to go through many variations of input arguments
%    note: only run-time errors checked, manual bug check necessary
%
% Syntax:
%    res = testLong_polyZonotope_plot
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

resvec = [];

% instantiate polynomial zonotope
c = rand(4,1)-0.5*ones(4,1);
G = rand(4,6)-0.5*ones(4,6);
ind = datasample(1:6,4,'Replace',false);
G(:,ind) = G(:,ind)./10;
GI = rand(4,2)-0.5*ones(4,2);
E = [eye(4), round(rand(4,2)*5)];
pZ = polyZonotope(c,G,GI,E);

try
    % try all variations in plotting
    figure;
    
    % one argument: object
    plot(pZ);
    resvec(end+1) = true;
    
    % two arguments: object, dimensions
    plot(pZ,1);
    plot(pZ,[1,2]);
    plot(pZ,[2,3]);
    resvec(end+1) = true;
    
    % three arguments: object, dimensions, linespec
    plot(pZ,[1,2],'r+');
    resvec(end+1) = true;
    
    % three arguments: object, dimensions, NVpairs
    plot(pZ,[1,2],'LineWidth',2);
    plot(pZ,[1,2],'Color',[.6 .6 .6],'LineWidth',2);
    resvec(end+1) = true;
    
    % three arguments: object, dimensions, NVpair 'Splits'
    plot(pZ,[1,2],'Splits',0);
    plot(pZ,[1,2],'Splits',6);
    plot(pZ,[1,2],'Splits',6,'LineWidth',2);
    plot(pZ,[1,2],'Splits',6,'EdgeColor','k','FaceColor',[.8 .8 .8]);
    resvec(end+1) = true;
    
    % four arguments: object, dimensions, linespec, NVpairs
    plot(pZ,[1,2],'FaceColor','r','LineWidth',2);
    plot(pZ,[1,2],'FaceColor','r','LineWidth',2,'EdgeColor',[.6 .6 .6]);
    resvec(end+1) = true;

    % the polyZonotope
    c = [4;4];
    G = [2 1 2; 0 2 2];
    E = [1 0 3;0 1 1];
    GI = [0.5;0];
    pZ = polyZonotope(c,G,GI,E);

    % check if plotted correctly
    ax = gca();
    colorOrder = ax.ColorOrder;
    
    % plot set
    plot(pZ,[1,2]);
    V = [8.0341 8.7500 4.0784 1.8035 2.6281 5.7534 8.0341;
        7.6591 7.0000 0.7659 2.7770 5.4406 6.2534 7.6591];
    % check points
    resvec(end+1) = compareMatrices(V, [ax.Children(1).XData;ax.Children(1).YData],1e-4,'subset',true);
    % test color
    resvec(end+1) = isequal(colorOrder(1,:), ax.Children(1).Color);

    % test line polyZonotope
    plot(polyZonotope([1;1], [1; 1]))
    ax = gca();
    resvec(end+1) = compareMatrices([0 2 0; 0 2 0], ...
        [ax.Children(1).XData; ax.Children(1).YData], 1e-8,"equal");

    % test almost line polyZonotope
    plot(polyZonotope([1;1], [1 0; 0 0.000000000000001]))
    ax = gca();
    resvec(end+1) = compareMatrices([0 2 2 0 0; 1 1 1 1 1], ...
        [ax.Children(1).XData; ax.Children(1).YData], 1e-8,"equal");
    
    % close figure
    close;
catch ME
    close;
    resvec(end+1) = false;
end

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
