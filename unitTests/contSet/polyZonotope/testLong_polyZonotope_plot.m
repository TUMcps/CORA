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
    
    % two arguments: object, dimensions
    plot(pZ,1);
    plot(pZ,[1,2]);
    plot(pZ,[2,3]);
    
    % three arguments: object, dimensions, linespec
    plot(pZ,[1,2],'r+');
    
    % three arguments: object, dimensions, NVpairs
    plot(pZ,[1,2],'LineWidth',2);
    plot(pZ,[1,2],'Color',[.6 .6 .6],'LineWidth',2);
    
    % three arguments: object, dimensions, NVpair 'Splits'
    plot(pZ,[1,2],'Splits',0);
    plot(pZ,[1,2],'Splits',6);
    plot(pZ,[1,2],'Splits',6,'LineWidth',2);
    plot(pZ,[1,2],'Splits',6,'EdgeColor','k','FaceColor',[.8 .8 .8]);
    
    % four arguments: object, dimensions, linespec, NVpairs
    plot(pZ,[1,2],'FaceColor','r','LineWidth',2);
    plot(pZ,[1,2],'FaceColor','r','LineWidth',2,'EdgeColor',[.6 .6 .6]);

    % 3d plot
    plot(pZ,[1,2,3],'Splits',4);
    plot(pZ,[1,2,3],'Splits',4,'FaceColor',CORAcolor("CORA:next"),'FaceAlpha',0.2);

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
    V = [ ...
     7.6320190429687500, 8.0000000000000000, 6.3209228515623499, 5.4114990234375000, 3.8965148925780624, 3.2341918945312500, 1.7291564941406250, 1.8480533676285997, 1.2078215279029001, 1.4117870330807998, 3.0229568481445312, 5.0312500000000000, 7.5898437500000000, 7.6320190429687500 ; ...
     7.3820190429687500, 6.0000000000000000, 3.7584228515624498, 2.5364990234375000, 0.5215148925781625, 1.7341918945312500, 2.2291564941406250, 3.2767874295594002, 3.9538536071778001, 4.6617870330807998, 5.6479568481445312, 6.0312500000000000, 7.3398437500000000, 7.3820190429687500 ; ...
    ];
    % check points
    assert(compareMatrices(V, [ax.Children(1).XData;ax.Children(1).YData],1e-4,'subset',true));
    % test color
    assert(isequal(colorOrder(1,:), ax.Children(1).Color));

    % test line polyZonotope
    plot(polyZonotope([1;1], [1; 1]))
    ax = gca();
    assert(compareMatrices([0 2 0; 0 2 0], ...
        [ax.Children(1).XData; ax.Children(1).YData], 1e-8,"equal"));

    % test almost line polyZonotope
    plot(polyZonotope([1;1], [1 0; 0 0.000000000000001]))
    ax = gca();
    % points should reflect that set is not just a single line
    assert(compareMatrices([0 2 2 0 0; 1 1 1 1 1], ...
        [ax.Children(1).XData; ax.Children(1).YData], 1e-8,"equal"));
    
    % close figure
    close;
catch ME
    close;
    rethrow(ME)
end

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
