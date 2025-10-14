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
    assert(compareMatrices(V, readVerticesFromFigure(ax.Children(1)),1e-4,'subset',true));
    % test color
    if CORA_PLOT_FILLED
        assert(isequal(colorOrder(1,:), ax.Children(1).EdgeColor));
        assert(isequal(colorOrder(1,:), ax.Children(1).FaceColor));
    else
        assert(isequal(colorOrder(1,:), ax.Children(1).Color));
    end

    % test line polyZonotope
    plot(polyZonotope([1;1], [1; 1]))
    ax = gca();
    assert(compareMatrices([0 2 0; 0 2 0], readVerticesFromFigure(ax.Children(1)), 1e-8,"equal"));

    % test almost line polyZonotope
    plot(polyZonotope([1;1], [1 0; 0 0.000000000000001]))
    ax = gca();
    % points should reflect that set is not just a single line
    assert(compareMatrices([0 2 2 0 0; 1 1 1 1 1], readVerticesFromFigure(ax.Children(1)), 1e-8,"equal"));

    % test polyZonotope with holes
    pZ = polyZonotope([0;0],[1 0 -1 0; 0 1 0 -1], 0.01*eye(2), [1,3,7,15]);
    plot(pZ,1:2,'Filled',false);
    V = readVerticesFromFigure(ax.Children(1));
    % one nan column per hole (+ the first value...)
    assert(nnz(all(isnan(V))) == 2 + 1);
    % check presence of subset of nodes ---
    nanIdx = find(all(isnan(V)));
    % outer line 
    V_outer = V(:,(nanIdx(1)+1):(nanIdx(2)-1));
    V_outer_true = [ 0.05250 0.09937 0.13062 0.17749 0.20873 0.25553 0.28667 0.33318 0.36396 0.40950 0.43920 0.48219 0.50930 0.54654 0.56817 0.59382 0.60490 0.60917 0.59362 0.58283 0.56037 0.54069 0.50324 0.47491 0.43592 0.36603 0.31100 0.21435 0.13953 0.02042 ; 0.01025 0.01131 0.01278 0.01659 0.02047 0.02874 0.03617 0.05062 0.06273 0.08508 0.10303 0.13497 0.15986 0.20295 0.23569 0.29080 0.33132 0.37630 0.44218 0.46284 0.49074 0.50612 0.52168 0.52385 0.51460 0.47829 0.43414 0.32726 0.22063 0.01003 ];
    assert(compareMatrices(V_outer_true,V_outer,1e-4,"subset"))
    % hole 1
    V_hole1 = V(:,(nanIdx(2)+1):(nanIdx(3)-1));
    V_hole1_true = [ -0.01000 0.07911 0.19461 0.26104 0.34624 0.39449 0.45508 0.50855 0.53698 0.57115 0.58886 0.60871 0.62177 0.62992 0.62990 0.62838 0.62510 0.62024 0.61398 0.60648 0.59787 0.58828 0.56663 0.54231 0.50223 0.47380 0.42952 0.39924 0.35319 0.32221 0.27553 0.24432 0.19749 0.16622 0.11937 0.08800 0.02125 -0.03597 -0.15831 -0.22885 -0.31963 -0.37121 -0.43632 -0.49247 -0.52337 -0.56083 -0.58048 -0.60293 -0.61398 -0.62947 -0.62991 -0.62839 -0.62511 -0.62025 -0.61399 -0.60648 -0.59787 -0.58828 -0.57783 -0.55475 -0.51599 -0.48815 -0.44445 -0.41443 -0.36861 -0.33771 -0.29111 -0.25993 -0.21311 -0.18185 -0.13500 -0.10366 -0.04125 ; 0.01005 0.17558 0.34821 0.42536 0.49883 0.52631 0.54414 0.54441 0.53788 0.51895 0.50204 0.47271 0.44054 0.37422 0.35420 0.33247 0.31125 0.29064 0.27073 0.25156 0.23317 0.21559 0.18286 0.15334 0.11488 0.09290 0.06500 0.04953 0.03055 0.02046 0.00869 0.00283 -0.00344 -0.00622 -0.00871 -0.00954 -0.01000 -0.09921 -0.29910 -0.39004 -0.47915 -0.51439 -0.54083 -0.54525 -0.54187 -0.52626 -0.51083 -0.48290 -0.46229 -0.39631 -0.35421 -0.33248 -0.31125 -0.29064 -0.27073 -0.25156 -0.23317 -0.21559 -0.19882 -0.16770 -0.12695 -0.10354 -0.07365 -0.05697 -0.03635 -0.02527 -0.01219 -0.00558 0.00166 0.00496 0.00807 0.00919 0.01000 ];
    assert(compareMatrices(V_hole1_true,V_hole1,1e-4,"subset"))
    % hole 2
    V_hole2 = V(:,(nanIdx(3)+1):end);
    V_hole2_true = [ -0.02125 -0.06812 -0.11500 -0.14625 -0.19311 -0.22434 -0.27111 -0.30221 -0.34860 -0.37924 -0.42443 -0.45378 -0.49593 -0.52224 -0.55773 -0.57774 -0.60006 -0.60815 -0.60651 -0.58861 -0.57622 -0.55103 -0.52931 -0.48845 -0.45614 -0.41421 -0.33941 -0.28074 -0.17804 -0.09875 ; -0.01003 -0.01048 -0.01195 -0.01382 -0.01838 -0.02288 -0.03225 -0.04052 -0.05641 -0.06960 -0.09374 -0.11299 -0.14704 -0.17344 -0.21891 -0.25327 -0.31071 -0.35253 -0.39850 -0.45266 -0.47264 -0.49884 -0.51247 -0.52418 -0.52049 -0.50586 -0.45854 -0.40456 -0.27804 -0.15420 ];
    assert(compareMatrices(V_hole2_true,V_hole2,1e-4,"subset"))

    % plot with face color
    plot(pZ,1:2,'FaceColor','b');
    V = readVerticesFromFigure(ax.Children(1));
    V_true = [ -0.01000 0.01005 ; 0.19461 0.34821 ; 0.34624 0.49883 ; 0.45508 0.54414 ; 0.53698 0.53788 ; 0.58886 0.50204 ; 0.62177 0.44054 ; 0.62990 0.35420 ; 0.62510 0.31125 ; 0.61398 0.27073 ; 0.59787 0.23317 ; 0.56663 0.18286 ; 0.50223 0.11488 ; 0.42952 0.06500 ; 0.35319 0.03055 ; 0.27553 0.00869 ; 0.19749 -0.00344 ; 0.11937 -0.00871 ; 0.02125 -0.01000 ; 0.09937 0.01131 ; 0.17749 0.01659 ; 0.25553 0.02874 ; 0.33318 0.05062 ; 0.40950 0.08508 ; 0.48219 0.13497 ; 0.54654 0.20296 ; 0.59382 0.29082 ; 0.60652 0.39850 ; 0.57622 0.47264 ; 0.52933 0.51247 ; 0.47490 0.52383 ; 0.36602 0.47826 ; 0.21433 0.32721 ; 0.02125 0.01003 ; -0.11987 -0.24187 ; -0.29129 -0.45485 ; -0.41621 -0.53500 ; -0.52337 -0.54189 ; -0.58048 -0.51084 ; -0.61368 -0.46229 ; -0.62992 -0.35422 ; -0.62511 -0.31126 ; -0.61400 -0.27073 ; -0.59788 -0.23317 ; -0.57784 -0.19882 ; -0.52936 -0.13976 ; -0.45922 -0.08295 ; -0.38396 -0.04266 ; -0.30667 -0.01611 ; -0.22871 -0.00043 ; -0.15059 0.00725 ; -0.07250 0.00980 ; -0.08375 -0.01082 ; -0.16187 -0.01508 ; -0.23994 -0.02562 ; -0.31772 -0.04533 ; -0.39442 -0.07705 ; -0.46812 -0.12362 ; -0.53468 -0.18780 ; -0.58634 -0.27164 ; -0.60933 -0.36948 ; -0.58861 -0.45266 ; -0.55102 -0.49884 ; -0.48838 -0.52418 ; -0.39094 -0.49391 ; -0.24856 -0.36919 ; -0.05562 -0.07776 ; -0.01000 0.01001 ]';
    assert(compareMatrices(V_true,V,1e-4,"subset"))

    % close figure
    close;
catch ME
    close;
    rethrow(ME)
end

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
