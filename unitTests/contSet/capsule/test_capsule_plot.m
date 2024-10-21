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
    
    % three arguments: object, dimensions, linespec
    plot(C,[1,2],'r+');
    
    % three arguments: object, dimensions, NVpairs
    plot(C,[1,2],'LineWidth',2);
    plot(C,[1,2],'Color',[.6 .6 .6],'LineWidth',2);
    plot(C,[1,2],'FaceColor',[.6 .6 .6]);
    
    % four arguments: object, dimensions, linespec, NVpairs
    plot(C,[1,2],'FaceColor','r','LineWidth',2);
    plot(C,[1,2],'FaceColor','r','LineWidth',2,'EdgeColor',[.6 .6 .6]);
   
    close;

    % check if plotted correctly
    figure; hold on;
    ax = gca();
    colorOrder = ax.ColorOrder;
    
    % plot set
    plot(C,[1,2]);
    V = [ ...
     1.5000000000000000, 1.4059835861602852, 1.1563036824392341, 0.8466895194061126, 0.5958618798622533, 0.5000000000000000, 0.5940164138397146, 0.8436963175607657, 1.1533104805938872, 1.4041381201377467, 1.5000000000000000 ; ...
     2.0000000000000000, -0.2918515509097633, -0.4749412162109486, -0.4759158502719480, -0.2944017320797011, 0.0000000000000001, 2.2918515509097634, 2.4749412162109485, 2.4759158502719480, 2.2944017320797014, 2.0000000000000000 ; ...
    ];

    % check points
    assert(compareMatrices(V, [ax.Children(1).XData;ax.Children(1).YData],1e-5,'subset',true));
    % test color
    assert(isequal(colorOrder(1,:), ax.Children(1).Color));

    % close figure
    close;
catch ME
    close;
    rethrow(ME)
end

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
