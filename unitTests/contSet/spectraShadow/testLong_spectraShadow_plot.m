function res = testLong_spectraShadow_plot
% testLong_spectraShadow_plot - unit test function of plot
%
% Syntax:
%    res = testLong_spectraShadow_plot
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
% See also: polytope/plot

% Authors:       Adrian Kulmburg
% Written:       14-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% instantiate polytope (via conversion from zonotope)
Z = zonotope([1;-1;2],[2 -1 3; 0 1 -1; -1 4 2]);
SpS = spectraShadow(Z);

splits = 8;

try
    % try all variations in plotting
    clf;%figure;
    hold on
    axis equal
    xlim([-3 0]);
    ylim([0 2]);

    % one argument: object
    % Try differently positioned objects
    SpS = spectraShadow(ellipsoid(eye(2))) & interval([-3;0],[0;2]);
    plot(SpS);
    SpS = spectraShadow(ellipsoid([0.7 0;0 1])) & interval([-3;0],[0;2]);
    plot(SpS);

    % two arguments: object, dimensions
    p  = plot(SpS,1,'Splits',splits);
    % Try to delete object
    delete(p);
    SpS = spectraShadow(interval([-2.93;1.5],[-2.83;1.9]));
    plot(SpS,[1,2],'Splits',splits);
    SpS = spectraShadow(ellipsoid(0.0005*eye(3)) + [0;-0.66;0.68]);
    plot(SpS,[2,3],'Splits',splits);

    % three arguments: object, dimensions, NVpairs
    SpS = spectraShadow(polytope([-0.877526 1.67753; -1 1.67753; -0.9 1.5; -0.7775 1.5]'));
    plot(SpS,[1 2], 'FaceColor', colorblind('blue'), 'Splits', splits)
    SpS = spectraShadow(polytope([-1.12247 1.67753; -1 1.67753; -1.12247 1.5; -1 1.5]'));
    plot(SpS,[1 2], 'FaceColor', colorblind('blue'), 'Splits', splits)
    SpS = spectraShadow(polytope([-0.5 1.9; -0.4 1.9; -0.7 1.5; -0.6 1.5]'));
    plot(SpS,[1 2], 'FaceColor', colorblind('blue'), 'Splits', splits)
    SpS = spectraShadow(polytope([-0.5 1.9; -0.4 1.9; -0.3 1.5; -0.2 1.5]'));
    plot(SpS,[1 2], 'FaceColor', colorblind('blue'), 'Splits', splits)
    SpS = spectraShadow(polytope([-0.55 1.7; -0.35 1.7; -0.55 1.6; -0.35 1.6]'));
    plot(SpS,[1 2], 'FaceColor', colorblind('blue'), 'Splits', splits)

    SpS = spectraShadow(polytope([-0.760028 0.418087; -0.7625 0.411607; -0.78 0.42; -0.78 0.41]'));
    plot(SpS, [1 2], 'FaceColor', 'k', 'Splits', splits)
    SpS = spectraShadow(polytope([-0.8 0.42; -0.8 0.43; -0.78 0.42; -0.78 0.41]'));
    plot(SpS, [1 2], 'FaceColor', 'k', 'Splits', splits)
    SpS = spectraShadow(polytope([-0.8 0.42; -0.8 0.43; -0.83 0.45]'));
    plot(SpS, [1 2], 'FaceColor', 'k', 'Splits', splits)

    SpS = spectraShadow(ellipsoid(0.0005*eye(3)) + [0;-0.66;0.68]);
    plot(SpS,[2,3],'FaceColor', 'k','Splits',splits);

    % four arguments: object, dimensions, linespec, NVpairs
    SpS = spectraShadow(ellipsoid(0.05*eye(2)) + [-1.4;1.7]);
    plot(SpS,[1,2],'LineWidth',2,'EdgeColor',colorblind('red'),'Splits',splits);
    SpS = convHull(spectraShadow(ellipsoid(0.01*eye(2), [-2.5;1.8])), zonotope([-2.4;1.5]));
    plot(SpS,[1,2],'LineWidth',2,'FaceColor',colorblind('red'),'Splits',2*splits);
    SpS = convHull(spectraShadow(ellipsoid(0.01*eye(2), [-2.3;1.8])), zonotope([-2.4;1.5]));
    plot(SpS,[1,2],'LineWidth',2,'FaceColor',colorblind('red'),'Splits',2*splits);
    SpS = spectraShadow(ellipsoid(0.05*eye(2)));
    SpS = SpS & (interval(SpS)-[0.1;0]);
    SpS = SpS + [-1.8;1.7];
    plot(SpS,[1,2],'LineWidth',2,'FaceColor',colorblind('blue'),'Splits',2*splits);
    SpS = spectraShadow(ellipsoid(0.02*eye(2)));
    SpS = SpS & (interval(SpS)-[0.1;0]);
    SpS = SpS + [-1.72;1.7];
    plot(SpS,[1,2],'LineWidth',3,'FaceColor','w','Splits',2*splits);
    SpS = spectraShadow(ellipsoid(0.015*eye(2)));
    SpS = convHull(SpS, (interval(SpS)-[0.1;0]));
    SpS = SpS + [-0.9;1.8];
    plot(SpS,[1,2],'LineWidth',2,'FaceColor',colorblind('blue'),'Splits',2*splits);
    SpS = spectraShadow(ellipsoid(0.015*eye(2)));
    SpS = convHull(SpS, (interval(SpS)-[0.1;0]));
    SpS = 0.4*SpS + [-0.9;1.8];
    plot(SpS,[1,2],'LineWidth',2,'FaceColor','w','Splits',2*splits);
    
    % close figure
    close;
catch ME
    close;
    rethrow(ME)
end

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
