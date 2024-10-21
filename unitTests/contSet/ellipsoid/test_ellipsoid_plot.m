function res = test_ellipsoid_plot
% test_ellipsoid_plot - unit test function of plot
%
% Syntax:
%    res = test_ellipsoid_plot
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

% Authors:       Victor Gassmann
% Written:       27-July-2021
% Last update:   09-May-2023 (TL, added plotted point checks)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init cases
E1 = ellipsoid([ 5.4387811500952807 12.4977183618314545 ; 12.4977183618314545 29.6662117284481646 ], [ -0.7445068341257537 ; 3.5800647524843665 ], 0.000001);
Ed1 = ellipsoid([ 4.2533342807136076 0.6346400221575308 ; 0.6346400221575309 0.0946946398147988 ], [ -2.4653656883489115 ; 0.2717868749873985 ], 0.000001);
E0 = ellipsoid([ 0.0000000000000000 0.0000000000000000 ; 0.0000000000000000 0.0000000000000000 ], [ 1.0986933635979599 ; -1.9884387759871638 ], 0.000001);
    
assert(aux_tryPlot(E1))
assert(aux_tryPlot(Ed1))
assert(aux_tryPlot(E0));

% check if plotted correctly
figure;
ax = gca();
colorOrder = ax.ColorOrder;

% plot set
plot(E1,[1,2]);
V = [-3.0497 -3.0597 0.6514 1.5522 -2.6948 -3.0279 -3.0497;
    -1.8645 -1.8570 7.5677 9.0267 -1.4353 -1.8649 -1.8645];
% check points
assert(compareMatrices(V, [ax.Children(1).XData;ax.Children(1).YData],1e-4,'subset',true));
% test color
assert(isequal(colorOrder(1,:), ax.Children(1).Color));
close;

% gather results
res = true;

end


% Auxiliary functions -----------------------------------------------------

function res = aux_tryPlot(E)

try
    % try all variations in plotting
    h = figure;
    
    % one argument: object
    plot(E);

    % two arguments: object and dimension
    plot(E,1);
    plot(E,[1,2]);
    
    % three arguments: object, dimensions, linespec
    plot(E,[1,2],'r+');
    
    % three arguments: object, dimensions, NVpairs
    plot(E,[1,2],'LineWidth',2);
    plot(E,[1,2],'Color',[.6 .6 .6],'LineWidth',2);
    
    % four arguments: object, dimensions, linespec, NVpairs
    plot(E,[1,2],'FaceColor','r','LineWidth',2);
    if isFullDim(project(E,[1,2]))
        plot(E,[1,2],'FaceColor','r','LineWidth',2,'EdgeColor',[.6 .6 .6]);
    end
    
    % close figure
    close(h);
catch ME
    close;
    rethrow(ME)
end

% gather results
res = true;

end

% ------------------------------ END OF CODE ------------------------------
