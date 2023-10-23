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

resvec = [];
load cases.mat E_c

for i=1:length(E_c)
    E1 = E_c{i}.E1;
    Ed1 = E_c{i}.Ed1;
    E0 = E_c{i}.E0;
    
    resvec(end+1) = aux_tryPlot(E1) && aux_tryPlot(Ed1) && aux_tryPlot(E0);

    % check if plotted correctly
    figure;
    ax = gca();
    colorOrder = ax.ColorOrder;
    
    % plot set
    plot(E1,[1,2]);
    V = [-3.0497 -3.0597 0.6514 1.5522 -2.6948 -3.0279 -3.0497;
        -1.8645 -1.8570 7.5677 9.0267 -1.4353 -1.8649 -1.8645];
    % check points
    resvec(end+1) = compareMatrices(V, [ax.Children(1).XData;ax.Children(1).YData],1e-4,'subset',true);
    % test color
    resvec(end+1) = isequal(colorOrder(1,:), ax.Children(1).Color);
    close;

    if ~resvec(end)
        break;
    end 
end

% gather results
res = all(resvec);

end


% Auxiliary functions -----------------------------------------------------

function res = aux_tryPlot(E)

resvec = [];
try
    % try all variations in plotting
    h = figure;
    
    % one argument: object
    plot(E);
    resvec(end+1) = true;

    % two arguments: object and dimension
    plot(E,1);
    plot(E,[1,2]);
    resvec(end+1) = true;
    
    % three arguments: object, dimensions, linespec
    plot(E,[1,2],'r+');
    resvec(end+1) = true;
    
    % three arguments: object, dimensions, NVpairs
    plot(E,[1,2],'LineWidth',2);
    plot(E,[1,2],'Color',[.6 .6 .6],'LineWidth',2);
    
    % four arguments: object, dimensions, linespec, NVpairs
    plot(E,[1,2],'FaceColor','r','LineWidth',2);
    if isFullDim(project(E,[1,2]))
        plot(E,[1,2],'FaceColor','r','LineWidth',2,'EdgeColor',[.6 .6 .6]);
    end
    resvec(end+1) = true;
    
    % close figure
    close(h);
catch ME
    close;
    resvec(end+1) = false;
end

% gather results
res = all(resvec);

end

% ------------------------------ END OF CODE ------------------------------
