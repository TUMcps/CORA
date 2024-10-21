function [xLim,yLim,zLim] = getUnboundedAxisLimits(V)
% getUnboundedAxisLimits - returns the axis limits to plot unbounded sets
%    Background: just using xlim/ylim is not sufficient as plotting at
%    axis limits might auto-adjust them.
%
% Syntax:
%    [xLim,yLim] = getUnboundedAxisLimits()
%    [xLim,yLim] = getUnboundedAxisLimits(V)
%    [xLim,yLim,zLim] = getUnboundedAxisLimits(V)
%
% Inputs:
%    V - vertices to be plotted
%
% Outputs:
%    xLim - x-axis limits
%    yLim - y-axis limits
%    zLim - z-axis limits
%

% Authors:       Tobias Ladner
% Written:       26-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% idea: plot at axis limits and see if they are auto adjusted

% parse input
plotDim = numel(axis)/2;
if nargin < 1
    V = zeros(plotDim,0);
end
plotDim = max(plotDim,size(V,1));
if size(V,1) < 1 || size(V,1) > 3
    throw(CORAerror('CORA:specialError','Specified vertices have to be 2- or 3-dimensional.'))
end

% append boundary
if plotDim == 2
    if size(V,1) == 1
        V = [V, xlim];
        V = [V;zeros(1,size(V,2))];
    else 
        V = [V, [xlim;ylim]];
    end
else % plotDim == 3
    if size(V,1) == 2 % && plotDim == 3
        V = [V;zeros(1,size(V,2))];
        V = [V, [xlim;ylim;zlim]];
    else % size(V,1) == 3 && plotDim == 3
        V = [V, [xlim;ylim;zlim]];
    end
end

% set mode to 'auto'
xmode = xlim("mode");
xlim('auto')
ymode = ylim("mode");
ylim('auto')
zmode = zlim("mode");
zlimOld = zlim();
zlim('auto');

% get color order
oldColorIndex = get(gca(),'ColorOrderIndex');

if size(V,1) == 2
    % plot at 2d boundary
    han = plot(V(1,:), V(2,:), 'HandleVisibility','off');
else
    % plot at 3d boundary
    han = plot3(V(1,:), V(2,:), V(3,:), 'HandleVisibility','off');
end

% read axis limits
xLim = xlim;
yLim = ylim;
zLim = zlim;

% delete plot
delete(han)
updateColorIndex(oldColorIndex-1)

% restore axis mode
xlim(xmode);
ylim(ymode);
zlim(zlimOld);
zlim(zmode);

end


% ------------------------------ END OF CODE ------------------------------
