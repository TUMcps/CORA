function han = plot2D(S,varargin)
% plot2D - plots a 2D projection of a contSet
%
% Syntax:
%    han = plot2D(S)
%    han = plot2D(S,NVpairs)
%
% Inputs:
%    S - projected contSet object
%    NVpairsPlot - (optional) plot settings (LineSpec and Name-Value pairs)
%    NVpairsPolygon - (optional) polygon computation settings
%
% Outputs:
%    han - handle to the graphics object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plot

% Authors:       Tobias Ladner
% Written:       14-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse additional parameters
[NVpairsPlot,NVpairsPolygon] = setDefaultValues({{},{}},varargin);

if ~isBounded(S) || isa(S,'interval')
    % compute vertices directly
    % (+ direct computation for (degenerate) intervals to speed up reachSet/plotOverTime+timepoint)
    % (TODO: find general solution to check if direct computation should be used)
    V = vertices_(S);
else
    % compute (outer-approximative) vertices via polygon
    % (obtains tighter results for some set representations due to splits
    %  or a result at all as e.g. ellipsoids/vertices is not feasible for 2D)
    pgon = polygon(S,NVpairsPolygon{:});
    
    % read vertices
    V = vertices_(pgon);
end

% add first vertex at the end to close the polygon
if size(V,2) > 1 && ~any(isnan(V),"all")
    V = [V,V(:,1)];
end

% plot vertices
han = plotPolygon(V,NVpairsPlot{:});

end

% ------------------------------ END OF CODE ------------------------------
