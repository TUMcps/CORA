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

if isBounded(S)
    % compute (outer-approximative) vertices via polygon
    pgon = polygon(S,NVpairsPolygon{:});
    
    % read vertices
    V = vertices_(pgon);
else
    % compute vertices directly
    V = vertices_(S);
end

% add first vertex at the end to close the polygon
if size(V,2) > 1 && ~any(isnan(V),"all")
    V = [V,V(:,1)];
end

% plot vertices
han = plotPolygon(V,NVpairsPlot{:});

end

% ------------------------------ END OF CODE ------------------------------
