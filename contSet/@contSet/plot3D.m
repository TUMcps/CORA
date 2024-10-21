function han = plot3D(S,varargin)
% plot3D - plots a 3D projection of a contSet
%
% Syntax:
%    han = plot3D(S)
%    han = plot3D(S,NVpairs)
%
% Inputs:
%    S - projected contSet object
%    NVpairsPlot - (optional) plot settings (LineSpec and Name-Value pairs)
%    NVpairsVertices - (optional) vertex computation settings
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
[NVpairsPlot,NVpairsVertices] = setDefaultValues({{},{}},varargin);

% compute vertices
V = vertices(S,NVpairsVertices{:});

% plot vertices
han = plotPolytope3D(V,NVpairsPlot{:});

end

% ------------------------------ END OF CODE ------------------------------
