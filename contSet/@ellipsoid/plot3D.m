function han = plot3D(E,varargin)
% plot3D - plots a 3D projection of a contSet
%
% Syntax:
%    han = plot3D(S)
%    han = plot3D(S,dims)
%    han = plot3D(S,dims,NVpairs)
%
% Inputs:
%    S - contSet object
%    dims - (optional) dimensions for projection
%    NVpairsPlot - (optional) plot settings (LineSpec and Name-Value pairs)
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
[NVpairsPlot] = setDefaultValues({{}},varargin);

% outer-approximate using zonotope
Z = zonotope(E,'outer:norm',100);

% plot zonotope enclosure
han = plot3D(Z,NVpairsPlot);

end

% ------------------------------ END OF CODE ------------------------------
