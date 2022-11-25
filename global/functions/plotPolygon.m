function han = plotPolygon(V,varargin)
% plotPolygon - plot a polygon defined by its vertices
%
% Syntax:  
%    han = plotPolygon(V,varargin)
%
% Inputs:
%    V - matrix storing the polygon vertices
%    varargin - plot settings specified as linespec and name-value pairs
%
% Outputs:
%    han - handle of graphics object
%
% Example: 
%    zono = zonotope.generateRandom(2);
%    V = vertices(zono);
%
%    plotPolygon(V,'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Author:       Niklas Kochdumper
% Written:      05-May-2020
% Last update:  15-July-2020 (MW, merge with plotFilledPolygon)
% Last revision:---

%------------- BEGIN CODE --------------

% default
linespec = 'b';
filled = false;
NVpairs = {};

% read plot options
if ~isempty(varargin)
    [linespec,NVpairs] = readPlotOptions(varargin);
    [NVpairs,filled] = readNameValuePair(NVpairs,'Filled','islogical');
end

% check if vertex array is not empty
if ~isempty(V)

    % compute convex hull for more than two vertices
    if size(V,2)>2

        try
            ind = convhull(V(1,:),V(2,:));
            V = V(:,ind);
        catch
            error('Plotting the set failed');
        end
    end

    % plot the constrained zonotope
    if filled
        han = fill(V(1,:), V(2,:), linespec, NVpairs{:});
    else
        han = plot(V(1,:), V(2,:), linespec, NVpairs{:});
    end
end

%------------- END OF CODE --------------