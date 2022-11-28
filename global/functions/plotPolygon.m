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
%    han - handle to the graphics object
%
% Example: 
%    zono = zonotope.generateRandom('Dimension',2);
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
NVpairs = {'Color',colorblind('b')};

% read plot options
if ~isempty(varargin)
    NVpairs = readPlotOptions(varargin);
end
% readout 'FaceColor' to decide plot/fill call where necessary
[~,facecolor] = readNameValuePair(NVpairs,'FaceColor');

% check if vertex array is not empty
if ~isempty(V)

    % compute convex hull for more than two vertices
    if size(V,2)>2

        try
            ind = convhull(V(1,:),V(2,:));
            V = V(:,ind);
        catch
            throw(CORAerror('CORA:specialError','Plotting the set failed'));
        end
    end

    % plot the constrained zonotope
    if isempty(facecolor) || strcmp(facecolor,'none') 
        % axis are reset even if they were set manually. Correct afterwards
        xmode = xlim("mode");
        x = xlim;
        ymode = ylim("mode");
        y = ylim;

        han = plot(V(1,:), V(2,:), NVpairs{:});

        if xmode == "manual"
            xlim(x);
        end
        if ymode == "manual"
            ylim(y);
        end
    else
        han = fill(V(1,:), V(2,:), facecolor, NVpairs{:});
    end
    
end

%------------- END OF CODE --------------