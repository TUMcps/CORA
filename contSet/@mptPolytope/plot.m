function han = plot(obj,varargin)
% plot - Plots 2-dimensional projection of a mptPolytope
%
% Syntax:  
%    han = plot(obj)
%    han = plot(obj,dims)
%    han = plot(obj,dims,type)
%
% Inputs:
%    obj - mptPolytope object
%    dims - (optional) dimensions of the projection
%    type - (optional) plot settings (LineSpec and name-value pairs)
%
% Outputs:
%    han - handle to the plotted graphics object
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      19-October-2010
% Last update:  24-March-2015
%               17-March-2017
%               20-April-2018 (exception for empty sets)
% Last revision:---

%------------- BEGIN CODE --------------

% default settings
dims = [1,2];
plotOptions{1} = 'b';

% parse input arguments
if nargin >= 2 && ~isempty(varargin{1})
	dims = varargin{1}; 
end

if nargin >= 3
	plotOptions = varargin(2:end); 
end

% check dimension
if length(dims) < 2
    error('At least 2 dimensions have to be specified!');
elseif length(dims) > 3
    error('Only up to 3 dimensions can be plotted!');
end
 
% check if polyhedron is bounded
if ~isBounded(obj.P)
    
    % get size of current plot
    xLim = get(gca,'Xlim');
    yLim = get(gca,'Ylim');
    
    % intersect with the current polytope
    if length(dims) == 2
        obj = project(obj,dims) & interval([xLim(1);yLim(1)], ...
                                           [xLim(2);yLim(2)]);
        dims = [1,2];
    else
        zLim = get(gca,'Zlim');
        obj = project(obj,dims) & interval([xLim(1);yLim(1);zLim(1)], ...
                                           [xLim(2);yLim(2),zLim(2)]);
        dims = [1,2,3];
    end
end
    
% compute vertices
V = vertices(obj);

% plot projected vertices
if length(dims) == 2
    han = plotPolygon(V(dims,:),plotOptions{:});
elseif length(dims) == 3
    han = plotPolytope3D(V(dims,:),plotOptions{:});
end

%------------- END OF CODE --------------