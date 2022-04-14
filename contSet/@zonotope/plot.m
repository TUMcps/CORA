function han = plot(Z,varargin)
% plot - plots 2-dimensional projection of a zonotope
%
% Syntax:  
%    h = plot(Z) plots the zonotope Z for the first two dimensions
%    h = plot(Z,dims) plots the zonotope Z for the two dimensions i,j:
%        "dims=[i,j]" and returns handle to line-plot object
%    h = plot(Z,dims,'Color','red',...) adds the standard plotting preferences
%
% Inputs:
%    Z - zonotope object
%    dims - (optional) dimensions onto which the zonotope is projected
%    type - (optional) plot settings (LineSpec and name-value pairs)
%
% Outputs:
%    han - handle of graphics object
%
% Example: 
%    Z=zonotope([1 1 0; 0 0 1]);
%    plot(Z)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polygon

% Author:       Matthias Althoff
% Written:      27-July-2016
% Last update:  14-July-2020 (merge with plotFilled)
% Last revision:---

%------------- BEGIN CODE --------------

% default values
dims = [1,2];
linespec = 'b';
filled = false;
height = [];
NVpairs = {};

% if two arguments are passed    
if nargin==2
    dims=varargin{1};
    
% if three or more arguments are passed
elseif nargin>=3
    dims = varargin{1};
    % parse plot options
    [linespec,NVpairs] = readPlotOptions(varargin(2:end));
    [NVpairs,filled] = readNameValuePair(NVpairs,'Filled','islogical',filled);
    [NVpairs,height] = readNameValuePair(NVpairs,'Height','isscalar');
end

% check dimension
if length(dims) < 2
    error('At least 2 dimensions have to be specified!');
elseif length(dims) > 3
    error('Only up to 3 dimensions can be plotted!');
end

% project zonotope
Z = project(Z,dims);

% 2D or 3D plot
if length(dims) == 2

    % convert zonotope to polygon
    p = polygon(Z);

    % plot and output the handle
    if filled
        if isempty(height) % no 3D plot
            han = fill(p(1,:),p(2,:),linespec,NVpairs{:});
        else
            zCoordinates = height*ones(length(p(1,:)),1); 
            han = fill3(p(1,:),p(2,:),zCoordinates,linespec,NVpairs{:}); 
        end
    else   
        if isempty(height) % no 3D plot
            han = plot(p(1,:),p(2,:),linespec,NVpairs{:});
        else
            zCoordinates = height*ones(length(p(1,:)),1); 
            han = plot3(p(1,:),p(2,:),zCoordinates,linespec,NVpairs{:}); 
        end
    end

else
    
    % compute vertices
    V = vertices(Z);
    
    % generate 3D plot
    if ~isempty(V)
        han = plotPolytope3D(V(dims,:),linespec,NVpairs{:},'Filled',filled);
    else
        han = [];
    end
end

%------------- END OF CODE --------------