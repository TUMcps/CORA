function han = plot(Z,varargin)
% plot - plots a projection of a zonotope
%
% Syntax:  
%    han = plot(Z)
%    han = plot(Z,dims)
%    han = plot(Z,dims,type)
%
% Inputs:
%    Z - zonotope object
%    dims - (optional) dimensions for projection
%    type - (optional) plot settings (LineSpec and Name-Value pairs),
%        including added pairs:
%          'Height', <height> height of z-coordinate
%
% Outputs:
%    han - handle to the graphics object
%
% Example: 
%    Z = zonotope([1;0],[0.4 -1 0.4; -1 0.3 0]);
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
%               25-May-2022 (TL: 1D Plotting)
% Last revision:---

%------------- BEGIN CODE --------------

% default values
dims = setDefaultValues({[1,2]},varargin);

% parse plot options
NVpairs = readPlotOptions(varargin(2:end));
[NVpairs,height] = readNameValuePair(NVpairs,'Height','isscalar');
% readout 'FaceColor' to decide plot/fill call where necessary
[NVpairs,facecolor] = readNameValuePair(NVpairs,'FaceColor');

% check dimension
if length(dims) < 1
    throw(CORAerror('CORA:plotProperties',1));
elseif length(dims) > 3
    throw(CORAerror('CORA:plotProperties',3));
end

% project zonotope
Z = project(Z,dims);

% 1D, 2D or 3D plot
if length(dims) == 1
    V = vertices(Z);
    han = plot(V,[0,0],NVpairs{:});

elseif length(dims) == 2

    % convert zonotope to polygon
    p = polygon(Z);

    % plot and output the handle
    if isempty(facecolor) || strcmp(facecolor,'none')
        if isempty(height) % no 3D plot
            han = plot(p(1,:),p(2,:),NVpairs{:});
        else
            zCoordinates = height*ones(length(p(1,:)),1); 
            han = plot3(p(1,:),p(2,:),zCoordinates,NVpairs{:}); 
        end
    else   
        if isempty(height) % no 3D plot
            han = fill(p(1,:),p(2,:),facecolor,NVpairs{:});
        else
            zCoordinates = height*ones(length(p(1,:)),1); 
            han = fill3(p(1,:),p(2,:),zCoordinates,facecolor,NVpairs{:}); 
        end
    end

else
    
    % compute vertices
    V = vertices(Z);
    
    % generate 3D plot
    if ~isempty(V)
        han = plotPolytope3D(V(dims,:),NVpairs{:});
    else
        han = [];
    end
end

if nargout == 0
    clear han;
end

%------------- END OF CODE --------------