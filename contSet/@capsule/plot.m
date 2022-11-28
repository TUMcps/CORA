function han = plot(C,varargin)
% plot - plots a projection of a capsule
%
% Syntax:  
%    han = plot(C)
%    han = plot(C,dims)
%    han = plot(C,dims,type)
%
% Inputs:
%    C - capsule object
%    dims - (optional) dimensions for projection
%    type - (optional) plot settings (LineSpec and Name-Value pairs)
%
% Outputs:
%    han - handle to the graphics object
%
% Example: 
%    C = capsule([1; 1; 0], [0; 1; 1], 0.5);
%    plot(C)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polygon

% Author:       Matthias Althoff
% Written:      04-March-2019
% Last update:  25-May-2022 (TL: 1D Plotting)
% Last revision:---

%------------- BEGIN CODE --------------

% parse input arguments
dims = setDefaultValues({[1,2]},varargin{:});
% check input arguments
inputArgsCheck({{C,'att','capsule','nonempty'};
                {dims,'att','numeric',{'nonempty','vector','integer','positive'}} });

% parse plot options
NVpairs = readPlotOptions(varargin(2:end));

% readout 'FaceColor' value to decide plot/fill call where necessary
[~,facecolor] = readNameValuePair(NVpairs,'FaceColor');

% check dimension
if length(dims) < 1
    throw(CORAerror('CORA:plotProperties',1));
elseif length(dims) > 3
    throw(CORAerror('CORA:plotProperties',3));
end


% 1D, 2D vs 3D plot
if length(dims) == 1
    han = plot(interval(C), 1, NVpairs{:});

elseif length(dims) == 2

    % project capsule
    C = project(C,dims);

    % underapproximate capsule by polygon
    p = polygon(C);

    % plot and output the handle
    if isempty(facecolor) || strcmp(facecolor,'none')
        han = plot(p(1,:),p(2,:),NVpairs{:});
    else
        han = fill(p(1,:),p(2,:),facecolor,NVpairs{:});
    end

else

    % generate 3D plot
    if isempty(C.g)         % generator is empty
        E = ellipsoid(eye(3)*C.r^2,C.c);
        Z = zonotope(E,100,'outer:norm');
        han = plot(Z,[1,2,3],NVpairs{:});
        
    else                    % generator is not empty
        
        % generate ellipsoids for spheres
        E1 = ellipsoid(eye(3)*C.r^2,C.c + C.g);
        E2 = ellipsoid(eye(3)*C.r^2,C.c - C.g);
        
        % enclose ellipsoid with a zonotope
        Z1 = zonotope(E1,100,'outer:norm');
        Z2 = zonotope(E2,100,'outer:norm');
        
        % compute vertices of the zonotopes
        V1 = vertices(Z1);
        V2 = vertices(Z2);
        
        % plot the convex hull
        han = plotPolytope3D([V1,V2],NVpairs{:});
    end
end

if nargout == 0
    clear han;
end

%------------- END OF CODE --------------