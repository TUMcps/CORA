function han = plot(C,varargin)
% plot - Plots 2-dimensional projection of a capsule
%
% Syntax:  
%    h = plot(C) plots the capsule for the first two dimensions
%    h = plot(C,dims) plots the capsule for the two dimensions i,j:
%                   "dims=[i,j]" and returns handle to line-plot object
%    h = plot(C,dims,'Color','red',...) adds the standard plotting preferences
%
% Inputs:
%    C - capsule object
%    dims - (optional) dimensions that should be projected
%    type - (optional) plot settings (LineSpec and name-value pairs)
%
% Outputs:
%    han - handle to graphics object
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
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % default
    dims=[1,2];
    linespec = 'b';
    filled = false;
    NVpairs = {};

    % if two arguments are passed    
    if nargin==2
        dims=varargin{1};

    % if three or more arguments are passed
    elseif nargin>=3
        dims=varargin{1};   
        % parse plot options
        [linespec,NVpairs] = readPlotOptions(varargin(2:end));
        [NVpairs,filled] = readNameValuePair(NVpairs,'Filled','islogical');
    end
    
    % check dimension
    if length(dims) < 2
        error('At least 2 dimensions have to be specified!');
    elseif length(dims) > 3
        error('Only up to 3 dimensions can be plotted!');
    end

    % project zonotope
    C = project(C,dims);

    % 2D vs 3D plot
    if length(dims) == 2
    
        % underapproximate capsule by polygon
        p = polygon(C);

        % plot and output the handle
        if filled
            han = fill(p(1,:),p(2,:),linespec,NVpairs{:});
        else
            han = plot(p(1,:),p(2,:),linespec,NVpairs{:});
        end
    
    else
        
        % generate 3D plot
        if isempty(C.g)         % generator is empty
            E = ellipsoid(eye(3)*C.r^2,C.c);
            Z = zonotope(E,100,'o:norm');
            han = plot(Z,[1,2,3],linespec,NVpairs{:},'Filled',filled);
        else                    % generator is not empty
            
            % generate ellipsoids for spheres
            E1 = ellipsoid(eye(3)*C.r^2,C.c + C.g);
            E2 = ellipsoid(eye(3)*C.r^2,C.c - C.g);
            
            % enclose ellipsoid with a zonotope
            Z1 = zonotope(E1,100,'o:norm');
            Z2 = zonotope(E2,100,'o:norm');
            
            % compute vertices of the zonotopes
            V1 = vertices(Z1);
            V2 = vertices(Z2);
            
            % plot the convex hull
            han = plotPolytope3D([V1,V2],linespec,NVpairs{:}, ...
                                 'Filled',filled);
        end
    end

%------------- END OF CODE --------------