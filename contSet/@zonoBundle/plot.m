function han = plot(obj,varargin)
% plot - plots 2-dimensional over-approximative projection of a zonotope bundle
%
% Syntax:  
%    han = plot(obj,dims,type)
%
% Inputs:
%    obj - zonoBundle object
%    dims - dimensions that should be projected (optional) 
%    type - plot options (LineSpec and Name-Value pairs)
%
% Outputs:
%    han - handle for the graphics object
%
% Example: 
%    Z{1} = zonotope([0 1 1 -1;0 0 2 1]);
%    Z{2} = zonotope([2 1 1 -1;1 0 -1 -1]);
%    zB = zonoBundle(Z);
%
%    plot(zB,[1,2],'r','LineWidth',2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      09-November-2010 
% Last update:  13-February-2012
%               19-October-2015 (NK, accelerate by using polyshape class)
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

% 2D vs. 3D plot
if length(dims) == 2

    w = warning;
    warning('off','all');

    % compute polytopes
    for i = 1:obj.parallelSets

        % delete zero generators
        Z = deleteZeros(obj.Z{i});

        % project zonotope
        Z = project(Z,dims);

        % convert to polyshape (Matlab build in class)
        temp = polygon(Z);
        V = temp(:,2:end);
        P{i} = polyshape(V(1,:),V(2,:));
    end

    % intersect polytopes
    Pint = P{1};
    for i = 2:obj.parallelSets
        Pint = intersect(Pint,P{i});
    end

    % get vertices
    V = Pint.Vertices';

    % reset warning state to previous setting
    warning(w);

    % plot polytope
    han = plotPolygon(V,plotOptions{:});

else
    
    % project to plotted dimensions
    obj = project(obj,dims);
    
    % compute vertices
    V = vertices(obj);
    
    % plot 3D polytope
    han = plotPolytope3D(V(dims,:),plotOptions{:}); 
end

%------------- END OF CODE --------------