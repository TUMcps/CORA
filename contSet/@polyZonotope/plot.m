function han = plot(pZ,varargin)
% plot - Plots 2-dimensional over-approximation of a polynomial zonotope
%
% Syntax:  
%    han = plot(pZ)
%    han = plot(pZ,dims,linespec)
%    han = plot(pZ,dims,linespec,'Splits',splits)
%
% Inputs:
%    pZ - polyZonotope object
%    dims - dimensions that should be projected
%    linespec - (optional) LineSpec properties
%    splits - (optional) number of splits for refinement
%    type - (optional) name-value pairs
%
% Outputs:
%    han - handle for the resulting graphics object
%
% Example: 
%    pZ = polyZonotope([0;0],[2 0 1;0 2 1],[0;0],[1 0 3;0 1 1]);
%     
%    figure; hold on;
%    plotRandPoint(pZ,[1,2],100000,'.r');
%    plot(pZ,[1,2],'b','Splits',3);
%
%    figure; hold on;
%    plotRandPoint(pZ,[1,2],100000,'.r');
%    plot(pZ,[1,2],'b','Splits',10);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plotRandPoint

% Author:       Niklas Kochdumper, Mark Wetzlinger
% Written:      29-March-2018
% Last update:  23-June-2020 (MW, harmonize with other plot functions)
%               14-July-2020 (MW, merge with plotFilled)
% Last revision:---

%------------- BEGIN CODE --------------

% default values for the optional input arguments
dims = [1,2];
linespec = 'b';
splits = 8;
NVpairs = {};
filled = 0;

% parse input arguments
if nargin > 1 && ~isempty(varargin{1})
    dims = varargin{1}; 
end
if nargin > 2 && ~isempty(varargin{2})
    % read additional name-value pairs
    [linespec,NVpairs] = readPlotOptions(varargin(2:end));
    [NVpairs,filled] = readNameValuePair(NVpairs,'Filled','islogical',filled);
    [NVpairs,splits] = readNameValuePair(NVpairs,'Splits','isscalar',splits);
end

% delete all zero-generators
pZ = deleteZeros(pZ);

% split the polynomial zonotope multiple times to obtain a better 
% over-approximation of the real shape
pZsplit{1} = project(pZ,dims);

for i=1:splits
    pZnew = [];
    for j=1:length(pZsplit)
        res = splitLongestGen(pZsplit{j});
        pZnew{end+1} = res{1};
        pZnew{end+1} = res{2};
    end
    pZsplit = pZnew;
end

% over-approximate all splitted sets with zonotopes, convert them to 2D
% polytopes (Matlab built-in) and compute the union
warOrig = warning;
warning('off','all');

for i = 1:length(pZsplit)
   
    % zonotope over-approximation
    zono = zonotope(pZsplit{i});
    
    % calculate vertices of zonotope
    vert = vertices(zono);
    
    % transform to 2D polytope
    pgon = polyshape(vert(1,:),vert(2,:));
    
    % calculate union with previous sets
    if i == 1
       polyAll = pgon; 
    else
       polyAll = union(polyAll,pgon); 
    end  
end

warning(warOrig);


% (previous method) plot the polygon
% if isempty(type)
%     han = plot(polyAll,'FaceAlpha',1,'FaceColor','none','EdgeColor',linespec);
% else
%     han = plot(polyAll,'FaceAlpha',1,'FaceColor','none','EdgeColor',linespec,type{:});
% end

% add first point to end to close polytope
xVals = [polyAll.Vertices(:,1);polyAll.Vertices(1,1)];
yVals = [polyAll.Vertices(:,2);polyAll.Vertices(1,2)];
if filled
    han = fill(xVals,yVals,linespec,NVpairs{:});
else
    han = plot(xVals,yVals,linespec,NVpairs{:});
end

%------------- END OF CODE --------------