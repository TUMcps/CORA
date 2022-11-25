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
splits = 10;
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

% check dimension
if length(dims) < 2
    error('At least 2 dimensions have to be specified!');
elseif length(dims) > 3
    error('Only up to 3 dimensions can be plotted!');
end

% delete all zero-generators
pZ = deleteZeros(pZ);

% project to desired dimensions
pZ = project(pZ,dims);

% 2D vs 3D plot
if length(dims) == 2
    
    warOrig = warning;
    warning('off','all');
    
    % split the polynomial zonotope multiple times to obtain a better 
    % over-approximation of the real shape
    pZsplit{1} = pZ;
    [polyPrev,V] = getPolygon(pZ);
    Vlist{1} = V;

    for i=1:splits
        
        pZnew = [];
        Vnew = [];
        polyAll = [];
        
        for j = 1:length(pZsplit)
            
            % only split the set if it is part of the boundary
            if any(sum(ismembertol(Vlist{j}',polyPrev.Vertices,1e-12),2) == 2)
                
                % split the polynomial zonotope
                res = splitLongestGen(pZsplit{j});
                
                % compute the corresponding polygons
                [poly1,V1] = getPolygon(res{1});
                pZnew{end+1} = res{1};
                Vnew{end+1} = V1;
                
                [poly2,V2] = getPolygon(res{2});
                pZnew{end+1} = res{2};
                Vnew{end+1} = V2;
                
                % unite the polygons
                poly = union(poly1,poly2);
                
            else
                
                % compute polygon
                [poly,V] = getPolygon(pZsplit{j});
                pZnew{end+1} = pZsplit{j};
                Vnew{end+1} = V;
            end
            
            % unite with previous polygons
            if isempty(polyAll)
                polyAll = poly;
            else
                polyAll = union(polyAll,poly);
            end
        end
        
        % update lists
        pZsplit = pZnew;
        Vlist = Vnew;
        polyPrev = polyAll;
    end

    warning(warOrig);

    % add first point to end to close polygon
    xVals = [polyAll.Vertices(:,1);polyAll.Vertices(1,1)];
    yVals = [polyAll.Vertices(:,2);polyAll.Vertices(1,2)];
    
    % plot the polygon
    if filled
        han = fill(xVals,yVals,linespec,NVpairs{:});
    else
        han = plot(xVals,yVals,linespec,NVpairs{:});
    end

else
    
    % split the polynomial zonotope multiple times to obtain a better 
    % over-approximation of the real shape
    pZsplit{1} = pZ;

    for i=1:splits
        pZnew = [];
        for j=1:length(pZsplit)
            res = splitLongestGen(pZsplit{j});
            pZnew{end+1} = res{1};
            pZnew{end+1} = res{2};
        end
        pZsplit = pZnew;
    end
    
    % loop over all parallel sets
    hold on;
    for i = 1:length(pZsplit{i})
        if filled
            han = plot(zonotope(pZsplit{i}),[1,2,3],linespec, ...
                       NVpairs{:},'Filled',true); 
        else
            han = plot(zonotope(pZsplit{i}),[1,2,3],linespec,NVpairs{:});  
        end
    end
end
end


% Auxiliary Functions -----------------------------------------------------

function [poly,V] = getPolygon(pZ)
% enclose polynomial zonotope with a polygon

    % zonotope over-approximation
    zono = zonotope(pZ);

    % calculate vertices of zonotope
    V = vertices(zono);

    % transform to 2D polytope
    poly = polyshape(V(1,:),V(2,:));
end

%------------- END OF CODE --------------