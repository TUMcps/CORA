function han = plot(pZ,varargin)
% plot - plots an over-approximative projection of a polynomial zonotope
%
% Syntax:  
%    han = plot(pZ)
%    han = plot(pZ,dims)
%    han = plot(pZ,dims,type)
%
% Inputs:
%    pZ - polyZonotope object
%    dims - (optional) dimensions for projection
%    type - (optional) plot settings (LineSpec and Name-Value pairs)
%           additional Name-Value pairs:
%               <'Splits',splits> - number of splits for refinement
%
% Outputs:
%    han - handle to the graphics object
%
% Example: 
%    pZ = polyZonotope([0;0],[2 0 1;0 2 1],[0;0],[1 0 3;0 1 1]);
%     
%    figure; hold on;
%    plotRandPoint(pZ);
%    plot(pZ,[1,2],'b','Splits',4);
%
%    figure; hold on;
%    plotRandPoint(pZ);
%    plot(pZ,[1,2],'b','Splits',10);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plotRandPoint

% Author:       Niklas Kochdumper, Mark Wetzlinger, Tobias Ladner
% Written:      29-March-2018
% Last update:  23-June-2020 (MW, harmonize with other plot functions)
%               14-July-2020 (MW, merge with plotFilled)
%               25-May-2022 (TL: 1D Plotting)
%               23-February-2023 (TL: enlarge polygon if 2 regions)
%               23-March-2023 (TL: bugfix Splits=0, clean up)
%               05-April-2023 (TL: clean up using plotPolygon)
% Last revision:---

%------------- BEGIN CODE --------------

% parse input

% default values for the optional input arguments
dims = setDefaultValues({[1,2]},varargin);

% read additional name-value pairs
NVpairs = readPlotOptions(varargin(2:end));
% read out 'Splits', default value given
[NVpairs,splits] = readNameValuePair(NVpairs,'Splits','isscalar',10);

% check input arguments
inputArgsCheck({{pZ,'att','polyZonotope'};
    {dims,'att','numeric',{'nonempty','integer','positive','vector'}}; ...
    {splits,'att','numeric',{'scalar','integer','nonnegative'}} ...
});

% check dimension
if length(dims) < 1
    throw(CORAerror('CORA:plotProperties',1));
elseif length(dims) > 3
    throw(CORAerror('CORA:plotProperties',3));
end

% delete all zero-generators
pZ = deleteZeros(pZ);

% project to desired dimensions
pZ = project(pZ,dims);
dims = 1:dim(pZ);
pZ = compact(pZ);

% check if is zonotope
if all(sum(pZ.expMat > 0, 2) == 1)
    han = plot(zonotope(pZ), dims, NVpairs{:});

% 1D, 2D vs 3D plot
elseif length(dims) == 1
    han = plot(interval(pZ, 'split'), 1, NVpairs{:});

elseif length(dims) == 2
    
    warOrig = warning;
    warning('off','all');
    
    % split the polynomial zonotope multiple times to obtain a better 
    % over-approximation of the real shape

    % init with no split
    pZsplit{1} = pZ;
    [polyUnion,V] = aux_getPolygon(pZ);
    polyPrev = polyUnion;
    Vlist{1} = V;

    for i=1:splits
        
        % construct new polygon
        scounter = 0;
        pZnew = cell(length(pZsplit), 1);
        Vnew = cell(length(pZsplit), 1);
        polyUnion = [];
        
        for j = 1:length(pZsplit)
            
            % only split the set if it is part of the boundary
            if any(sum(ismembertol(Vlist{j}',polyPrev.Vertices,1e-12),2) == 2)
                % split set if on boundary

                % split the polynomial zonotope
                res = splitLongestGen(pZsplit{j});
                
                % compute the corresponding polygons
                [poly1,V1] = aux_getPolygon(res{1});
                scounter = scounter + 1;
                pZnew{scounter} = res{1};
                Vnew{scounter} = V1;
                
                [poly2,V2] = aux_getPolygon(res{2});
                scounter = scounter + 1;
                pZnew{scounter} = res{2};
                Vnew{scounter} = V2;
                
                % unite the polygons
                poly = union(poly1,poly2);
                
            else
                % no splitting if not on boundary

                % compute polygon
                [poly,V] = aux_getPolygon(pZsplit{j});
                scounter = scounter + 1;
                pZnew{scounter} = pZsplit{j};
                Vnew{scounter} = V;

            end
            
            % unite with previous polygons
            if isempty(polyUnion)
                polyUnion = poly;
            else
                polyUnion = union(polyUnion,poly);
            end

            if polyUnion.NumRegions >= 2
                % might be due to numeric instability
                % enlargen polygon slightly
                polyUnion_ = polybuffer(polyUnion, 1e-8);
                if polyUnion_.NumRegions < polyUnion.NumRegions
                    polyUnion = polyUnion_;
                end
            end
        end
        
        % update lists
        pZsplit = pZnew(1:scounter);
        Vlist = Vnew(1:scounter);
        polyPrev = polyUnion;
    end

    warning(warOrig);

    if size(polyUnion.Vertices,1) > 0
        % add first point to end to close polygon
        xVals = [polyUnion.Vertices(:,1);polyUnion.Vertices(1,1)];
        yVals = [polyUnion.Vertices(:,2);polyUnion.Vertices(1,2)];
    else
        % if e.g. all generators are collinear
        V = cell2mat(Vlist);
        xVals = V(1,:);
        yVals = V(2,:);
    end
    
    % plot the polygon
    han = plotPolygon([xVals,yVals]', NVpairs{:});

else
    % plot 3d
    
    % split the polynomial zonotope multiple times to obtain a better 
    % over-approximation of the real shape
    pZsplit{1} = pZ;

    for i=1:splits
        pZnew = cell(2*length(pZsplit));
        for j=1:length(pZsplit)
            res = splitLongestGen(pZsplit{j});
            pZnew{2*j-1} = res{1};
            pZnew{2*j} = res{2};
        end
        pZsplit = pZnew;
    end
    
    % loop over all parallel sets
    hold on;
    for i = 1:length(pZsplit{i})
        han = plot(zonotope(pZsplit{i}),[1,2,3],NVpairs{:}); 
    end
end

if nargout == 0
    clear han;
end

end


% Auxiliary Functions -----------------------------------------------------

function [poly,V] = aux_getPolygon(pZ)
% enclose polynomial zonotope with a polygon

    % zonotope over-approximation
    Z = zonotope(pZ);

    % calculate vertices of zonotope
    V = vertices(Z);

    % transform to 2D polytope
    poly = polyshape(V(1,:),V(2,:));
end

%------------- END OF CODE --------------