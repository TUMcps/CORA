function pgon = polygon(pZ,varargin)
% polygon - creates a polygon enclosure of a 2-dimensional polyZonotope
%
% Syntax:
%    pgon = polygon(pZ)
%    pgon = polygon(pZ,splits)
%
% Inputs:
%    pZ - polyZonotope object
%    splits - number of splits for refinement (optional)
%
% Outputs:
%    pgon - polygon object
%
% Example: 
%    pZ = polyZonotope([0;0],[2 0 1;0 2 1],[0;0],[1 0 3;0 1 1]);
%    pgon = polygon(pZ,8);
%
%    plot(pgon);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plot

% Authors:       Niklas Kochdumper
% Written:       08-April-2020
% Last update:   29-June-2024 (TL, bug fix during numeric issue check)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
splits = setDefaultValues({8},varargin);

% check input arguments
inputArgsCheck({{pZ,'att','polyZonotope'};
                {splits,'att','numeric','nonnan'}});

% quick return if pZ is empty
if representsa_(pZ,'emptySet',eps)
    pgon = polygon(); return
end

% quick return if pZ represents a zonotope
[res,Z] = representsa(pZ,'zonotope');
if res
    pgon = polygon(Z); return;
end

% check if polynomial zonotope is two dimensional
if dim(pZ) ~= 2
    throw(CORAerror('CORA:noExactAlg',pZ,...
        'Method "polygon" is only applicable for 2D polyZonotopes!'));
end

% split the polynomial zonotope multiple times to obtain a better 
% over-approximation of the real shape

% turn off warning from polygon
warOrig = warning;
warning('off','all');

% init with no split
pZsplit{1} = pZ;
[polyUnion,V] = aux_getPolygon(pZ);
polyPrev = polyUnion;
Vlist{1} = V;
polylist{1} = polyUnion;

% loop over number of splits
for i=1:splits
    
    % construct new polygon
    scounter = 0;
    pZnew = cell(1,2*length(pZsplit));
    Vnew = cell(1,2*length(pZsplit));
    polysnew = cell(1,2*length(pZsplit));
    polyUnion = [];
    
    % loop over each subset
    for j = 1:length(pZsplit)
        
        % only split the subset if it is part of the boundary
        if any(sum(ismembertol(Vlist{j}',polyPrev.Vertices,1e-12),2) == 2)
            % split set if on boundary

            % split the polynomial zonotope
            res = splitLongestGen(pZsplit{j});
            
            % compute the corresponding polygons
            [poly1,V1] = aux_getPolygon(res{1});
            scounter = scounter + 1;
            pZnew{scounter} = res{1};
            Vnew{scounter} = V1;
            polysnew{scounter} = poly1;
            
            [poly2,V2] = aux_getPolygon(res{2});
            scounter = scounter + 1;
            pZnew{scounter} = res{2};
            Vnew{scounter} = V2;
            polysnew{scounter} = poly2;
            
            % unite the polygons
            poly = union(poly1,poly2);
            
        else
            % no splitting if not on boundary

            % reuse previous result
            scounter = scounter + 1;
            pZnew{scounter} = pZsplit{j};
            Vnew{scounter} = Vlist{j};
            poly = polylist{j};
            polysnew{scounter} = poly;
        end
        
        % unite with previous polygons
        polyUnion = aux_unitePolygons(polyUnion,poly);
    end
    
    % update lists
    pZsplit = pZnew(1:scounter);
    Vlist = Vnew(1:scounter);
    polylist = polysnew(1:scounter);
    polyPrev = polyUnion;
end

% restore warning
warning(warOrig);

% obtain final polygon
pgon = polygon(polyUnion);

end


% Auxiliary functions -----------------------------------------------------

function [poly,V] = aux_getPolygon(pZ)
    % enclose polynomial zonotope with a polygon

    % zonotope over-approximation
    Z = zonotope(pZ);

    % calculate vertices of zonotope
    V = vertices_(Z);

    % transform to 2D polyshape (zonotope/vertices are already 'simple')
    poly = polyshape(V','Simplify',false);
end

function polyUnion = aux_unitePolygons(poly1,poly2)
% unites two polygons

if isempty(poly1)
    polyUnion = poly2;
elseif isempty(poly2)
    polyUnion = poly1;
else
    polyUnion = union(poly1,poly2);
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

% ------------------------------ END OF CODE ------------------------------
