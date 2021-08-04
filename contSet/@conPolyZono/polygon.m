function pgon = polygon(cPZ,varargin)
% polygon - compute polygon enclosure of a 2D constrained poly. zonotope
%
% Syntax:  
%    pgon = polygon(cPZ)
%    pgon = polygon(cPZ,splits)
%
% Inputs:
%    cPZ - conPolyZono object
%    splits - number of splits for refinement
%
% Outputs:
%    pgon - polygon object
%
% Example: 
%    c = [0;0];
%    G = [2 1; 0 1];
%    expMat = [1 0; 0 1; 0 0];
%    A = [1 1 -0.25];
%    b = 0.75;
%    expMat_ = [2 0 0; 0 2 0; 0 0 1];
%    cPZ = conPolyZono(c,G,expMat,A,b,expMat_);
%
%    pgon = polygon(cPZ,15);
%    
%    plot(pgon,[1,2],'b','Filled',true,'EdgeColor','none');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polygon, conPolyZonotope, zonotope

% Author:       Niklas Kochdumper
% Written:      19-January-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

     % parse input arguments
    splits = 10;
    
    if nargin > 1 
       splits = varargin{1}; 
    end

    % check dimension of set
    if dim(cPZ) ~= 2
       error('Method "conPolyZono/polygon" requires a 2D set!'); 
    end
    
    % remove analytical constraints
    L = analyticalConstraintElimination(cPZ);
    
    if length(L) > 1
       pgon = [];
       for i = 1:length(L)
          poly = polygon(L{i},splits);
          pgon = unite(pgon,poly);
       end
       return;
    else
       cPZ = L{1}; 
    end
    
    % remove redundant constraints
    cPZ = reduceConstraints(cPZ);
    
    % transform to equivalent higher-dimensional polynomial zonotope
    c = [cPZ.c; -cPZ.b];
    G = blkdiag(cPZ.G,cPZ.A);
    expMat = [cPZ.expMat,cPZ.expMat_];
    
    pZ = polyZonotope(c,G,[],expMat);
    
    % split the constrained polynomial zonotope multiple times to obtain a 
    % better over-approximation of the real shape
    warOrig = warning;
    warning('off','all');
    
    p = length(pZ.id); temp = ones(p,1);
    dom = interval(-temp,temp);
    
    [polyPrev,V] = getPolygon(pZ,cPZ.Grest);
    list{1}.set = pZ;
    list{1}.dom = dom;
    list{1}.V = V;

    for i = 1:splits
        
        list_ = [];
        polyAll = [];
        
        for j = 1:length(list)
                
            % only split the set if it is at the boundary
            if any(sum(ismembertol(list{j}.V',polyPrev.set.Vertices,1e-12),2) == 2)
            
                % split the polynomial zonotope
                res = splitDomain(list{j});

                % compute the corresponding polygons
                if intersectsNullSpace(res{1})
                    [poly,V] = getPolygon(res{1}.set,cPZ.Grest);
                    list_{end+1} = res{1};
                    list_{end}.V = V;
                    polyAll = unite(polyAll,poly);
                end

                if intersectsNullSpace(res{2})
                    [poly,V] = getPolygon(res{2}.set,cPZ.Grest);
                    list_{end+1} = res{2};
                    list_{end}.V = V;
                    polyAll = unite(polyAll,poly);
                end
            
            else
                list_{end+1} = list{j};
                polyAll = unite(polyAll,polygon(list{j}.V(1,:),list{j}.V(2,:)));
            end
        end
        
        % update list
        list = list_;
        polyPrev = polyAll;
    end
    
    warning(warOrig);
    
    % assign output arguments
    if ~isempty(polyAll)
        pgon = polyAll;
    else
        error('Set is empty!');
    end
end


% Auxiliary Functions -----------------------------------------------------

function [poly,V] = getPolygon(set,Grest)
% construct the polygon that corresponds to the splitted set

    % convert to zonotope
    zono = zonotope(set);
    zono = project(zono,[1,2]);
    zono = zonotope([zono.Z,Grest]);
    
    % convert zonotope to polygon
    V = vertices(zono);
    poly = polygon(V(1,:),V(2,:));
    
    % catch the case when the zonotope is degenerate
    if isempty(poly.set.Vertices)
        zono = zono + zonotope([0;0],eye(2)*1e-5);
        V = vertices(zono);
        poly = polygon(V(1,:),V(2,:));
    end
end

function res = unite(set1,set2)
% compute union of two polygons (also works for empty sets contrary to the
% original polygon/union function)

    if isempty(set1)
        res = set2;
    elseif isempty(set2)
        res = set1;
    else
        res = set1 | set2;
    end
end

function res = splitDomain(obj)
% split the factor domain of the constrained polynomial zonotope

   % split constraint polynomial zonotope along the longest generator
   [pZsplit,factor] = splitLongestGen(obj.set);
   
   % split the domain
   temp = split(obj.dom,factor);
   
   % construct output varible
   res = cell(2,1);
   
   res{1}.set = pZsplit{1}; res{1}.dom = temp{2};
   res{2}.set = pZsplit{2}; res{2}.dom = temp{1};
end

function res = intersectsNullSpace(obj)
% test if the split set violates the constraints (if it not intersects any
% of the hyperplanes)

    res = 1;
    n = length(obj.set.c);

    % loop over all constraint dimensions
    for i = 3:n
       
        c = zeros(n,1); c(i) = 1;
        hs = conHyperplane(c,0);
        
        if ~isIntersecting(hs,zonotope(obj.set))
           res = 0;
           return;
        end
    end
end

%------------- END OF CODE --------------