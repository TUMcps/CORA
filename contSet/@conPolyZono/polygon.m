function pgon = polygon(cPZ,varargin)
% polygon - compute polygon enclosure of a two-dimensional constrained
%    polynomial zonotope
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
%    E = [1 0; 0 1; 0 0];
%    A = [1 1 -0.25];
%    b = 0.75;
%    EC = [2 0 0; 0 2 0; 0 0 1];
%    cPZ = conPolyZono(c,G,E,A,b,EC);
%
%    pgon = polygon(cPZ,15);
%    
%    plot(pgon,[1,2],'FaceColor','b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polygon, conPolyZono, zonotope

% Authors:       Niklas Kochdumper
% Written:       19-January-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % parse input arguments
    splits = setDefaultValues({10},varargin);

    % check input arguments
    inputArgsCheck({{cPZ,'att','conPolyZono'};
                    {splits,'att','numeric','nonnan'}});

    % check dimension of set
    if dim(cPZ) ~= 2
        throw(CORAerror('CORA:wrongValue','first','be a 2D set'));
    end
    
    % remove analytical constraints
    L = analyticalConstraintElimination(cPZ);
    
    if length(L) > 1
       pgon = [];
       for i = 1:length(L)
          poly = polygon(L{i},splits);
          pgon = aux_unite(pgon,poly);
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
    E = [cPZ.E,cPZ.EC];
    
    pZ = polyZonotope(c,G,[],E);
    
    % split the constrained polynomial zonotope multiple times to obtain a 
    % better over-approximation of the real shape
    warOrig = warning;
    warning('off','all');
    
    p = length(pZ.id); temp = ones(p,1);
    dom = interval(-temp,temp);
    
    [polyPrev,V] = aux_getPolygon(pZ,cPZ.GI);
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
                res = aux_splitDomain(list{j});

                % compute the corresponding polygons
                if aux_intersectsNullSpace(res{1})
                    [poly,V] = aux_getPolygon(res{1}.set,cPZ.GI);
                    list_{end+1} = res{1};
                    list_{end}.V = V;
                    polyAll = aux_unite(polyAll,poly);
                end

                if aux_intersectsNullSpace(res{2})
                    [poly,V] = aux_getPolygon(res{2}.set,cPZ.GI);
                    list_{end+1} = res{2};
                    list_{end}.V = V;
                    polyAll = aux_unite(polyAll,poly);
                end
            
            else
                list_{end+1} = list{j};
                polyAll = aux_unite(polyAll,polygon(list{j}.V(1,:),list{j}.V(2,:)));
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
        throw(CORAerror('CORA:emptySet'));
    end

end


% Auxiliary functions -----------------------------------------------------

function [poly,V] = aux_getPolygon(set,GI)
% construct the polygon that corresponds to the splitted set

    % convert to zonotope
    Z = zonotope(set);
    Z = project(Z,[1,2]);
    Z = zonotope(Z.c, [Z.G,GI]);
    
    % convert zonotope to polygon
    V = vertices(Z);
    poly = polygon(V(1,:),V(2,:));
    
    % catch the case when the zonotope is degenerate
    if isempty(poly.set.Vertices)
        Z = Z + zonotope([0;0],eye(2)*1e-5);
        V = vertices(Z);
        poly = polygon(V(1,:),V(2,:));
    end
end

function res = aux_unite(S1,S2)
% compute union of two polygons (also works for empty sets contrary to the
% original polygon/union function)

    if representsa_(S1,'emptySet',eps)
        res = S2;
    elseif representsa_(S2,'emptySet',eps)
        res = S1;
    else
        res = S1 | S2;
    end
end

function res = aux_splitDomain(obj)
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

function res = aux_intersectsNullSpace(obj)
% test if the split set violates the constraints (if it not intersects any
% of the hyperplanes)

    res = true;
    n = length(obj.set.c);

    % loop over all constraint dimensions
    for i = 3:n
       
        c = zeros(n,1); c(i) = 1;
        hs = conHyperplane(c,0);
        
        if ~isIntersecting_(hs,zonotope(obj.set),'exact')
           res = false;
           return;
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
