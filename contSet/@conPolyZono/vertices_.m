function V = vertices_(cPZ,splits)
% vertices_ - computes the vertices of a two-dimensional constrained
%
% Syntax:
%    V = vertices_(cPZ)
%    V = vertices_(cPZ,splits)
%
% Inputs:
%    cPZ - conPolyZono object
%    splits - number of splits for refinement
%
% Outputs:
%    V - numeric, vertices
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/vertices

% Authors:       Niklas Kochdumper
% Written:       19-January-2020
% Last update:   13-March-2024 (TL, avoid timeout by reducing splitted set)
%                19-July-2024 (TL, keep collinear points in polygon)
%                11-October-2024 (TL, renamed to vertices_)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check dimension of set
if dim(cPZ) ~= 2
    throw(CORAerror('CORA:wrongValue','first','be a 2D set'));
end

% check if only point
[res,V] = representsa(cPZ,'point');
if res
    return;
end

% remove analytical constraints
L = priv_analyticalConstraintElimination(cPZ);

if length(L) > 1
   pgon = [];
   for i = 1:length(L)
      poly = polygon(L{i},splits);
      pgon = aux_unite(pgon,poly);
   end
   V = vertices_(pgon);
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
        if any(sum(ismembertol(list{j}.V,vertices(polyPrev),1e-12),1) == 2)
        
            % split the polynomial zonotope
            res = aux_splitDomain(list{j});

            % compute the corresponding polygons
            for k=1:numel(res)
                if aux_intersectsNullSpace(res{k})
                    [poly,V] = aux_getPolygon(res{k}.set,cPZ.GI);
                    list_{end+1} = res{k};
                    list_{end}.V = V;
                    polyAll = aux_unite(polyAll,poly);
                end
            end
        
        else
            list_{end+1} = list{j};
            polyAll = aux_unite(polyAll,polygon(list{j}.V,'KeepCollinearPoints',true));
        end
    end
    
    % update list
    list = list_;
    polyPrev = polyAll;
end

warning(warOrig);

% assign output arguments
if ~representsa_(polyAll,'emptySet',eps)
    pgon = polyAll;
else
    throw(CORAerror('CORA:emptySet'));
end

V = vertices_(pgon);

end


% Auxiliary functions -----------------------------------------------------

function [poly,V] = aux_getPolygon(set,GI)
% construct the polygon that corresponds to the splitted set

    % convert to zonotope
    Z = zonotope(set);
    Z = project(Z,[1,2]);
    Z = zonotope(Z.c, [Z.G,GI]);

    % reduce order for later union computations
    % 2*10000 generators should be enough for a 2d zonotope...
    Z = reduce(Z,'girard',10000); 

    % compact generators to remove any aligned/zero generators
    Z = compact_(Z,'all',1e-8);
    
    % convert zonotope to polygon (zonotope/vertices are already 'simple')
    V = vertices(Z);
    poly = polygon(V(1,:),V(2,:),'KeepCollinearPoints',true);
    
    % catch the case when the zonotope is degenerate
    if isempty(vertices_(poly))
        Z = Z + zonotope([0;0],eye(2)*1e-5);
        V = vertices(Z);
        poly = polygon(V(1,:),V(2,:),'KeepCollinearPoints',true);
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
        res = or(S1,S2,'KeepCollinearPoints',true);
    end
end

function res = aux_splitDomain(obj)
% split the factor domain of the constrained polynomial zonotope

   % split constraint polynomial zonotope along the longest generator
   [pZsplit,factor] = splitLongestGen(obj.set);
   
   % split the domain
   temp = split(obj.dom,factor);
   
   % construct output variable
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
        P = polytope([],[],unitvector(i,n)',0);
        if ~isIntersecting_(P,zonotope(obj.set),'exact',1e-8)
            res = false;
            return;
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
