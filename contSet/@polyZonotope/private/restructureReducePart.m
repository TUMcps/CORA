function pZ = restructureReducePart(pZ, order, method)
% restructureReducePart - Calculate a new representation of a polynomial
%    zonotope through partial reduction of the independent generators 
%
% Syntax:
%    res = restructureReducePart(pZ, order, method)
%
% Inputs:
%    pZ - polyZonotope object
%    order - desired zonotope order of the dependent factors for the
%            resulting polynomial zonotope 
%    method - reduction technique for linear zonotopes
%             (see zonotope/reduce)
%
% Outputs:
%    pZ - polyZonotope object over-approximating input polynomial zonotope
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/restructure, zonotope/reduce

% Authors:       Niklas Kochdumper
% Written:       18-January-2020 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if isempty(pZ.GI)
    return; 
end

% initialization
n = dim(pZ); p = length(pZ.id); 

% reduce independent generators to order 1
zono = zonotope(zeros(dim(pZ),1),pZ.GI);
zono_ = reduce(zono,method,1);

% check if upper bound for dependent factors is violated
if p/n + 1 < order
    
    % order not violated -> define all independent generators as new
    % dependent genertors
    G = [pZ.G generators(zono_)];
    E = blkdiag(pZ.E,eye(n));
    id = [pZ.id; (max(pZ.id)+1:max(pZ.id) + n)'];
    
    pZ = polyZonotope(c,G,[],E,id);
    
else
    
    % half the generator length for exponents that are all even
    Gtemp = pZ.G;         
    temp = prod(ones(size(pZ.E))-mod(pZ.E,2),1);
    ind = find(temp == 1);
    Gtemp(:,ind) = 0.5 * Gtemp(:,ind);
    
    % compute cost (= length of generators) of removing each dependent
    % factor
    cost = zeros(p,1);
    
    for i = 1:p
       cost(i) = sum(sum(Gtemp(:,pZ.E(i,:) > 0).^2,1));
    end
    
    % compare with cost (= length of genertors) of independent
    % generators
    GI = zono_.G;
    costInd = sum(GI.^2,1)';
    
    [~,ind] = sort([cost;costInd],'descend');
    
    pmax = floor(order*p);
    indKeep = ind(1:pmax); indRem = ind(pmax+1:end);
    
    % keep all dependent genertors
    if all(indKeep <= p)                
        pZ = pZ;
        
    % keep all dependent generators + redefined some independent gens.
    elseif all(indRem > p) 
        
        index = indKeep(indKeep > p) - p;
        
        G = [pZ.G GI(:,index)];
        E = blkdiag(pZ.E,eye(length(index)));
        id = [pZ.id; (max(pZ.id)+1:max(pZ.id) + length(index))'];
    
        pZ = polyZonotope(pZ.c,G,GI(:,indRem-p),E,id);
        
    % remove smallest dependent generators    
    else
        
        % find dependent generators that are removed
        index = indRem(indRem <= p);
        ind = [];
        
        for i = index
            ind = [ind, find(pZ.E(i,:) > 0)];
        end
        
        ind = unique(ind);
        ind_ = setdiff(1:size(pZ.G,2),ind);
        
        % enclose generators that are removed with a zonotope
        pZ_ = polyZonotope(pZ.c,pZ.G(:,ind),pZ.GI,pZ.E(:,ind));
        zono = zonotope(pZ_);
        zono_ = reduce(zono,method,1);
        
        % construct resulting zonotope
        temp = pmax - length(index);
        id = [pZ.id(index); (max(pZ.id)+1:max(pZ.id) + temp)'];
        
        pZ = polyZonotope(zono_.c, pZ.G(:,ind_), zono_.G, pZ.E(:,ind_),id);
    end
end

% ------------------------------ END OF CODE ------------------------------
