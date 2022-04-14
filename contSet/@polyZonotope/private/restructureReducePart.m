function res = restructureReducePart(pZ, order, method)
% restructureReducePar - Calculate a new representation of a polynomial
%                        zonotope through partial reduction of the 
%                        independent generators 
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
%    res - resuling polyZonotope object which over-approximates pZ
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/restructure, zonotope/reduce

% Author:       Niklas Kochdumper
% Written:      18-January-2020 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    if isempty(pZ.Grest)
       res = pZ; return; 
    end
    
    % initialization
    n = dim(pZ); p = length(pZ.id); 

    % reduce independent generators to order 1
    zono = zonotope(zeros(dim(pZ),1),pZ.Grest);
    zono_ = reduce(zono,method,1);
    
    % check if upper bound for dependent factors is violated
    if p/n + 1 < order
        
        % order not violated -> define all independent generators as new
        % dependent genertors
        G = [pZ.G generators(zono_)];
        expMat = blkdiag(pZ.expMat,eye(n));
        id = [pZ.id; (max(pZ.id)+1:max(pZ.id) + n)'];
        
        res = polyZonotope(c,G,[],expMat,id);
        
    else
        
        % half the generator length for exponents that are all even
        Gtemp = pZ.G;         
        temp = prod(ones(size(pZ.expMat))-mod(pZ.expMat,2),1);
        ind = find(temp == 1);
        Gtemp(:,ind) = 0.5 * Gtemp(:,ind);
        
        % compute cost (= length of generators) of removing each dependent
        % factor
        cost = zeros(p,1);
        
        for i = 1:p
           cost(i) = sum(sum(Gtemp(:,pZ.expMat(i,:) > 0).^2,1));
        end
        
        % compare with cost (= length of genertors) of independent
        % generators
        Grest = generators(zono_);
        costInd = sum(Grest.^2,1)';
        
        [~,ind] = sort([cost;costInd],'descend');
        
        pmax = floor(order*p);
        indKeep = ind(1:pmax); indRem = ind(pmax+1:end);
        
        % keep all dependent genertors
        if all(indKeep <= p)                
            res = pZ;
            
        % keep all dependent generators + redefined some independent gens.
        elseif all(indRem > p) 
            
            index = indKeep(indKeep > p) - p;
            
            G = [pZ.G Grest(:,index)];
            expMat = blkdiag(pZ.expMat,eye(length(index)));
            id = [pZ.id; (max(pZ.id)+1:max(pZ.id) + length(index))'];
        
            res = polyZonotope(pZ.c,G,Grest(:,indRem-p),expMat,id);
            
        % remove smallest dependent generators    
        else
            
            % find dependent generators that are removed
            index = indRem(indRem <= p);
            ind = [];
            
            for i = index
                ind = [ind, find(pZ.expMat(i,:) > 0)];
            end
            
            ind = unique(ind);
            ind_ = setdiff(1:size(pZ.G,2),ind);
            
            % enclose generators that are removed with a zonotope
            pZ_ = polyZonotope(pZ.c,pZ.G(:,ind),pZ.Grest,pZ.expMat(:,ind));
            zono = zonotope(pZ_);
            zono_ = reduce(zono,method,1);
            
            % construct resulting zonotope
            temp = pmax - length(index);
            id = [pZ.id(index); (max(pZ.id)+1:max(pZ.id) + temp)'];
            
            res = polyZonotope(center(zono_),pZ.G(:,ind_), ...
                               generators(zono_),pZ.expMat(:,ind_),id);
        end
    end
end

%------------- END OF CODE --------------