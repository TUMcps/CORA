function pZ = reduce(pZ,option,order,varargin)
% reduce - Reduces the order of a polynomial zonotope
%
% Syntax:  
%    pZ = reduce(pZ,option,order)
%
% Inputs:
%    pZ - polyZonotope object
%    option - reduction algorithm (see zonotope/reduce)
%    order - order of reduced polynomial zonotope
%
% Outputs:
%    pZ - reduced polynomial zonotope
%
% Example: 
%    pZ = polyZonotope([0;0],[2 0 1;0 1 1],[0.1,-0.4;0.2,0.3],[1 0 3;0 1 1]);
%    pZred = reduce(pZ,'girard',2);
%
%    hold on
%    plot(pZred,[1,2],'b','Filled',true,'EdgeColor','none');
%    plot(pZ,[1,2],'r','Filled',true,'EdgeColor','none');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/reduce

% Author:       Niklas Kochdumper
% Written:      23-March-2018 
% Last update:  06-July-2021 (MW, add adaptive)
% Last revision: ---

%------------- BEGIN CODE --------------

    % adaptive order reduction
    if strcmp(option,'adaptive')
        % note: var 'order' is not an order here
        pZ = reduceAdaptive(pZ,order);
        return;
    end

    global idv;
    if contains(option,'approxdep_')
        % remove independent generators
        pZ = polyZonotope(pZ.c,pZ.G,zeros(length(pZ.c),0),pZ.expMat,pZ.id);
        if ~isempty(idv)
            method = erase(option,'approxdep_');
            [pZ_v,pZ_nv] = onlyId(pZ,idv);
            o_v = ceil(order/2);
            o_nv = order-o_v;
            pZr_v = reduce(pZ_v,method,o_v);
            pZr_nv = reduce(pZ_nv,method,o_nv);
            pZ = exactPlus(pZr_v,pZr_nv);
            pZ = polyZonotope(pZ.c,pZ.G,zeros(length(pZ.c),1),pZ.expMat,pZ.id);
            return;
        end
    end
    % extract dimensions
    N = length(pZ.c);
    P = size(pZ.G,2);
    Q = size(pZ.Grest,2);

    % number of generators that stay unreduced (N generators are added again
    % after reduction)
    K = max(0,floor(N*order - N));

    % check if it is necessary to reduce the order
    if P + Q > N*order && K >= 0
        
        % concatenate all generators
        G = [pZ.G,pZ.Grest];

        % half the generator length for exponents that are all even
        temp = prod(ones(size(pZ.expMat))-mod(pZ.expMat,2),1);
        ind = find(temp == 1);
        G(:,ind) = 0.5 * G(:,ind);

        % calculate the length of the generator vectors with a special metric
        len = sum(G.^2,1);

        % determine the smallest generators (= generators that are removed)
        [~,ind] = sort(len,'descend');
        ind = ind(K+1:end);

        % split the indizes into the ones for dependent and independent
        % generators
        indDep = ind(ind <= P);
        indInd = ind(ind > P);
        indInd = indInd - P * ones(size(indInd));

        % construct a zonotope from the generators that are removed
        Grem = pZ.G(:,indDep);
        Erem = pZ.expMat(:,indDep);

        GrestRem = pZ.Grest(:,indInd);

        pZtemp = polyZonotope(zeros(N,1),Grem,GrestRem,Erem);

        zono = zonotope(pZtemp);    % zonotope over-approximation

        % reduce the constructed zontope with the reduction techniques for
        % linear zonotopes
        zonoRed = reduce(zono,option,1,varargin{:});

        % remove the generators that got reduced from the generator matrices
        pZ.G(:,indDep) = [];
        pZ.expMat(:,indDep) = [];
        pZ.Grest(:,indInd) = [];

        % add the reduced generators as new independent generators 
        pZ.c = pZ.c + center(zonoRed);
        pZ.Grest = [pZ.Grest, generators(zonoRed)];

    end
    if contains(option,'approxdep_')
        pZ = polyZonotope(pZ.c,pZ.G,zeros(length(pZ.c),0),pZ.expMat,pZ.id);
    else
        % remove all exponent vector dimensions that have no entries
        temp = sum(pZ.expMat,2);
        ind = find(temp > 0);
        pZ.expMat = pZ.expMat(ind,:);
        pZ.id = pZ.id(ind);
    end
end

%------------- END OF CODE --------------