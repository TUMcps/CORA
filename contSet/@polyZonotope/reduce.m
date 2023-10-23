function pZ = reduce(pZ,option,order,varargin)
% reduce - reduces the order of a polynomial zonotope
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
%    figure; hold on;
%    plot(pZred,[1,2],'FaceColor','b');
%    plot(pZ,[1,2],'FaceColor','r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/reduce

% Authors:       Niklas Kochdumper
% Written:       23-March-2018 
% Last update:   06-July-2021 (MW, add adaptive)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % adaptive order reduction
    if strcmp(option,'adaptive')
        % note: var 'order' is not an order here
        pZ = reduceAdaptive(pZ,order);
        return;
    end

    if contains(option,'approxdep_')
        % remove independent generators
        pZ.GI = zeros(dim(pZ),0);
    end
    % extract dimensions
    N = length(pZ.c);
    P = size(pZ.G,2);
    Q = size(pZ.GI,2);

    % number of generators that stay unreduced (N generators are added again
    % after reduction)
    K = max(0,floor(N*order - N));

    % check if it is necessary to reduce the order
    if P + Q > N*order && K >= 0
        
        % concatenate all generators
        G = [pZ.G,pZ.GI];

        % half the generator length for exponents that are all even
        ind = ~any(mod(pZ.E,2),1);
        G(:,ind) = 0.5 * G(:,ind);

        % calculate the length of the generator vectors with a special metric
        len = sum(G.^2,1);

        % determine the smallest generators (= generators that are removed)
        [~,ind] = sort(len,'descend');
        ind = ind(K+1:end);

        % split the indices into the ones for dependent and independent
        % generators
        indDep = ind(ind <= P);
        indInd = ind(ind > P);
        indInd = indInd - P * ones(size(indInd));

        % construct a zonotope from the generators that are removed
        Grem = pZ.G(:,indDep);
        Erem = pZ.E(:,indDep);

        GIrem = pZ.GI(:,indInd);

        pZtemp = polyZonotope(zeros(N,1),Grem,GIrem,Erem);

        zono = zonotope(pZtemp);    % zonotope over-approximation

        % reduce the constructed zonotope with the reduction techniques for
        % linear zonotopes
        zonoRed = reduce(zono,option,1,varargin{:});

        % remove the generators that got reduced from the generator matrices
        pZ.G(:,indDep) = [];
        pZ.E(:,indDep) = [];
        pZ.GI(:,indInd) = [];

        % add the reduced generators as new independent generators 
        pZ.c = pZ.c + center(zonoRed);
        pZ.GI = [pZ.GI, generators(zonoRed)];

    end
    
    if contains(option,'approxdep_')
        % again remove rest generators
        pZ.GI = zeros(dim(pZ),0);
    end
    
    % remove all exponent vector dimensions that have no entries
    ind = sum(pZ.E,2)>0;
    pZ.E = pZ.E(ind,:);
    pZ.id = pZ.id(ind);
    
end

% ------------------------------ END OF CODE ------------------------------
