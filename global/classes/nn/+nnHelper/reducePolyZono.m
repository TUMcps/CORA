function [c, G, GI, E, id, d] = reducePolyZono(c, G, GI, E, id, nrGen,S)
% reducePolyZono - reduce the number of generators of a polynomial zonotope, 
%    where we exploit that an interval remainder is added when reducing 
%    with Girards method
%
% Syntax:
%    [c, G, GI, E, id] = nnHelper.reducePolyZono(c, G, GI, E, id, id_, nrGen, S)
%
% Inputs:
%    c - center of polyZonotope
%    G - dep. generator of polyZonotope
%    GI - indep. generator of polyZonotope
%    E - exponential matrix of polyZonotope
%    id - ids
%    nrGen - nrGen
%    S - sensitivity
%
% Outputs:
%    [c, G, GI, E, id, d] - reduced polynomial zonotope
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Niklas Kochdumper, Tobias Ladner
% Written:       17-September-2021
% Last update:   01-August-2023 (TL, sensitivity aware order reduction)
% Last revision: 28-March-2022 (TL)
%                24-July-2023 (TL)

% ------------------------------ BEGIN CODE -------------------------------

if nargin < 7 || isempty(S)
    S = 1;
end

d = zeros(size(c));

if nrGen < size(G,2) + size(GI,2)
    % extract dimensions
    N = length(c);
    P = size(G,2);
    Q = size(GI,2);
    order = nrGen/N;

    % number of generators that stay unreduced (N generators are added again
    % after reduction)
    K = max(0,floor(N*order - N));

    % check if it is necessary to reduce the order
    if P + Q > N*order && K >= 0
        
        % concatenate all generators, weighted by sensitivity
        SG = S * [G,GI];

        % half the generator length for exponents that are all even
        indEven = ~any(mod(E,2),1);
        SG(:,indEven) = 0.5 * SG(:,indEven);

        % calculate the length of the generator vectors with a special metric
        len = sum(SG.^2,1);

        % determine the smallest generators (= generators that are removed)
        [~,indSmallest] = sort(len,'descend');
        indSmallest = indSmallest(K+1:end);

        % split the indices into the ones for dependent and independent
        % generators
        indDep = indSmallest(indSmallest <= P);
        indInd = indSmallest(indSmallest > P);
        indInd = indInd - P * ones(size(indInd));

        % construct a zonotope from the generators that are removed
        Grem = G(:,indDep);
        GIRem = GI(:,indInd);
        cRed = zeros(N,1);

        % half generators with all even exponents
        indEvenRem = indEven(indDep);
        Grem(:,indEvenRem) = 0.5 * Grem(:,indEvenRem);
        cRed = cRed + sum(0.5 * Grem(:,indEvenRem),2);

        % remove the generators that got reduced from the generator matrices
        G(:,indDep) = [];
        E(:,indDep) = [];
        GI(:,indInd) = [];

        % add shifted center
        c = c + cRed;

        % box over-approximation as approx error
        d = sum(abs([Grem,GIRem]),2);
    end
    
    % remove all exponent vector dimensions that have no entries
    ind = sum(E,2)>0;
    E = E(ind,:);
    id = id(ind);
end

end

% ------------------------------ END OF CODE ------------------------------
