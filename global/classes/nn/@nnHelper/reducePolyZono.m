function [c, G, Grest, expMat, id, d] = reducePolyZono(c, G, Grest, expMat, ...
    id, d, nrGen)
% reducePolyZono - reduce the number of generators of a polynomial zonotope, 
%    where we exploit that an interval remainder is added when reducing 
%    with Girards method
%
% Syntax:
%    [c, G, Grest, expMat, id] = nnHelper.reducePolyZono(c, G, Grest, expMat, id, id_, nrGen)
%
% Inputs:
%    c - center of polyZonotope
%    G - dep. generator of polyZonotope
%    Grest - indep. generator of polyZonotope
%    expMat - exponential matrix of polyZonotope
%    id - ids
%    d - error bound interval
%    nrGen - nrGen
%
% Outputs:
%    [c, G, Grest, expMat, id, d] - reduced polynomial zonotope
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:        Niklas Kochdumper, Tobias Ladner
% Written:       17-September-2021
% Last update:   ---
% Last revision: 28-March-2022 (TL)

%------------- BEGIN CODE --------------

n = length(c);

if ~isempty(nrGen) && nrGen < size(G, 2) + size(Grest, 2)
    temp = reduce(polyZonotope(c, G, Grest, expMat, id), 'girard', nrGen/n);
    c = temp.c;
    G = temp.G;
    expMat = temp.expMat;
    id = temp.id;
    m = max(1, size(temp.Grest, 2)-n);
    d = d + rad(interval(zonotope(zeros(n, 1), temp.Grest(:, m:end))));
    Grest = temp.Grest(:, 1:m-1);
end
end

%------------- END OF CODE --------------