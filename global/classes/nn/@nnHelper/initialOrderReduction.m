function [c, G, Grest, expMat, id] = initialOrderReduction(c, G, Grest, ...
    expMat, id, id_, nrGen)
% initialOrderReduction - order reduction in the first layer, where we 
%    potentially redefine large independent generators as new dependent 
%    generators
%
% Syntax:
%    [c, G, Grest, expMat, id] = nnHelper.initialOrderReduction(c, G, Grest, expMat, id, id_, nrGen)
%
% Inputs:
%    c - center of polyZonotope
%    G - dep. generator of polyZonotope
%    Grest - indep. generator of polyZonotope
%    expMat - exponential matrix of polyZonotope
%    id - ids
%    id_ - max id
%    nrGen - nrGen
%
% Outputs:
%    [c, G, Grest, expMat, id] - reduced polynomial zonotope
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

temp = reduce(polyZonotope(c, G, Grest, expMat, id), ...
    'girard', nrGen/length(c));

if isempty(temp.G)
    c = temp.c;
    G = temp.Grest;
    Grest = zeros(size(c));
    expMat = eye(size(G, 2));
    id = id_ + (1:size(G, 2))';
else
    c = temp.c;
    G = temp.G;
    Grest = temp.Grest;
    expMat = temp.expMat;
    id = temp.id;
end
end

%------------- END OF CODE --------------