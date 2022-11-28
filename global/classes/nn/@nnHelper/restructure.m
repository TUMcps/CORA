function [c, G, Grest, expMat, id] = restructure(c, G, Grest, ...
    expMat, id, id_, nrGen)
% restructure - Calculate a new over-approxmiating representation of a 
%    polynomial zonotope in such a way that there remain no independent 
%    generators
%
% Syntax:
%    [c, G, Grest, expMat, id] = nnHelper.restructure(c, G, Grest, expMat, id, id_, nrGen)
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
%    [c, G, Grest, expMat, id] - restructured polynomial zonotope
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:        Tobias Ladner
% Written:       31-May-2022
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

pZ = polyZonotope(c, G, Grest, expMat, id);
temp = restructure(pZ, 'reduceGirard', nrGen/length(c));

c = temp.c;
G = temp.G;
Grest = temp.Grest;
expMat = temp.expMat;
id = temp.id;
end

%------------- END OF CODE --------------