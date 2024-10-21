function [c, G, GI, E, id] = restructurePolyZono(c, G, GI, E, id, id_, nrGen)
% restructurePolyZono - Calculate a new over-approxmiating representation 
%    of a polynomial zonotope in such a way that there remain no 
%    independent generators
%
% Syntax:
%    [c, G, GI, E, id] = nnHelper.restructurePolyZono(c, G, GI, E, id, id_, nrGen)
%
% Inputs:
%    c - center of polyZonotope
%    G - dep. generator of polyZonotope
%    GI - indep. generator of polyZonotope
%    E - exponential matrix of polyZonotope
%    id - ids
%    id_ - max id
%    nrGen - nrGen
%
% Outputs:
%    [c, G, GI, E, id] - restructured polynomial zonotope
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       31-May-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% restructure
pZ = polyZonotope(c, G, GI, E, id);
pZ = restructure(pZ, 'reduceGirard', nrGen/length(c));

% read properties
c = pZ.c;
G = pZ.G;
GI = pZ.GI;
E = pZ.E;
id = pZ.id;

end

% ------------------------------ END OF CODE ------------------------------
