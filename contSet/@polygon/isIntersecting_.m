function res = isIntersecting_(pgon, S, type, tol, varargin)
% isIntersecting_ - check if two polygon objects intersect
%
% Syntax:
%    res = isIntersecting_(pgon, S)
%
% Inputs:
%    pgon - polygon
%    S - contSet
%
% Outputs:
%    res - logical
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/isIntersecting

% Authors:       Niklas Kochdumper, Tobias Ladner
% Written:       13-March-2020
% Last update:   11-October-2024 (contSet integration, generalized)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% ensure that numeric is second input argument
[pgon,S] = reorderNumeric(pgon,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < pgon.precedence
    res = isIntersecting_(S,pgon,type,tol);
    return
end

% numeric case: check containment
if isnumeric(S)
    res = contains_(pgon,S,type,tol,0,false,false);
    return
end

% sets must not be empty
if representsa_(pgon,'emptySet',0) || representsa_(S,'emptySet',0)
    res = false;
    return
end

% convert to polygon
S = polygon(S);

% expand boundaries
pgon = expandBoundaries(pgon,tol);

% compute intersection
pgon_intersection = polygon(intersect(pgon.set, S.set));
res = ~representsa_(pgon_intersection,"emptySet",tol);

end

% ------------------------------ END OF CODE ------------------------------
