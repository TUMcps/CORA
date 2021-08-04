function res = isIntersectingMixed(E,S)
% isIntersectingMixed - checks whether E intersects contSet S
%
% Syntax:  
%    res = isIntersectingMixed(E,S)
%
% Inputs:
%    E - ellipsoid object
%    S - set
%
% Outputs:
%    res - boolean
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Victor Gassmann
% Written:      18-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
if E.isdegenerate
    % cannot do exactly, as e.g. intersection of zonotope with hyperplane
    % also not exactly possible
    % thus: simply bloat degenerate dimensions
    [T,D,~] = svd(E.Q);
    msv = max(diag(D));
    r = rank(E);
    n_rem = dim(E)-r;
    d = diag(D);
    D = diag([d(1:r);2*E.TOL*msv*ones(n_rem,1)]);
    E = ellipsoid(T*D*T',E.q);
end
% compute interval bounds for ellipsoid equation (x-q)'Q^-1(x-q) 
S = S + (-E.q);
temp = quadMap(S,{inv(E.Q)});
int = interval(temp);

% check if \exists x: 0 < (x-q)'Q^-1(x-q) < 1 -> intersection
if ~isempty(int & interval(0,1))
  res = true; 
else
  res = false; 
end   
%------------- END OF CODE --------------