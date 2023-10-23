function res = isIntersectingMixed(E,S)
% isIntersectingMixed - checks whether an ellipsoid intersects another set
%
% Syntax:
%    res = isIntersectingMixed(E,S)
%
% Inputs:
%    E - ellipsoid object
%    S - contSet object
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Victor Gassmann
% Written:       18-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if ~isFullDim(E)
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
I = interval(temp);

% check if \exists x: 0 < (x-q)'Q^-1(x-q) < 1 -> intersection
if ~representsa_(and_(I,interval(0,1),'exact'),'emptySet',eps)
    res = true; 
else
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
