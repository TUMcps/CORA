function res = priv_isIntersectingMixed(E,S)
% priv_isIntersectingMixed - checks whether an ellipsoid intersects another
%    set
%
% Syntax:
%    res = priv_isIntersectingMixed(E,S)
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
    r = rank(E);
    n_rem = dim(E) - r;
    d = diag(D);
    D = diag([d(1:r); 2*E.TOL*max(d)*ones(n_rem,1)]);
    E = ellipsoid(T*D*T',E.q);
end
% compute interval bounds for ellipsoid equation (x-q)'Q^-1(x-q) 
S = S + (-E.q);
I = interval(quadMap(S,{inv(E.Q)}));

% check if \exists x: 0 < (x-q)'Q^-1(x-q) < 1 -> intersection
res = ~representsa_(and_(I,interval(0,1),'exact'),'emptySet',eps);

% ------------------------------ END OF CODE ------------------------------
