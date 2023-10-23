function res = isBounded(P)
% isBounded - determines whether a polytope is bounded
%
% Syntax:
%    res = isBounded(P)
%
% Inputs:
%    P - polytope object
%
% Outputs:
%    res - true/false whether the polytope is bounded
%
% Example: 
%    A = [2 1; 1 3; -1 2; -4 1];
%    b = ones(4,1);
%    P = polytope(A,b);
%    res = isBounded(P)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Viktor Kotsev, Mark Wetzlinger
% Written:       31-May-2022
% Last update:   30-November-2022 (MW, add 1D case)
%                08-December-2022 (MW, quick exit in nD case)
%                27-July-2023 (MW, fix 1D case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% check if property is set
if ~isempty(P.bounded.val)
    res = P.bounded.val;
end

% if polytope has V-rep then it is bounded (1D might still have Inf)
if ~isempty(P.V.val) && ~any(any(isinf(P.V.val)))
    P.bounded.val = true;
    return
end

% dimension of ambient space
n = dim(P);

% 1D case
if n == 1
    % compute vertices (and save)
    V = vertices_(P,'lcon2vert');
    P.V.val = V;
    % empty vertices or non-Inf vertices -> bounded
    if isempty(V)
        res = true;
        P.emptySet.val = true;
    elseif any(isinf(V))
        res = false;
        P.emptySet.val = false;
    else
        res = true;
        P.emptySet.val = false;
    end
    P.bounded.val = res;
    return
end

% loop over all dimensions
for i=1:n
    for e=[-1,1]
        % check every axis-aligned direction
        dir = [zeros(i-1,1);e;zeros(n-i,1)];
        % evaluate support function
        val = supportFunc_(P,dir,'upper');
        if val == Inf
            % set is unbounded
            res = false;
            P.bounded.val = false;
            P.emptySet.val = false;
            return
        elseif val == -Inf
            % set is empty
            res = true;
            P.bounded.val = true;
            P.emptySet.val = true;
            return
        end
    end
end

% code reaches this part: set is neither unbounded nor empty
res = true;
P.bounded.val = true;
P.emptySet.val = false;

% ------------------------------ END OF CODE ------------------------------
