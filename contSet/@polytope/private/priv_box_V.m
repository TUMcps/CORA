function [A,b,empty,fullDim,bounded] = priv_box_V(V,n)
% priv_box_V - computes the halfspace representation of the box enclosure
%    given a vertex representation
%
% Syntax:
%    [A,b,empty,fullDim,bounded] = priv_box_V(V,n)
%
% Inputs:
%    V - vertex representation
%    n - dimension of polytope
%
% Outputs:
%    A - inequality constraint matrix
%    b - inequality constraint offset
%    empty - true/false whether result is the empty set
%    fullDim - true/false on degeneracy
%    bounded - true/false on boundedness
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       03-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check for emptiness
empty = isempty(V);
if isempty(V)
    A = []; b = [];
    fullDim = false; bounded = true;
    return;
end

% compute lower and upper bound
ub = max(V,[],2);
lb = min(V,[],2);

% construct constraint matrix and offset
A = [eye(n); -eye(n)];
b = [ub; -lb];

% boundedness and degeneracy
bounded = ~any(isinf(b));
fullDim = ~any(withinTol(lb,ub,1e-10));

% ------------------------------ END OF CODE ------------------------------
