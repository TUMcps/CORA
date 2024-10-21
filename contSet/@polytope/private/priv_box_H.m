function [A,b,empty,fullDim,bounded] = priv_box_H(A,b,Ae,be,n)
% priv_box_H - computes the halfspace representation of the box enclosure
%    given a halfspace representation
%
% Syntax:
%    [A,b,empty] = priv_box_H(A,b,Ae,be,n)
%
% Inputs:
%    A - inequality constraint matrix
%    b - inequality constraint offset
%    Ae - equality constraint matrix
%    be - equality constraint offset
%    n - dimension of polytope
%
% Outputs:
%    A - inequality constraint matrix
%    b - inequality constraint offset
%    empty - true/false whether result is the empty set
%    fullDim - true/false on degeneracy
%    bounded - true/false on boundednesst
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

% init bounds
ub = Inf(n,1);
lb = -Inf(n,1);

% loop over all 2n positive/negative basis vectors
for i = 1:n
    % i-th basis vector
    e_i = unitvector(i,n);
    % maximize
    ub(i) = priv_supportFunc(A,b,Ae,be,e_i,'upper');
    if ub(i) == -Inf
        empty = true;
        A = []; b = [];
        fullDim = false; bounded = true;
        return
    end
    % minimize
    lb(i) = priv_supportFunc(A,b,Ae,be,e_i,'lower');
end

% construct output arguments
A = [eye(n); -eye(n)];
b = [ub; -lb];

% emptiness, boundedness, and degeneracy
empty = false;
bounded = ~any(isinf(b));
fullDim = ~any(withinTol(lb,ub,1e-10));

% ------------------------------ END OF CODE ------------------------------
