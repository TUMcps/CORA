function res = test_polytope_reduce
% test_polytope_reduce - unit test function of reduce
%
% Syntax:
%    res = test_polytope_reduce
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: test_polytope_hausdorffDist

% Authors:       Tobias Ladner
% Written:       08-December-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tol = 1e-8;

% 1D, empty
A = [1; -1]; b = [2; -5];
P = polytope(A,b);
P_reduce = reduce(P,'rand',2);
assert(representsa(P_reduce,'emptySet'))

% 1D, bounded
A = [1; -1]; b = [1; 1];
P = polytope(A,b);
P_reduce = reduce(P,'rand',2);
assert(contains(P_reduce,P,'exact',tol))

% 1D, unbounded
A = 1; b = 1;
P = polytope(A,b);
P_reduce = reduce(P,'rand',2);
assert(contains(P_reduce,P,'exact',tol))

% 2D, empty
A = [1 0; -1 1; -1 -1; 0 1]; b = [1;1;1;-4];
P = polytope(A,b);
P_reduce = reduce(P,'rand',2);
assert(contains(P_reduce,P,'exact',tol))

% 2D, bounded
A = [1 0;-1 0; 0 1; 0 -1]; b = [6;-5;1;1];
P = polytope(A,b);
P_reduce = reduce(P,'rand',2);
assert(contains(P_reduce,P,'exact',tol))

% 2D, bounded & degenerate
A = [1 0; 0 1; -1 0; 0 -1]; b = [1;1;1;-1];
P = polytope(A,b);
P_reduce = reduce(P,'rand',2);
assert(contains(P_reduce,P,'exact',tol))

% 2D, unbounded
A = [1 0;0 1;-1 0]; b = [1;1;1];
P = polytope(A,b);
P_reduce = reduce(P,'rand',2);
assert(contains(P_reduce,P,'exact',tol))

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
