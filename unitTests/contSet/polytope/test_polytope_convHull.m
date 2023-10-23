function res = test_polytope_convHull
% test_polytope_convHull - unit test function of convHull
%
% Syntax:
%    res = test_polytope_convHull
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
% See also: none

% Authors:       Viktor Kotsev, Mark Wetzlinger
% Written:       30-October-2022
% Last update:   27-July-2023 (MW, more tests)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% 1D, only inequalities, bounded; only inequalities, bounded
P1 = polytope([2;-1],[4;1]);
P2 = polytope([1;-2],[5;1]);
% compute convex hull and compare to true result
P = convHull(P1,P2);
P_ = polytope([1;-1],[5;1]);
res(end+1,1) = P == P_;

% 1D, only inequalities, bounded; only equalities, bounded
P1 = polytope([2;-1],[4;1]);
P2 = polytope([],[],1,8);
% compute convex hull and compare to true result
P = convHull(P1,P2);
P_ = polytope([-1;1],[1;8]);
res(end+1,1) = P == P_;

% 1D, only inequalites, unbounded; only equalities, bounded
P1 = polytope(2,4);
P2 = polytope([],[],1,5);
% compute convex hull and compare to true result
P = convHull(P1,P2);
P_ = polytope(1,5);
res(end+1,1) = P == P_;


% 2D, bounded, bounded
P1 = polytope([1 0; 0 1;-1 0;0 -1],[1;1;0;1]);
P2 = polytope([1 0; 0 1;-1 0;0 -1],[6;1;-5;1]);
% compute convex hull and compare to true result
P = convHull(P1,P2);
P_ = polytope([1 0;0 1;-1 0;0 -1],[6;1;0;1]);
res(end+1,1) = P == P_;

% 2D, bounded, bounded degenerate
P1 = polytope([1 1; -1 1; -1 -1; 1 -1],ones(4,1));
P2 = polytope([0 1; 0 -1],[2;2],[1 0],5);
% compute convex hull and compare to true result
P = convHull(P1,P2);
V_ = [-1 0; 0 1; 5 2; 5 -2; 0 -1]';
P_ = polytope(V_);
res(end+1,1) = contains(P,P_,'exact',1e-8) && contains(P_,P,'exact',1e-8); % equivalent to "=="

% 2D, bounded, empty
P1 = polytope([1 1; -1 1; 0 -1],[1;1;1]);
P2 = polytope([0 1; -1 -1; 1 -1],[-1;0.1;0.1]);
% compute convex hull and check emptiness
P = convHull(P1,P2);
res(end+1,1) = representsa(P,'emptySet');


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
