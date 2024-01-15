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
A = [2;-1]; b = [4;1];
P1 = polytope(A,b);
A = [1;-2]; b = [5;1];
P2 = polytope(A,b);
% compute convex hull and compare to true result
P_conv = convHull(P1,P2);
A_true = [1;-1]; b_true = [5;1];
P_true = polytope(A_true,b_true);
res(end+1,1) = P_conv == P_true;

% 1D, only inequalities, bounded; only equalities, bounded
P1 = polytope([2;-1],[4;1]);
P2 = polytope([],[],1,8);
% compute convex hull and compare to true result
P_conv = convHull(P1,P2);
P_true = polytope([-1;1],[1;8]);
res(end+1,1) = P_conv == P_true;

% 1D, only inequalites, unbounded; only equalities, bounded
A = 2; b = 4;
P1 = polytope(A,b);
Ae = 1; be = 5;
P2 = polytope([],[],Ae,be);
% compute convex hull and compare to true result
P_conv = convHull(P1,P2);
A_true = 1; b_true = 5;
P_true = polytope(A_true,b_true);
res(end+1,1) = P_conv == P_true;


% 2D, fully empty (fullspace); bounded
A = zeros(0,2); b = zeros(0,0);
P1 = polytope(A,b);
A = [1 0; -1 1; -1 -1]; b = [1;1;1];
P2 = polytope(A,b);
P_conv = convHull(P1,P2);
res(end+1,1) = representsa(P_conv,'fullspace');
P_conv = convHull(P2,P1);
res(end+1,1) = representsa(P_conv,'fullspace');

% 2D, bounded, bounded
A = [1 0; 0 1;-1 0;0 -1]; b = [1;1;0;1];
P1 = polytope(A,b);
A = [1 0; 0 1;-1 0;0 -1]; b = [6;1;-5;1];
P2 = polytope(A,b);
% compute convex hull and compare to true result
P_conv = convHull(P1,P2);
A_true = [1 0;0 1;-1 0;0 -1]; b_true = [6;1;0;1];
P_true = polytope(A_true,b_true);
res(end+1,1) = P_conv == P_true;

% 2D, bounded, bounded degenerate
A = [1 1; -1 1; -1 -1; 1 -1]; b = ones(4,1);
P1 = polytope(A,b);
A = [0 1; 0 -1]; b = [2;2]; Ae = [1 0]; be = 5;
P2 = polytope(A,b,Ae,be);
% compute convex hull and compare to true result
P_conv = convHull(P1,P2);
V_ = [-1 0; 0 1; 5 2; 5 -2; 0 -1]';
P_true = polytope(V_);
res(end+1,1) = contains(P_conv,P_true,'exact',1e-8) ...
    && contains(P_true,P_conv,'exact',1e-8);

% 2D, bounded, empty
A = [1 1; -1 1; 0 -1]; b = [1;1;1];
P1 = polytope(A,b);
A = [0 1; -1 -1; 1 -1]; b = [-1;0.1;0.1];
P2 = polytope(A,b);
% compute convex hull and check emptiness
P_conv = convHull(P1,P2);
res(end+1,1) = representsa(P_conv,'emptySet');


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
