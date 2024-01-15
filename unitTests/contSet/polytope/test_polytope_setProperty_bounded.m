function res = test_polytope_setProperty_bounded
% test_polytope_setProperty_bounded - unit test function to check whether
%    the internally-used set property 'bounded' is changed correctly
%    following different set operations on a polytope
%
% Syntax:
%    res = test_polytope_setProperty_bounded
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

% Authors:       Mark Wetzlinger
% Written:       01-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% --- polytope ------------------------------------------------------------
% init set via vertex representation

% 1D, unbounded, non-degenerate
V = [-Inf, 2];
P = polytope(V);
res(end+1,1) = ~isempty(P.bounded.val) && ~P.bounded.val;

% 1D, bounded, non-degenerate
V = [-3 -1 4 5];
P = polytope(V);
res(end+1,1) = ~isempty(P.bounded.val) && P.bounded.val;

% 1D, bounded, degenerate
V = 2;
P = polytope(V);
res(end+1,1) = ~isempty(P.bounded.val) && P.bounded.val;

% 2D, bounded, non-degenerate
V = [2 1; -1 4; -4 0; -1 -2; 3 -1]';
P = polytope(V);
res(end+1,1) = ~isempty(P.bounded.val) && P.bounded.val;

% 2D, bounded, degenerate
V = [-1 1; 2 0]';
P = polytope(V);
res(end+1,1) = ~isempty(P.bounded.val) && P.bounded.val;

% 3D, bounded (with redundancies)
V = [1 1 -1; -3 2 1; 0 4 2; 2 3 1; 2 -2 1; 4 3 2; -1 1 1]';
P = polytope(V);
res(end+1,1) = ~isempty(P.bounded.val) && P.bounded.val;


% --- and -----------------------------------------------------------------

% 2D, bounded, non-degenerate & bounded, non-degenerate
P1 = polytope([1 1; 1 -1; -1 0],[1;1;0.5]);
P2 = polytope([-1 1; -1 -1; 1 0],[1;1;0.5]);
% determine boundedness
isBounded(P1); isBounded(P2);
P = P1 & P2;
% intersection also known to be bounded
res(end+1,1) = ~isempty(P.bounded.val) && P.bounded.val;

% 2D, bounded, nondegenerate & empty
P1 = polytope([1 1; 1 -1; -1 0],[1;1;0.5]);
P2 = polytope([1 0; -1 0],[2;-3]);
% determine boundedness
isBounded(P1); isBounded(P2);
P = P1 & P2;
% intersection also known to be bounded
res(end+1,1) = ~isempty(P.bounded.val) && P.bounded.val;


% --- box -----------------------------------------------------------------

% 2D, empty set
P = polytope([1 0; -1 0],[2;-3]);
P_ = box(P);
res(end+1,1) = ~isempty(P.bounded.val) && P.bounded.val;
res(end+1,1) = ~isempty(P_.bounded.val) && P_.bounded.val;

% 2D, unbounded, non-degenerate, non-empty
P = polytope([1 0 0; 0 1 0],[1;-3]);
P_ = box(P);
res(end+1,1) = ~isempty(P.bounded.val) && ~P.bounded.val;
res(end+1,1) = ~isempty(P_.bounded.val) && ~P_.bounded.val;

% 2D, bounded, degenerate, non-empty
P = polytope([1 0 0; 0 1 0; -1 -1 0],[1;2;2],[0,0,1],0);
P_ = box(P);
res(end+1,1) = ~isempty(P.bounded.val) && P.bounded.val;
res(end+1,1) = ~isempty(P_.bounded.val) && P_.bounded.val;


% --- compact -------------------------------------------------------------

% 1D, empty
P = polytope([1;-1],[1;-2]);
P_ = compact(P);
res(end+1,1) = ~isempty(P.bounded.val) && P.bounded.val ...
    && ~isempty(P_.bounded.val) && P_.bounded.val;
% 3D, empty
P = polytope([1 1 0; -1 1 0; 0 -1 0; 1 0 0],[1;1;1;-3],[0 0 1],3);
P_ = compact(P);
res(end+1,1) = ~isempty(P.bounded.val) && P.bounded.val ...
    && ~isempty(P_.bounded.val) && P_.bounded.val;


% --- convHull ------------------------------------------------------------

% 2D, bounded, empty
P1 = polytope([1 1; -1 1; 0 -1],[1;1;1]);
P2 = polytope([0 1; -1 -1; 1 -1],[-1;0.1;0.1]);
P = convHull(P1,P2);
% P2 and P are bounded
res(end+1,1) = ~isempty(P2.bounded.val) && P2.bounded.val ...
    && ~isempty(P.bounded.val) && P.bounded.val;

% 1D, bounded, bounded
P1 = polytope([],[],1,2);
P2 = polytope([1;-1],[0;5]);
P = convHull(P1,P2);
% result is bounded
res(end+1,1) = ~isempty(P.bounded.val) && P.bounded.val;


% --- empty ---------------------------------------------------------------

n = 3;
P = polytope.empty(n);
res(end+1,1) = ~isempty(P.bounded.val) && P.bounded.val;


% --- Inf -----------------------------------------------------------------

n = 3;
P = polytope.Inf(n);
res(end+1,1) = ~isempty(P.bounded.val) && ~P.bounded.val;


% --- isBounded -----------------------------------------------------------

% 2D, empty
P = polytope([0 1; -1 -1; 1 -1],[-1;0.1;0.1]);
% determine boundedness
isBounded(P);
res(end+1,1) = ~isempty(P.bounded.val) && P.bounded.val;

% 3D, unbounded
P = polytope([],[],[1 1 1; -1 1 1],[1;1]);
% determine boundedness
isBounded(P);
res(end+1,1) = ~isempty(P.bounded.val) && ~P.bounded.val;


% --- isFullDim -----------------------------------------------------------

% 1D, unbounded
P = polytope(1,1);
% determine boundedness
isFullDim(P);
res(end+1,1) = ~isempty(P.bounded.val) && ~P.bounded.val;


% --- plus ----------------------------------------------------------------

% 2D, bounded + vector
A = [1 0; -1 1; -1 -1]; b = ones(3,1);
P = polytope(A,b);
isBounded(P);
v = [-1;1];
P_sum = P + v;
% resulting polytope is also bounded
res(end+1,1) = ~isempty(P_sum.bounded.val) && P_sum.bounded.val;


% --- polytope ------------------------------------------------------------
% copy constructor

% 2D, only inequalities, bounded
P = polytope([1 1; -1 1; 0 -1],ones(3,1));
% determine boundedness
isBounded(P);
% copy polytope, property should also be copied
P_ = polytope(P);
res(end+1,1) = ~isempty(P_.bounded.val) && P_.bounded.val;


% --- project -------------------------------------------------------------

% 3D, empty
P = polytope([1 0 0; -1 0 0],[2;-3]);
P_ = project(P,[1,2]);
res(end+1,1) = ~isempty(P_.bounded.val) && P_.bounded.val;


% --- lift ----------------------------------------------------------------

% 2D, non-empty
P = polytope([1 1; -1 1; 0 -1],[1;1;1]);
% project to higher-dimensional space
P_ = lift(P,5,[2,3]);
% higher-dimensional polytope always unbounded
res(end+1,1) = ~isempty(P_.bounded.val) && ~P_.bounded.val;


% --- representsa ---------------------------------------------------------

% 2D, empty
P = polytope([1 1; -1 1; 0 -1],[1;1;-2]);
% determine boundedness via emptiness
representsa(P,'emptySet');
res(end+1,1) = ~isempty(P.bounded.val) && P.bounded.val;

% 3D, unconstrained
P = polytope([0,0,0],4);
% polytope is R^3, thus unbounded
representsa(P,'fullspace');
res(end+1,1) = ~isempty(P.bounded.val) && ~P.bounded.val;


% --- vertices ------------------------------------------------------------

% 2D, empty set
P = polytope([1 0; -1 0],[2;-3]);
V = vertices(P);
res(end+1,1) = ~isempty(P.bounded.val) && P.bounded.val;

% 2D, bounded, non-degenerate
P = polytope([1 0; 0 1; -1 -1],[2;1;2]);
V = vertices(P);
res(end+1,1) = ~isempty(P.bounded.val) && P.bounded.val;

% 2D, bounded, degenerate
P = polytope([1 0; -1 0],[2;1],[0 1],3);
V = vertices(P);
res(end+1,1) = ~isempty(P.bounded.val) && P.bounded.val;


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
