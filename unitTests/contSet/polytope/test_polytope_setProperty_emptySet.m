function res = test_polytope_setProperty_emptySet
% test_polytope_setProperty_emptySet - unit test function to check whether
%    the internally-used set property 'emptySet' is changed correctly
%    following different set operations on a polytope
%
% Syntax:
%    res = test_polytope_setProperty_emptySet
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
% Written:       31-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% --- polytope ------------------------------------------------------------
% init set via vertex representation -> always non-empty

% 1D, unbounded, non-degenerate
V = [-Inf, 2];
P = polytope(V);
res(end+1,1) = ~isempty(P.emptySet.val) && ~P.emptySet.val;

% 1D, bounded, non-degenerate
V = [-3 -1 4 5];
P = polytope(V);
res(end+1,1) = ~isempty(P.emptySet.val) && ~P.emptySet.val;

% 1D, bounded, degenerate
V = 2;
P = polytope(V);
res(end+1,1) = ~isempty(P.emptySet.val) && ~P.emptySet.val;

% 2D, bounded, non-degenerate
V = [2 1; -1 4; -4 0; -1 -2; 3 -1]';
P = polytope(V);
res(end+1,1) = ~isempty(P.emptySet.val) && ~P.emptySet.val;

% 2D, bounded, degenerate
V = [-1 1; 2 0]';
P = polytope(V);
res(end+1,1) = ~isempty(P.emptySet.val) && ~P.emptySet.val;


% --- box -----------------------------------------------------------------

% 2D, empty set
P = polytope([1 0; -1 0],[2;-3]);
P_ = box(P);
res(end+1,1) = ~isempty(P.emptySet.val) && P.emptySet.val;
res(end+1,1) = ~isempty(P_.emptySet.val) && P_.emptySet.val;

% 2D, unbounded, non-degenerate, non-empty
P = polytope([1 0 0; 0 1 0],[1;-3]);
P_ = box(P);
res(end+1,1) = ~isempty(P.emptySet.val) && ~P.emptySet.val;
res(end+1,1) = ~isempty(P_.emptySet.val) && ~P_.emptySet.val;

% 2D, bounded, degenerate, non-empty
P = polytope([1 0 0; 0 1 0; -1 -1 0],[1;2;2],[0,0,1],0);
P_ = box(P);
res(end+1,1) = ~isempty(P.emptySet.val) && ~P.emptySet.val;
res(end+1,1) = ~isempty(P_.emptySet.val) && ~P_.emptySet.val;


% --- vertices ------------------------------------------------------------

% 2D, empty set
P = polytope([1 0; -1 0],[2;-3]);
V = vertices(P);
res(end+1,1) = ~isempty(P.emptySet.val) && P.emptySet.val;

% 2D, bounded, non-degenerate
P = polytope([1 0; 0 1; -1 -1],[2;1;2]);
V = vertices(P);
res(end+1,1) = ~isempty(P.emptySet.val) && ~P.emptySet.val;

% 2D, bounded, degenerate
P = polytope([1 0; -1 0],[2;1],[0 1],3);
V = vertices(P);
res(end+1,1) = ~isempty(P.emptySet.val) && ~P.emptySet.val;


% --- and -----------------------------------------------------------------

% 2D, bounded, non-degenerate & bounded, non-degenerate
P1 = polytope([1 1; 1 -1; -1 0],[1;1;0.5]);
P2 = polytope([-1 1; -1 -1; 1 0],[1;1;0.5]);
P = P1 & P2;
% emptiness unknown
res(end+1,1) = isempty(P.emptySet.val);

% 2D, bounded, nondegenerate & empty
P1 = polytope([1 1; 1 -1; -1 0],[1;1;0.5]);
P2 = polytope([1 0; -1 0],[2;-3]);
% get knowledge about emptiness
representsa(P2,'emptySet');
% compute intersection
P = P1 & P2;
% emptiness known via emptiness of P2
res(end+1,1) = ~isempty(P.emptySet.val) && P.emptySet.val;


% --- compact -------------------------------------------------------------

% 1D, empty
P = polytope([1;-1],[1;-2]);
P_ = compact(P);
res(end+1,1) = ~isempty(P.emptySet.val) && P.emptySet.val ...
    && ~isempty(P_.emptySet.val) && P_.emptySet.val;
% 3D, empty
P = polytope([1 1 0; -1 1 0; 0 -1 0; 1 0 0],[1;1;1;-3],[0 0 1],3);
P_ = compact(P);
res(end+1,1) = ~isempty(P.emptySet.val) && P.emptySet.val ...
    && ~isempty(P_.emptySet.val) && P_.emptySet.val;


% --- convHull ------------------------------------------------------------

% 2D, bounded, empty
P1 = polytope([1 1; -1 1; 0 -1],[1;1;1]);
P2 = polytope([0 1; -1 -1; 1 -1],[-1;0.1;0.1]);
P = convHull(P1,P2);
% P2 and P are empty
res(end+1,1) = ~isempty(P2.emptySet.val) && P2.emptySet.val ...
    && ~isempty(P.emptySet.val) && P.emptySet.val;

% 1D, bounded, bounded
P1 = polytope([],[],1,2);
P2 = polytope([1;-1],[0;5]);
P = convHull(P1,P2);
% result is non-empty
res(end+1,1) = ~isempty(P.emptySet.val) && ~P.emptySet.val;


% --- empty ---------------------------------------------------------------

n = 3;
P = polytope.empty(n);
res(end+1,1) = ~isempty(P.emptySet.val) && P.emptySet.val;


% --- Inf -----------------------------------------------------------------

n = 3;
P = polytope.Inf(n);
res(end+1,1) = ~isempty(P.emptySet.val) && ~P.emptySet.val;


% --- isBounded -----------------------------------------------------------

% 2D, empty
P = polytope([0 1; -1 -1; 1 -1],[-1;0.1;0.1]);
% determine boundedness (and emptiness)
isBounded(P);
res(end+1,1) = ~isempty(P.emptySet.val) && P.emptySet.val;


% --- isFullDim -----------------------------------------------------------

% 2D, empty
P = polytope([0 1; -1 -1; 1 -1],[-1;0.1;0.1]);
% determine full-dimensionality (and emptiness)
isFullDim(P);
res(end+1,1) = ~isempty(P.emptySet.val) && P.emptySet.val;


% --- plus ----------------------------------------------------------------

% 2D, bounded + vector
A = [1 0; -1 1; -1 -1]; b = ones(3,1);
P = polytope(A,b);
representsa(P,'emptySet');
v = [-1;1];
P_sum = P + v;
% resulting polytope is also non-empty
res(end+1,1) = ~isempty(P_sum.emptySet.val) && ~P_sum.emptySet.val;


% --- polytope ------------------------------------------------------------
% copy constructor

% 2D, only inequalities, non-empty
P = polytope([1 1; -1 1; 0 -1],ones(3,1));
% determine emptiness
representsa(P,'emptySet');
% copy polytope, property should also be copied
P_ = polytope(P);
res(end+1,1) = ~isempty(P_.emptySet.val) && ~P_.emptySet.val;


% --- project -------------------------------------------------------------

% 3D, empty
P = polytope([1 0 0; -1 0 0],[2;-3]);
P_ = project(P,[1,2]);
res(end+1,1) = ~isempty(P_.emptySet.val) && P_.emptySet.val;


% --- lift ----------------------------------------------------------------

% 2D, non-empty
P = polytope([1 1; -1 1; 0 -1],[1;1;1]);
% get knowledge about emptiness
representsa(P,'emptySet');
% project to higher-dimensional space
P_ = lift(P,5,[2,3]);
% higher-dimensional polytope also non-empty
res(end+1,1) = ~isempty(P_.emptySet.val) && ~P_.emptySet.val;


% --- representsa ---------------------------------------------------------

P = polytope([1 1; -1 1; 0 -1],ones(3,1));
% determine emptiness
representsa(P,'emptySet');
res(end+1,1) = ~isempty(P.emptySet.val) && ~P.emptySet.val;


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
