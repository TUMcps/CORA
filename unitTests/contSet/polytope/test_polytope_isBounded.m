function res = test_polytope_isBounded
% test_polytope_isBounded - unit test function for determination of
%    boundedness
%
% Syntax:
%    res = test_polytope_isBounded
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
% Written:       29-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% 1D, unbounded
A = 1; b = 1;
P = polytope(A,b);
res(end+1,1) = ~isBounded(P) && ~isempty(P.bounded.val) && ~P.bounded.val;

% 1D, bounded
A = [1;-1]; b = [1;1];
P = polytope(A,b);
res(end+1,1) = isBounded(P) && ~isempty(P.bounded.val) && P.bounded.val;

% 1D, single point
Ae = 5; be = 2;
P = polytope([],[],Ae,be);
res(end+1,1) = isBounded(P) && ~isempty(P.bounded.val) && P.bounded.val;

% 1D, fully empty
A = zeros(0,1); b = zeros(0,0);
P = polytope(A,b);
res(end+1,1) = ~isBounded(P) && ~isempty(P.bounded.val) && ~P.bounded.val;


% 2D, unbounded, only equality constraints (axis-aligned)
Ae = [1 0]; be = 2;
P = polytope([],[],Ae,be);
res(end+1,1) = ~isBounded(P) && ~isempty(P.bounded.val) && ~P.bounded.val;

% 2D, bounded, only equality constraints (axis-aligned, single point)
Ae = [1 0; 0 1]; be = [2; -1];
P = polytope([],[],Ae,be);
res(end+1,1) = isBounded(P) && ~isempty(P.bounded.val) && P.bounded.val;

% 2D, unbounded, only equality constraints
Ae = [1 1]; be = 2;
P = polytope([],[],Ae,be);
res(end+1,1) = ~isBounded(P) && ~isempty(P.bounded.val) && ~P.bounded.val;

% 2D, unbounded (towards -x2)
A = [2 1; 1 3; -1 2; -4 1]; b = ones(4,1);
P = polytope(A,b);
res(end+1,1) = ~isBounded(P) && ~isempty(P.bounded.val) && ~P.bounded.val;

% 2D, bounded
A = [1 1; -2 1; -4 -2; 2 -3]; b = ones(4,1);
P = polytope(A,b) + [10;5];
res(end+1,1) = isBounded(P) && ~isempty(P.bounded.val) && P.bounded.val;

% 2D, unbounded, degenerate
A = [1 0; 0 1; 0 -1]; b = [1; 1; -1];
P = polytope(A,b);
res(end+1,1) = ~isBounded(P) && ~isempty(P.bounded.val) && ~P.bounded.val;

% 2D, bounded, degenerate using equality constraint
A = [0 1; 0 -1]; b = [3;-1]; Ae = [1 0]; be = -2;
P = polytope(A,b,Ae,be);
res(end+1,1) = isBounded(P) && ~isempty(P.bounded.val) && P.bounded.val;


% 3D, rotated unit cube
n = 3; [M,~,~] = svd(randn(n));
P = M * polytope([eye(n); -eye(n)],ones(2*n,1));
res(end+1,1) = isBounded(P) && ~isempty(P.bounded.val) && P.bounded.val;


% Sequence of functions
% 2D, bounded
A = [1 1; -2 1; -4 -2; 2 -3]; b = ones(4,1);
P = polytope(A,b);
isBounded(P);

A = [1 0; -1 0; 0 1; 0 -1]; b = [1;-1;1;-1];
P2 = polytope(A,b);
isBounded(P2);

% intersection should still be bounded (single point)
P = P & P2;
res(end+1,1) = ~isempty(P.bounded.val) && P.bounded.val;

% taking the Minkowski difference of two bounded polytopes should result in bounded polytope
P3 = polytope.generateRandom('Dimension', 2, 'isBounded', true);
P = minkDiff(P,P3);
res(end+1,1) = ~isempty(P.bounded.val) && P.bounded.val;

% taking the Minkowski sum of two bounded polytopes should result in bounded polytope
P = P + P2;
res(end+1,1) = ~isempty(P.bounded.val) && P.bounded.val;

% taking the Minkowski sum of two polytopes, one of which is unbounded
% results in unbounded polytope
P4 = polytope.generateRandom('Dimension',2,'isBounded',false);
P = P + P4;
res(end+1,1) = ~isempty(P.bounded.val) && ~P.bounded.val;

% no change
P = normalizeConstraints(P);
res(end+1,1) = ~isempty(P.bounded.val) && ~P.bounded.val;

% linear map with invertible matrix keeps properties
A = [2 1; -1 0];
P = A*P;
res(end+1,1) = ~isempty(P.bounded.val) && ~P.bounded.val;

% Intersect with bounded set again and lift to higher dim
% -> not bounded anymore
P = P & P2;
P = lift(P,10,[4,5]);
res(end+1,1) = ~isempty(P.bounded.val) && ~P.bounded.val;


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
