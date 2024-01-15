function res = test_polytope_and
% test_polytope_and - unit test function of intersection between a polytope
%    and another set or vector
%
% Syntax:
%    res = test_polytope_and
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

% Authors:       Viktor Kotsev
% Written:       09-May-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% TODOs
% intersection with fully-empty polytope
% equality constraints
% unbounded?


% 1D (convertible to intervals) -------------------------------------------

% empty intersection
A = [1;-1]; b = [2;-1];
P1 = polytope(A,b);
A = [1;-1]; b = [5; -4];
P2 = polytope(A,b);
% compute intersection and check emptiness
P_and = P1 & P2;
res = representsa(P_and,'emptySet');


% intersection at only one point
A = [1;-1]; b = [2;-1];
P1 = polytope(A,b);
A = [1;-1]; b = [3;-2];
P2 = polytope(A,b);
% check boundedness for property test (should be true)
isBounded(P1);
isBounded(P2);
% compute intersection, compare to true result and check boundedness
P_and = P1 & P2;
Ae_true = 1; be_true = 2;
P_true = polytope([],[],Ae_true,be_true);
res(end+1,1) = (P_and == P_true) && ~isempty(P_and.bounded.val) && P_and.bounded.val;

% intersection with polytope without constraints
A = [1; -1]; b = [2; -1];
P1 = polytope(A,b);
A = zeros(0,1); b = zeros(0,0);
P2 = polytope(A,b);
P_and = P1 & P2;
res(end+1,1) = P_and == P1;


% intersection is a full-dimensional polytope
A = [1;-1]; b = [3;-1];
P1 = polytope(A,b);
A = [1;-1]; b = [4;-2];
P2 = polytope(A,b);
% compute intersection and compare to true result
P_and = P1 & P2;
A_true = [1;-1]; b_true = [3;-2];
P_true = polytope(A_true,b_true);
res(end+1,1) = P_and == P_true;

% emptiness property test
A = [1;-1]; b = [2;-3];
P1 = polytope(A,b);
A = [1;-1]; b = [3;-2];
P2 = polytope(A,b);
% determine emptiness of P1
representsa(P1,"emptySet");
% compute intersection and check emptiness
P_and = P1 & P2;
res(end+1,1) = ~isempty(P_and.emptySet.val) && P_and.emptySet.val;

% intersection with numeric (non-empty)
A = [1; -1]; b = [2; -1];
P = polytope(A,b);
p = 1.5;
P_and = P & p;
Ae_true = 1; be_true = 1.5;
P_true = polytope([],[],Ae_true,be_true);
res(end+1,1) = P_and == P_true;

% intersection with numeric (empty)
A = [1; -1]; b = [2; -1];
P = polytope(A,b);
p = 2.5;
P_and = P & p;
res(end+1,1) = representsa(P_and,'emptySet') && dim(P_and) == 1;


% 2D ----------------------------------------------------------------------

% bounded & unbounded -> bounded
A = [-1 1; 0 1; 1 0; -1 -2];
b = [0; 1; 2; 0];
P1 = polytope(A,b);
% determine boundedness
isBounded(P1);
% unbounded polytope
A = [-1 -1]/sqrt(2); b = -sqrt(2);
P2 = polytope(A,b);
% compute intersection, compare to true result (by hand) and check
% boundedness
P_and = P1 & P2;
A_true = [-1/sqrt(2) -1/sqrt(2); 0 1; 1 0]; b_true = [-sqrt(2); 1; 2];
P_true = polytope(A_true,b_true);
res(end+1,1) = P_and == P_true && ~isempty(P_and.bounded.val) && P_and.bounded.val;


% bounded & unbounded -> bounded
A = [-1 1; 0 1; 1 0; -1 -2];
b = [0; 1; 2; 0];
P1 = polytope(A,b);
% determine boundedness
isBounded(P1);
% unbounded polytope
A = [-1 2]; b = 0;
P2 = polytope(A,b);
% compute intersection and compare to true result
P_and = P1 & P2;
A_true = [-1 2; -1 -2; 1 0]; b_true = [0; 0; 2];
P_true = polytope(A_true,b_true);
res(end+1,1) = P_and == P_true;


% bounded & fullspace -> bounded
A = [1 0; -1 1; -1 -1]; b = [1;1;1];
P1 = polytope(A,b);
A = zeros(0,2); b = zeros(0,0);
P2 = polytope(A,b);
P_and = P1 & P2;
res(end+1,1) = P_and == P1;


% bounded & bounded -> bounded
A = [-1 1; 0 1; 1 0; -1 -2];
b = [0; 1; 2; 0];
P1 = polytope(A,b);
% determine boundedness
isBounded(P1);
% another bounded polytope
A = [-1/sqrt(2) -1/sqrt(2); -1 2]; b = [-sqrt(2); 0];
P2 = polytope(A,b);
% compute intersection, compare to true result
P_and = P1 & P2;
A_true = [-1/sqrt(2) -1/sqrt(2); -1 2; 1 0]; b_true = [-sqrt(2); 0; 2];
P_true = polytope(A_true,b_true);
% intersection
res(end+1,1) = P_and == P_true;


% unbounded & unbounded -> bounded
A = [1 1; 1 -1]; b = [2; 2];
P1 = polytope(A,b);
A = [-1 1;-1 -1]; b = [2; 2];
P2 = polytope(A,b);
P_and = P1 & P2;
A_true = [1 1; 1 -1;-1 1;-1 -1]; b_true = [2; 2; 2; 2];
P_true = polytope(A_true,b_true);
res(end+1,1) = P_and == P_true;


% empty & bounded -> empty
A = [-1 -1; -1 1; 1 0];
b = [-2; -2; -1];
P1 = polytope(A,b);
A = [-1 1; 0 1; 1 0; -1 -2];
b = [0; 1; 2; 0];
P2 = polytope(A,b);
% determine emptiness of P1
representsa(P1,"emptySet");
% compute intersection and check emptiness value
P_and = P1 & P2;
res(end+1,1) = ~isempty(P_and.emptySet.val) && P_and.emptySet.val;


% polytope & zonotope
A = [1 0; -1 1; -1 -1]; b = [1;1;1];
P = polytope(A,b);
c = [1;2]; G = [1 1; -1 0];
Z = zonotope(c,G);
P_and = P & Z;
A_true = [1 0; -1 1; -1 -1]; b_true = [1;1;-2];
P_true = polytope(A_true,b_true);
res(end+1,1) = P_and == P_true;

% polytope & interval
A = [1 0; -1 1; -1 -1]; b = [1;1;1];
P = polytope(A,b);
lb = [-1;-2]; ub = [0;0];
I = interval(lb,ub);
P_and = P & I;
A_true = [0 1; 1 0; -1 -1]; b_true = [0; 0; 1];
P_true = polytope(A_true,b_true);
res(end+1,1) = P_and == P_true;


% nD ----------------------------------------------------------------------

% 3D: degenerate & non-degenerate
% square in 3D
A = [0 0 1; 0 0 -1; 1 1 0; -1 1 0; 1 -1 0; -1 -1 0];
b = [1;-1;2;2;2;2];
P1 = polytope(A,b);
% obtain information about degeneracy
isFullDim(P1);
% init non-degenerate polytope
A = [0 0 1; 0 0 -1; 1 1 0; -1 1 0; 1 -1 0; -1 -1 0];
b = [1;-1;3;-1;3;-1];
P2 = polytope(A,b);
% compute intersection
P_and = P1 & P2;
% true solution (by hand)
A_true = [0 0 1; 0 0 -1; 1 1 0; -1 1 0; 1 -1 0; -1 -1 0];
b_true = [1; -1; 2; -1; 2; -1];
P_true = polytope(A_true,b_true);
% compare, also P_ should be degenerate since P1 is degenerate
res(end+1,1) = P_and == P_true && ~isempty(P_and.fullDim.val) && ~P_and.fullDim.val;


% combine result
res = all(res);

% ------------------------------ END OF CODE ------------------------------
