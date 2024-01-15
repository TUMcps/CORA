function res = test_polytope_eq
% test_polytope_eq - unit test function of eq
%
% Syntax:
%    res = test_polytope_eq
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
% Written:       03-December-2022
% Last update:   16-December-2023 (MW, more tests)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% 1D, unbounded + unbounded
A = 3; b = 6;
P1 = polytope(A,b);
A = 1; b = 2;
P2 = polytope(A,b);
res(end+1,1) = P1 == P2;

% 1D, bounded + bounded
A = [1;-1]; b = [2;0.5];
P1 = polytope(A,b);
A = [1;-1]; b = [1;0.5];
P2 = polytope(A,b);
res(end+1,1) = ~(P1 == P2);

% 1D, degenerate + degenerate
A = [1; -1]; b = [1; -1];
P1 = polytope(A,b);
Ae = 1; be = 1;
P2 = polytope([],[],Ae,be);
res(end+1,1) = P1 == P2;

% 1D, bounded + point
A = [1; -2]; b = [4; -2];
P = polytope(A,b);
p = 1;
res(end+1,1) = ~(P == p);

% 1D, degenerate + point
Ae = 1; be = 3;
P = polytope([],[],Ae,be);
p = 3;
res(end+1,1) = P == p;


% 2D, bounded + bounded
A = [1 1; -1 1; -1 -1; 1 -1] / sqrt(2); b = ones(4,1) / sqrt(2);
P = polytope(A,b);
% same polytope via vertex instantiation
V = [1 0; 0 1; -1 0; 0 -1]';
P_ = polytope(V);
res(end+1,1) = P == P_;

% 2D, polytope is a zonotope
A = [1 1; -1 1; -1 -1; 1 -1] / sqrt(2); b = ones(4,1) / sqrt(2);
P = polytope(A,b);
% init equivalent polytope via conversion from interval
Z = zonotope(zeros(2,1),0.5*[1 1; -1 1]);
P_ = polytope(Z);
res(end+1,1) = P == P_;


% 3D, unit cube
n = 3; A = [eye(n); -eye(n)]; b = [ones(2*n,1)];
P = polytope(A,b);
% init equivalent polytope via conversion from interval
I = interval(-ones(n,1),ones(n,1));
P_ = polytope(I);
res(end+1,1) = P == P_;


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
