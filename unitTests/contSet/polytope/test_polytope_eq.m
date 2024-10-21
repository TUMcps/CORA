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

% 1D, unbounded + unbounded
A = 3; b = 6;
P1 = polytope(A,b);
A = 1; b = 2;
P2 = polytope(A,b);
assert(P1 == P2);

% 1D, bounded + bounded
A = [1;-1]; b = [2;0.5];
P1 = polytope(A,b);
A = [1;-1]; b = [1;0.5];
P2 = polytope(A,b);
assert(~(P1 == P2));

% 1D, degenerate + degenerate
A = [1; -1]; b = [1; -1];
P1 = polytope(A,b);
Ae = 1; be = 1;
P2 = polytope([],[],Ae,be);
assert(P1 == P2);

% 1D, bounded + point
A = [1; -2]; b = [4; -2];
P = polytope(A,b);
p = 1;
assert(~(P == p));

% 1D, degenerate + point
Ae = 1; be = 3;
P = polytope([],[],Ae,be);
p = 3;
assert(P == p);


% 2D, bounded + bounded
A = [1 1; -1 1; -1 -1; 1 -1] / sqrt(2); b = ones(4,1) / sqrt(2);
P = polytope(A,b);
% same polytope via vertex instantiation
V = [1 0; 0 1; -1 0; 0 -1]';
P_ = polytope(V);
assert(P == P_);

% 2D, polytope is a zonotope
A = [1 1; -1 1; -1 -1; 1 -1] / sqrt(2); b = ones(4,1) / sqrt(2);
P = polytope(A,b);
% init equivalent polytope via conversion from interval
Z = zonotope(zeros(2,1),0.5*[1 1; -1 1]);
P_ = polytope(Z);
assert(P == P_);

% 2D, V-polytope + V-polytope
V1 = [1 0; 0 1; -1 1; -1 -1]';
P1 = polytope(V1);
V2 = [1 0; 0 1; -1 1; -1 -1; 0 -1]';
P2 = polytope(V2);
assert(P1 == P1);
assert(~(P1 == P2));

% 2D, V-polytope + point
V1 = [1 0; -1 1]';
P1 = polytope(V1);
p = [1;0];
assert(~(P1 == p));

% 2D, H-polytope (single point) + point
Ae = [1 0; 0 1]; be = [1; 2];
P1 = polytope([],[],Ae,be);
p1 = [1;2]; p2 = [1;3];
assert(P1 == p1);
assert(~(P1 == p2));

% 2D, H-polytope (full-dimensional) + point
A = [1 0; -1 1; -1 -1]; b = [1; 1; 1];
P1 = polytope(A,b);
p = [1; -1];
assert(~(P1 == p));


% 3D, unit cube
n = 3; A = [eye(n); -eye(n)]; b = [ones(2*n,1)];
P = polytope(A,b);
% init equivalent polytope via conversion from interval
I = interval(-ones(n,1),ones(n,1));
P_ = polytope(I);
assert(P == P_);


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
