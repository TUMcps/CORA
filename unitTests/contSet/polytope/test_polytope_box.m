function res = test_polytope_box
% test_polytope_box - unit test function of box
%
% Syntax:
%    res = test_polytope_box
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
% See also: -

% Authors:       Viktor Kotsev, Mark Wetzlinger
% Written:       16-May-2022
% Last update:   27-July-2023 (MW, more tests)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% 1D ----------------------------------------------------------------------
% all box conversions are exact since 1D polytopes are intervals

% bounded
A = [1;-1]; b = [1;1];
P = polytope(A,b);
P_box = box(P);
res(end+1,1) = P_box == P;

% unbounded
A = 3; b = 6;
P = polytope(A,b);
P_box = box(P);
res(end+1,1) = P_box == P;

% fully empty (fullspace)
A = zeros(0,1); b = zeros(0,0);
P = polytope(A,b);
P_box = box(P);
res(end+1,1) = P_box == P;

% single point (inequalities)
A = [2;-1]; b = [1;-0.5];
P = polytope(A,b);
P_box = box(P);
res(end+1,1) = P_box == P;

% single point (equality)
A = zeros(0,1); b = zeros(0,0); Ae = 1; be = 3;
P = polytope(A,b,Ae,be);
P_box = box(P);
res(end+1,1) = P_box == P;

% empty (inequalities)
A = [2;-1]; b = [1;-4];
P = polytope(A,b);
P_box = box(P);
res(end+1,1) = representsa_(P_box,'emptySet');

% empty (equalities)
A = zeros(0,1); b = zeros(0,0); Ae = [2;1]; be = [5;2];
P = polytope(A,b,Ae,be);
P_box = box(P);
res(end+1,1) = representsa_(P_box,'emptySet');


% 2D ----------------------------------------------------------------------

% bounded 
A = [1 1; 1 -1;-1 1;-1 -1];
b = [2; 2; 2; 2];
P = polytope(A,b);
% compute box and compare to true result
P_box = box(P);
A_true = [1 0; 0 1;-1 0; 0 -1]; b_true = [2; 2; 2; 2];
P_true = polytope(A_true,b_true);
res(end+1,1) = P_box == P_true;

% bounded and degenerate
A = [1 0; -2 0]; b = [2; 5];
Ae = [-1 1]; be = 2;
P = polytope(A,b,Ae,be);
% compute box and compare to true result
P_box = box(P);
A_true = [1 0; 0 1;-1 0; 0 -1]; b_true = [2; 4; 2.5; 0.5];
P_true = polytope(A_true,b_true);
res(end+1,1) = P_box == P_true;

% unbounded
A = [1 1; 1 -1]; b = [2; 2];
P = polytope(A,b);
P_box = box(P);
% compute box and compare to true result
A_true = [1 0]; b_true = 2;
P_true = polytope(A_true,b_true);
res(end+1,1) = P_true == P_box;

% unbounded and degenerate
Ae = [1 1]; be = 2;
P = polytope([],[],Ae,be);
P_box = box(P);
res(end+1,1) = isemptyobject(P_box);

% unbounded and degenerate (axis-aligned)
Ae = [0 3]; be = 2;
P = polytope([],[],Ae,be);
P_box = box(P);
% compute box and compare to true result
Ae_true = [0 1]; be_true = 2/3;
P_true = polytope([],[],Ae_true,be_true);
res(end+1,1) = P_true == P_box;

% bounded by vertex instantiation
V = [1 0; -1 0; 0 1; 0 -1]';
P = polytope(V);
% compute box and compare to true result
P_box = box(P);
A_true = [1 0;-1 0;0 1;0 -1]; b_true = [1;1;1;1];
P_true = polytope(A_true,b_true);
res(end+1,1) = P_true == P_box;


% nD ----------------------------------------------------------------------

% 3D, degenerate
A = [0 0 1; 0 0 -1; 1 1 0; -1 1 0; 1 -1 0; -1 -1 0];
b = [1;-1;2;2;2;2];
P = polytope(A,b);
% compute box and compute to true result
P_box = box(P);
A_true = [0 0 1; 0 0 -1; 1 0 0; -1 0 0; 0 1 0; 0 -1 0]; b_true = [1;-1;2;2;2;2];
P_true = polytope(A_true,b_true);
res(end+1,1) = P_true == P_box;

% 5D, single point
[Ae,~,~] = svd(randn(5)); be = ones(5,1);
P = polytope([],[],Ae',be);
% compute box and check if its a point
P_box = box(P);
res(end+1,1) = representsa(P_box,'point');


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
