function res = test_polytope_dim
% test_polytope_dim - unit test function of dim
%
% Syntax:
%    res = test_polytope_dim
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

% empty polytope
P_empty = polytope.empty(2);
res = dim(P_empty) == 2;


% 1D, bounded, degenerate
A = [1;-1]; b = [1;-1];
P = polytope(A,b);
res(end+1,1) = dim(P) == 1;

% 1D, empty
A = [1;-1]; b = [2;-3];
P = polytope(A,b);
res(end+1,1) = dim(P) == 1;

% 1D, unbounded
A = 1; b = 1;
P = polytope(A,b);
res(end+1,1) = dim(P) == 1;

% 1D, single point, only equality constraints
Ae = 1; be = 3;
P = polytope([],[],Ae,be);
res(end+1,1) = dim(P) == 1;


% 2D, bounded, non-degenerate
A = [-1 -1; 1 0;-1 0; 0 1; 0 -1];
b = [2;3;2;3;2];
P = polytope(A,b);
res(end+1,1) = dim(P) == 2;

% 2D, bounded, degenerate
A = [1 0; -1 0]; b = [1; 1]; Ae = [0 1]; be = 1;
P = polytope(A,b,Ae,be);
res(end+1,1) = dim(P) == 2;

% 2D, unbounded, non-degenerate 
A = [1 0; 0 1; -1 0]; b = [2; 2; 2];
P = polytope(A,b);
res(end+1,1) = dim(P) == 2;

% 2D, unbounded, degenerate
A = [1 0; -1 0]; b = [1; -1];
P = polytope(A,b);
res(end+1,1) = dim(P) == 2;

% 2D, empty, only equality constraints
Ae = [0 1; 1 0; 1 1]; be = [1;1;1];
P = polytope([],[],Ae,be);
res(end+1,1) = dim(P) == 2;


% 3D, fully empty
A = zeros(0,3); b = zeros(0,0);
P = polytope(A,b);
res(end+1,1) = dim(P) == 3;
Ae = zeros(0,3); be = zeros(0,0);
P = polytope([],[],Ae,be);
res(end+1,1) = dim(P) == 3;


% 4D, box
A = [eye(4); -eye(4)]; b = ones(8,1);
P = polytope(A,b);
res(end+1,1) = dim(P) == 4;


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
