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
% See also: -

% Authors:       Viktor Kotsev
% Written:       09-May-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty polytope
P_empty = polytope();
res = dim(P_empty) == 0;


% 1D, bounded, non-degenerate
P = polytope([1; -1],[1; -1]);
res(end+1,1) = dim(P) == 1;

% 1D, empty
P = polytope([1;-1],[2;-3]);
res(end+1,1) = dim(P) == 1;

% 1D, unbounded
P = polytope(1,1);
res(end+1,1) = dim(P) == 1;

% 1D, single point, only equality constraints
P = polytope([],[],1,3);
res(end+1,1) = dim(P) == 1;


% 2D, bounded, non-degenerate
A = [-1 -1; 1 0;-1 0; 0 1; 0 -1];
b = [2;3;2;3;2];
P = polytope(A,b);
res(end+1,1) = dim(P) == 2;

% 2D, bounded, degenerate
P = polytope([1 0; -1 0],[1; 1],[0 1], 1);
res(end+1,1) = dim(P) == 2;

% 2D, unbounded, non-degenerate 
P = polytope([1 0; 0 1; -1 0],[2; 2; 2]);
res(end+1,1) = dim(P) == 2;

% 2D, unbounded, degenerate
P = polytope([1 0; -1 0],[1; -1]);
res(end+1,1) = dim(P) == 2;

% 2D, empty, only equality constraints
P = polytope([],[],[0 1; 1 0; 1 1],[1;1;1]);
res(end+1,1) = dim(P) == 2;


% 3D, fully empty
P = polytope(zeros(0,3),[]);
res(end+1,1) = dim(P) == 3;
P = polytope([],[],zeros(0,3),[]);
res(end+1,1) = dim(P) == 3;


% 4D, box
P = polytope([eye(4); -eye(4)],ones(8,1));
res(end+1,1) = dim(P) == 4;


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
