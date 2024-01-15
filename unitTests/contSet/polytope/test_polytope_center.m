function res = test_polytope_center
% test_polytope_center - unit test function of center
%
% Syntax:
%    res = test_polytope_center
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
% Written:       25-April-2022
% Last update:   27-July-2023 (MW, more cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% 1D, only inequalities, bounded
A = [2;-1]; b = [6;1];
P = polytope(A,b);
c = center(P);
c_true = 1;
res(end+1,1) = all(withinTol(c,c_true));

% 1D, only equalities, single point
Ae = 3; be = 5;
P = polytope([],[],Ae,be);
c = center(P);
c_true = 5/3;
res(end+1,1) = all(withinTol(c,c_true));

% 1D, only inequalities, unbounded
A = [3;2;4]; b = [5;2;-3];
P = polytope(A,b);
c = center(P);
res(end+1,1) = all(isnan(c));

% 1D, only inequalities, empty
Ae = [1;4]; be = [2;-5];
P = polytope([],[],Ae,be);
c = center(P);
res(end+1,1) = isempty(c);

% 1D, inequalities and equalities, empty
A = [1;-4]; b = [4;-2]; Ae = 5; be = 100;
P = polytope(A,b,Ae,be);
c = center(P);
res(end+1,1) = isempty(c);

% 1D, fully empty
A = zeros(0,1); b = zeros(0,0);
P = polytope(A,b);
c = center(P);
res(end+1,1) = isscalar(c) && isnan(c);


% 2D, only inequalities, bounded
A = [1 1; -1 1; 1 -1; -1 -1]; b = ones(4,1);
P = polytope(A,b);
c = center(P);
c_true = [0; 0];
res(end+1,1) = all(withinTol(c,c_true));

% 2D, only inequalities, empty
A = [1 0; -1 0]; b = [-1; -1];
P = polytope(A,b);
c = center(P);
res(end+1,1) = isempty(c);

% 2D, only equalities, empty
Ae = [1 0; 0 1; 0 1]; be = [1 -1 0];
P = polytope([],[],Ae,be);
c = center(P);
res(end+1,1) = isempty(c);

% 2D, only equalities, single point
Ae = [1 0; 0 1]; be = [0;0];
P = polytope([],[],Ae,be);
c = center(P);
c_true = [0; 0];
res(end+1,1) = all(withinTol(c,c_true));

% 2D, inequalities and equalities, unbounded
A = [1 0]; b = 1; Ae = [0 1]; be = 1;
P = polytope(A,b,Ae,be);
c = center(P);
res(end+1,1) = all(isnan(c));

% 2D, inequalities and equalities, bounded
A = [1 0; -1 0]; b = [1; 1]; Ae = [0 1]; be = 1;
P = polytope(A,b,Ae,be);
c = center(P);
c_true = [0; 1];
res(end+1,1) = all(withinTol(c,c_true));

% 2D, inequalities and equalities, empty
A = [2 1; -1 2; 0 -1]; b = ones(3,1); Ae = [1 1]; be = 10;
P = polytope(A,b,Ae,be);
c = center(P);
res(end+1,1) = isempty(c);

% 2D, fully empty
Ae = zeros(0,2); be = zeros(0,0);
P = polytope([],[],Ae,be);
c = center(P);
res(end+1,1) = all(size(c) == [2,1]) && all(isnan(c));


% 3D, only inequalities, bounded
A = [0 1 0; 0 0 1; 0 -1 0; 0 0 -1; 1 0 0; -1 0 0]; b = ones(6,1);
P = polytope(A,b);
c = center(P);
c_true = [0; 0; 0];
res(end+1,1) = all(withinTol(c,c_true));

% 3D, inequalities and equalities, bounded, degenerate
A = [0 1 0; 0 0 1; 0 -1 0; 0 0 -1]; b = ones(4,1); Ae = [1 0 0]; be = 2;
P = polytope(A,b,Ae,be);
c = center(P);
c_true = [2; 0; 0];
res(end+1,1) = all(withinTol(c,c_true));

% 3D, only inequalities, unbounded
A = [1 0 0; 0 1 0]; b = [0; 0];
P = polytope(A,b);
c = center(P);
res(end+1,1) = all(isnan(c));

% 3D, only equalities, unbounded
Ae = [1 0 0; 0 1 0]; be = [0;0];
P = polytope([],[],Ae,be);
c = center(P);
res(end+1,1) = all(isnan(c));

% 3D, inequalities and equalities, unbounded
A = [1 1 0; -1 0 0]; b = [1; 1]; Ae = [0 1 1]; be = 1;
P = polytope(A,b,Ae,be);
c = center(P);
res(end+1,1) = all(isnan(c));

% 3D, single point
Ae = [1 2 -1; 0 1 1; -1 2 1]; be = [1;1;1];
P = polytope([],[],Ae,be);
c = center(P);
c_true = [0.5;0.5;0.5];
res(end+1,1) = all(withinTol(c,c_true));


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
