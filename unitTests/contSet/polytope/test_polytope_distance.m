function res = test_polytope_distance
% test_polytope_distance - unit test function of shortest distance
%    computation
%
% Syntax:
%    res = test_polytope_distance
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
% See also: test_polytope_hausdorffDist

% Authors:       Viktor Kotsev, Mark Wetzlinger
% Written:       07-November-2022
% Last update:   18-December-2023 (MW, more tests)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);
tol = 1e-8;

% 1D, empty, bounded
A = [1; -1]; b = [2; -5];
P1 = polytope(A,b);
Ae = 1; be = 2;
P2 = polytope([],[],Ae,be);
dist = distance(P1,P2);
dist_true = Inf;
res(end+1,1) = withinTol(dist,dist_true);

% 1D, empty, unbounded
A = [1; -1]; b = [2; -5];
P1 = polytope(A,b);
A = 1; b = 2;
P2 = polytope(A,b);
dist = distance(P1,P2);
dist_true = Inf;
res(end+1,1) = withinTol(dist,dist_true);

% 1D, bounded, bounded
A = [1; -1]; b = [1; 1];
P1 = polytope(A,b);
A = [1; -1]; b = [4;-2];
P2 = polytope(A,b);
dist = distance(P1,P2);
dist_true = 1;
res(end+1,1) = withinTol(dist,dist_true);

% 1D, unbounded, bounded
A = 1; b = 1;
P1 = polytope(A,b);
A = [1; -1]; b = [4;-2];
P2 = polytope(A,b);
dist = distance(P1,P2);
dist_true = 1;
res(end+1,1) = withinTol(dist,dist_true);

% 1D, unbounded, unbounded (non-intersecting)
A = 1; b = 1;
P1 = polytope(A,b);
A = -1; b = -3;
P2 = polytope(A,b);
dist = distance(P1,P2);
dist_true = 2;
res(end+1,1) = withinTol(dist,dist_true);

% 1D, unbounded, unbounded (intersecting)
A = 1; b = 1;
P1 = polytope(A,b);
A = -1; b = 3;
P2 = polytope(A,b);
dist = distance(P1,P2);
dist_true = 0;
res(end+1,1) = withinTol(dist,dist_true);

% 1D, fullspace, unbounded
A = zeros(0,1); b = zeros(0,0);
P1 = polytope(A,b);
A = 1; b = 2;
P2 = polytope(A,b);
dist = distance(P1,P2);
dist_true = 0;
res(end+1,1) = withinTol(dist,dist_true);

% 2D, empty, bounded
A = [1 0; -1 1; -1 -1; 0 1]; b = [1;1;1;-4];
P1 = polytope(A,b);
A = [1 0; -1 1; -1 -1]; b = [1;1;1];
P2 = polytope(A,b);
dist = distance(P1,P2);
dist_true = Inf;
res(end+1,1) = withinTol(dist,dist_true);

% 2D, bounded, bounded
A = [1 0;-1 0; 0 1; 0 -1]; b = [6;-5;1;1];
P1 = polytope(A,b);
A = [1 0;-1 0; 0 1; 0 -1]; b = [1;1;1;1];
P2 = polytope(A,b);
dist = distance(P1,P2);
dist_true = 4;
res(end+1,1) = withinTol(dist,dist_true);

% 2D, bounded & degenerate, bounded & degenerate
A = [1 0; 0 1; -1 0; 0 -1]; b = [1;1;1;-1];
P1 = polytope(A,b);
A = [1 0; 0 1; -1 0; 0 -1]; b = [6;3;-4;-3];
P2 = polytope(A,b);
dist = distance(P1,P2);
dist_true = 3.605551275;
res(end+1,1) = withinTol(dist,dist_true,tol);

% 2D, bounded, unbounded (non-intersecting)
A = [1 0; -1 1; -1 -1]; b = [1;1;1];
P1 = polytope(A,b);
A = [-1 0]; b = -5;
P2 = polytope(A,b);
dist = distance(P1,P2);
dist_true = 4;
res(end+1,1) = withinTol(dist,dist_true);

% 2D, bounded, unbounded (intersecting)
A = [1 0; -1 1; -1 -1]; b = [1;1;1];
P1 = polytope(A,b);
A = [0 1]; b = 0;
P2 = polytope(A,b);
dist = distance(P1,P2);
dist_true = 0;
res(end+1,1) = withinTol(dist,dist_true);

% 2D, unbounded, unbounded (non-intersecting)
A = [1 0;0 1;-1 0]; b = [1;1;1];
P1 = polytope(A,b);
A = [1 0;0 -1;-1 0]; b = [6;-3;-4];
P2 = polytope(A,b);
dist = distance(P1,P2);
dist_true = 3.605551275;
res(end+1,1) = withinTol(dist,dist_true,tol);

% 2D, unbounded, unbounded (intersecting)
A = [1 1; -1 1]; b = [1;1];
P1 = polytope(A,b);
A = [0 -1]; b = -0.5;
P2 = polytope(A,b);
dist = distance(P1,P2);
dist_true = 0;
res(end+1,1) = withinTol(dist,dist_true,tol);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
