function res = test_polytope_contains()
% test_polytope_contains - unit test function of contains
%
% Syntax:
%    res = test_polytope_contains
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
% Written:       25-April-2022
% Last update:   19-July-2023 (MW, many more cases)
%                10-July-2024 (MW, test containment of point clouds)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% 1D, unbounded (fully empty) >= unbounded (fully empty)
A = zeros(0,1); b = zeros(0,0);
P1 = polytope(A,b);
res(end+1,1) = contains(P1,P1);

% 1D, unbounded (fully empty) >= unbounded
A = zeros(0,1); b = zeros(0,0);
P1 = polytope(A,b);
A = 1; b = 1;
P2 = polytope(A,b);
res(end+1,1) = contains(P1,P2);

% 1D, unbounded >= unbounded (same set)
A = 3; b = 6;
P1 = polytope(A,b);
A = 1; b = 2;
P2 = polytope(A,b);
res(end+1,1) = contains(P1,P2);
res(end+1,1) = contains(P2,P1);

% 1D, unbounded >= unbounded
A = 1; b = 1;
P1 = polytope(A,b);
A = 1; b = 0;
P2 = polytope(A,b);
res(end+1,1) = contains(P1,P2);

% 1D, unbounded >= bounded
A = 1; b = 1;
P1 = polytope(A,b);
A = [1;-1]; b = [0.5;0.75];
P2 = polytope(A,b);
res(end+1,1) = contains(P1,P2);

% 1D, unbounded >= degenerate
A = 1; b = 1;
P1 = polytope(A,b);
Ae = 1; be = 0;
P2 = polytope([],[],Ae,be);
res(end+1,1) = contains(P1,P2);

% 1D, bounded >= unbounded
A = [1; -1]; b = [1; 1];
P1 = polytope(A,b);
A = [2; 3]; b = [7; 2];
P2 = polytope(A,b);
res(end+1,1) = ~contains(P1,P2);

% 1D, degenerate >= unbounded
A = [1; -1]; b = [1; -1];
P1 = polytope(A,b);
A = [2; 3]; b = [7; 2];
P2 = polytope(A,b);
res(end+1,1) = ~contains(P1,P2);

% 1D, degenerate >= degenerate
A = [1; -1]; b = [1; -1];
P1 = polytope(A,b);
Ae = 1; be = 0;
P2 = polytope([],[],Ae,be);
res(end+1,1) = ~contains(P1,P2);

% 1D, unbounded >= point cloud
A = [1; 2]; b = [1; 1];
P1 = polytope(A,b);
S = [-4, -2, 0];
res(end+1,1) = all(contains(P1,S));

% 1D, bounded >= point cloud
A = [1; -1]; b = [1; 1];
P1 = polytope(A,b);
S = [-3, 0.5, 0];
res(end+1,1) = ~all(contains(P1,S));

% 1D, degenerate >= point cloud
Ae = 1; be = -3;
P1 = polytope([],[],Ae,be);
S = -3;
res(end+1,1) = contains(P1,S);

% 1D, degenerate >= point cloud
Ae = 1; be = -3;
P1 = polytope([],[],Ae,be);
S = [-3, -2];
res(end+1,1) = ~all(contains(P1,S));

% 1D, polytope >= non-polytope
A = [1;-1]; b = [2;3];
P1 = polytope(A,b);
P2 = conHyperplane(1,1);
res(end+1,1) = contains(P1,P2);


% 2D, bounded, non-degenerate >= bounded, non-degenerate
A = [1 1; 1 -1; -1 1; 1 1]; b = [1; 1; 1; 1];
P1 = polytope(A,b);
A = [1 1; 1 -1; -1 1; 1 1]; b = [2; 2; 2; 2];
P2 = polytope(A,b);
p1 = [0; 0]; p2 = [5; 5];
res(end+1,1) = ~contains(P1,P2) && contains(P2,P1);
res(end+1,1) = contains(P1,p1) && ~contains(P2,p2);

% 2D, unbounded, non-degenerate >= degenerate?
A = [1 0; 0 1; -1 0]; b = [2; 2; 2];
P1 = polytope(A,b);
A = [1 0; 0 1; -1 0; 0 -1]; b = [1;1;1;-1];
P2 = polytope(A,b);
res(end+1,1) = contains(P1,P2);

% 2D, degenerate >= degenerate?
A = [1 0; 0 1; -1 0; 0 -1]; b = [3;1;3;-1];
P1 = polytope(A,b);
A = [1 0; 0 1; -1 0; 0 -1]; b = [1;1;1;-1];
P2 = polytope(A,b);
res(end+1,1) = contains(P1,P2);

% 2D, unbounded >= unbounded?
A = [1 0; 0 1; -1 0]; b = [2; 2; 2];
P1 = polytope(A,b);
A = [1 0; 0 1; -1 0]; b = [1; 1; 1];
P2 = polytope(A,b);
res(end+1,1) = contains(P1,P2);

% 2D, instantiation via vertices
V = [1 0; -1 0; 0 1; 0 -1]';
P1 = polytope(V);
V = [2 0; -2 0; 0 2; 0 -1]';
P2 = polytope(V);
res(end+1,1) = contains(P2,P1);
res(end+1,1) = ~contains(P1,P2);

% 2D, V-polytope >= H-polyhedron?
V = [1 0; -1 1; -1 -1]';
P1 = polytope(V);
A = [1 0]; b = 0;
P2 = polytope(A,b);
res(end+1,1) = ~contains(P1,P2);

% 2D, V-polytope >= point cloud?
V = [1 0; -1 1; -1 -1]';
P1 = polytope(V);
S = [0 0; 0.5 0.1; -0.8 -0.6; -0.5 0.5]';
res(end+1,1) = all(contains(P1,S));

% 2D, H-polyhedron >= point cloud?
A = [1 0; 0 1]; b = [1; 0];
P1 = polytope(A,b);
S = [-2 -1; 0 0; 0 -4]';
res(end+1,1) = all(contains(P1,S));


% 3D, bounded >= degenerate?
A = [1 0 0; 0 1 0; 0 0 1; -1 0 0; 0 -1 0; 0 0 -1]; b = ones(6,1);
P1 = polytope(A,b);
A = [1 0 0; -1 1 0; -1 -1 0]; b = [0.5;0.25;0.25]; Ae = [0 0 1]; be = 0;
P2 = polytope(A,b,Ae,be);
res(end+1,1) = contains(P1,P2);

% 3D, bounded >= unbounded?
[M,~,~] = svd(randn(3));
A = [1 0 0; 0 1 0; 0 0 1; -1 0 0; 0 -1 0; 0 0 -1]; b = ones(6,1);
P1 = M * polytope(A,b);
Ae = [1 0 0; 0 0 1]; be = [1;2];
P2 = polytope([],[],Ae,be);
res(end+1,1) = ~contains(P1,P2);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
