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
%    res - btrue/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Viktor Kotsev, Mark Wetzlinger
% Written:       25-April-2022
% Last update:   19-July-2023 (MW, many more cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% 1D, unbounded >= unbounded (same set)
P1 = polytope(3,6);
P2 = polytope(1,2);
res(end+1,1) = contains(P1,P2);
res(end+1,1) = contains(P2,P1);

% 1D, unbounded >= unbounded
P1 = polytope(1,1);
P2 = polytope(1,0);
res(end+1,1) = contains(P1,P2);

% 1D, unbounded >= bounded
P1 = polytope(1,1);
P2 = polytope([1;-1],[0.5;0.75]);
res(end+1,1) = contains(P1,P2);

% 1D, unbounded >= degenerate
P1 = polytope(1,1);
P2 = polytope([],[],1,0);
res(end+1,1) = contains(P1,P2);

% 1D, bounded >= unbounded
P1 = polytope([1; -1],[1; 1]);
P2 = polytope([2; 3],[7; 2]);
res(end+1,1) = ~contains(P1,P2);

% 1D, degenerate >= unbounded
P1 = polytope([1; -1],[1; -1]);
P2 = polytope([2; 3],[7; 2]);
res(end+1,1) = ~contains(P1,P2);

% 1D, degenerate >= degenerate
P1 = polytope([1; -1],[1; -1]);
P2 = polytope([],[],1,0);
res(end+1,1) = ~contains(P1,P2);


% 2D, bounded, non-degenerate >= bounded, non-degenerate
P1 = polytope([1 1; 1 -1; -1 1; 1 1],[1; 1; 1; 1]);
P2 = polytope([1 1; 1 -1; -1 1; 1 1],[2; 2; 2; 2]);
p1 = [0; 0];
p2 = [5; 5];
res(end+1,1) = ~contains(P1,P2) && contains(P2,P1);
res(end+1,1) = contains(P1,p1) && ~contains(P2,p2);

% 2D, unbounded, non-degenerate >= degenerate?
P1 = polytope([1 0; 0 1; -1 0],[2; 2; 2]);
P2 = polytope([1 0; 0 1; -1 0; 0 -1],[1;1;1;-1]);
res(end+1,1) = contains(P1,P2);

% 2D, degenerate >= degenerate?
P1 = polytope([1 0; 0 1; -1 0; 0 -1],[3;1;3;-1]);
P2 = polytope([1 0; 0 1; -1 0; 0 -1],[1;1;1;-1]);
res(end+1,1) = contains(P1,P2);

% 2D, unbounded >= unbounded?
P1 = polytope([1 0; 0 1; -1 0],[2; 2; 2]);
P2 = polytope([1 0; 0 1; -1 0],[1; 1; 1]);
res(end+1,1) = contains(P1,P2);

% 2D, instantiation via vertices
V1 = [1 0; -1 0; 0 1; 0 -1]';
V2 = [2 0; -2 0; 0 2; 0 -1]';
P1 = polytope(V1);
P2 = polytope(V2);
res(end+1,1) = contains(P2,P1);
res(end+1,1) = ~contains(P1,P2);


% 3D, bounded >= degenerate?
P1 = polytope([1 0 0; 0 1 0; 0 0 1; -1 0 0; 0 -1 0; 0 0 -1],ones(6,1));
P2 = polytope([1 0 0; -1 1 0; -1 -1 0],[0.5;0.25;0.25],[0 0 1],0);
res(end+1,1) = contains(P1,P2);

% 3D, bounded >= unbounded?
[M,~,~] = svd(randn(3));
P1 = M * polytope([1 0 0; 0 1 0; 0 0 1; -1 0 0; 0 -1 0; 0 0 -1],ones(6,1));
P2 = polytope([],[],[1 0 0; 0 0 1],[1;2]);
res(end+1,1) = ~contains(P1,P2);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
