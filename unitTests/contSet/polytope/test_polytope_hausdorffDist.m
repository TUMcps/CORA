function res = test_polytope_hausdorffDist
% test_polytope_hausdorffDist - unit test function of the computation of
%    the Hausdorff distance
%
% Syntax:
%    res = test_polytope_hausdorffDist
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
% Written:       04-December-2022
% Last update:   16-December-2023 (MW, more tests)
%                18-December-2023 (MW, more tests)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);


% 1D, unbounded, point
A = 1; b = 1;
P = polytope(A,b);
p = 2;
dH = hausdorffDist(P,p);
dH_true = Inf;
res(end+1,1) = withinTol(dH,dH_true);

% 1D, degenerate, point
Ae = 1; be = 3;
P = polytope([],[],Ae,be);
p = 2;
dH = hausdorffDist(P,p);
dH_true = 1;
res(end+1,1) = withinTol(dH,dH_true);

% 1D, empty, point
Ae = [1;1]; be = [1;2];
P = polytope([],[],Ae,be);
p = 2;
dH = hausdorffDist(P,p);
dH_true = Inf;
res(end+1,1) = withinTol(dH,dH_true);

% 1D, empty, empty point
Ae = [1;1]; be = [1;2];
P = polytope([],[],Ae,be);
p = zeros(1,0);
dH = hausdorffDist(P,p);
dH_true = 0;
res(end+1,1) = withinTol(dH,dH_true);

% 1D, bounded, bounded
A = [1; -1]; b = [3; -2];
P1 = polytope(A,b);
A = [2; -6]; b = [3; 2];
P2 = polytope(A,b);
dH = hausdorffDist(P1,P2);
dH_true = 2 + 1/3;
res(end+1,1) = withinTol(dH,dH_true);

% 1D, bounded, empty
A = [1; -1]; b = [3; -2];
P1 = polytope(A,b);
Ae = [1;1]; be = [1;2];
P2 = polytope([],[],Ae,be);
dH = hausdorffDist(P1,P2);
dH_true = Inf;
res(end+1,1) = withinTol(dH,dH_true);

% 1D, empty, empty
A = [1; -1]; b = [3; -4];
P1 = polytope(A,b);
Ae = [1;1]; be = [1;2];
P2 = polytope([],[],Ae,be);
dH = hausdorffDist(P1,P2);
dH_true = 0;
res(end+1,1) = withinTol(dH,dH_true);

% 1D, unbounded, empty
A = -1; b = 3;
P1 = polytope(A,b);
Ae = [1;1]; be = [1;2];
P2 = polytope([],[],Ae,be);
dH = hausdorffDist(P1,P2);
dH_true = Inf;
res(end+1,1) = withinTol(dH,dH_true);

% 1D, bounded, unbounded
A = [1; -1]; b = [3; -2];
P1 = polytope(A,b);
A = 2; b = 3;
P2 = polytope(A,b);
dH = hausdorffDist(P1,P2);
dH_true = Inf;
res(end+1,1) = withinTol(dH,dH_true);

% 1D, unbounded, unbounded (non-identical)
A = 1; b = 2;
P1 = polytope(A,b);
A = -1; b = 3;
P2 = polytope(A,b);
dH = hausdorffDist(P1,P2);
dH_true = Inf;
res(end+1,1) = withinTol(dH,dH_true);

% 1D, unbounded, unbounded (identical)
A = 1; b = 2;
P1 = polytope(A,b);
dH = hausdorffDist(P1,P1);
dH_true = 0;
res(end+1,1) = withinTol(dH,dH_true);

% 1D, fully empty, unbounded
A = zeros(0,1); b = zeros(0,0);
P1 = polytope(A,b);
A = 1; b = 1;
P2 = polytope(A,b);
dH = hausdorffDist(P1,P2);
dH_true = Inf;
res(end+1,1) = withinTol(dH,dH_true);


% 2D, bounded, empty
A = [1 0; -1 1; -1 1]; b = ones(3,1);
P1 = polytope(A,b);
A = [1 1]; b = 1; Ae = [1 1]; be = 2;
P2 = polytope(A,b,Ae,be);
dH = hausdorffDist(P1,P2);
dH_true = Inf;
res(end+1,1) = withinTol(dH,dH_true);

% 2D, unbounded, empty
A = [1 0; -1 1]; b = ones(2,1);
P1 = polytope(A,b);
A = [1 1]; b = 1; Ae = [1 1]; be = 2;
P2 = polytope(A,b,Ae,be);
dH = hausdorffDist(P1,P2);
dH_true = Inf;
res(end+1,1) = withinTol(dH,dH_true);

% 2D, box, point
A = [1 0; 0 1; -1 0; 0 -1]; b = ones(4,1);
P = polytope(A,b);
% first point
p = [2;0];
dH = hausdorffDist(P,p);
dH_true = 1;
res(end+1,1) = withinTol(dH,dH_true);
% second point
p = [2; 2];
dH = hausdorffDist(P,p);
dH_true = sqrt(2);
res(end+1,1) = withinTol(dH,dH_true);
% point in the set (dist = 0)
p = [0;0];
dH = hausdorffDist(P,p);
dH_true = 0;
res(end+1,1) = withinTol(dH,dH_true);

% 2D, rotated box, point
A = [1 1; -1 1; -1 -1; 1 -1]; b = ones(4,1);
P = polytope(A,b);
p = [1;1];
dH = hausdorffDist(P,p);
dH_true = sqrt(2)/2;
res(end+1,1) = withinTol(dH,dH_true);
p = [1 1; 1 0; -0.5 1; 0 -0.5]';
dH = hausdorffDist(P,p);
dH_true = sqrt(2)/2;
res(end+1,1) = withinTol(dH,dH_true);

% 2D, degenerate, point
P1 = polytope([1 0; 0 1;-1 0;0 -1],[1;1;1;-1]);
P2 = polytope([1 0;0 1;-1 0;0 -1],[6;1;-4;-1]);
dH = hausdorffDist(P1,P2);
dH_true = 5;
res(end+1,1) = withinTol(dH,dH_true);

% 2D, unbounded, point
P = polytope([1 0; 0 1;-1 0],[1;1;1]);
p = [2;-100];
dH = hausdorffDist(P,p);
dH_true = Inf;
res(end+1,1) = withinTol(dH,dH_true,1e-4);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
