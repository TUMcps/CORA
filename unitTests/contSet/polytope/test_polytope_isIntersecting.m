function res = test_polytope_isIntersecting
% test_polytope_isIntersecting - unit test function of intersection check
%
% Syntax:
%    res = test_polytope_isIntersecting
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
% Written:       30-November-2022
% Last update:   27-July-2023 (MW, more tests)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% 1D, unbounded & unbounded
A = 1; b = 1;
P1 = polytope(A,b);
A = -1; b = 5;
P2 = polytope(A,b);
res(end+1,1) = isIntersecting(P1,P2);

% 1D, unbounded & degenerate
A = 1; b = 1;
P1 = polytope(A,b);
Ae = 1; be = -5;
P2 = polytope([],[],Ae,be);
res(end+1,1) = isIntersecting(P1,P2);

% 1D, bounded & unbounded, intersection in a single point
A = [1;-1]; b = [1;1];
P1 = polytope(A,b);
Ae = -1; be = 1;
P2 = polytope([],[],Ae,be);
res(end+1,1) = isIntersecting(P1,P2);

% 1D, unbounded & point
A = 1; b = 1;
P1 = polytope(A,b);
p = 5;
res(end+1,1) = ~isIntersecting(P1,p);

% 1D, fully empty & point
A = zeros(0,1); b = zeros(0,0);
P1 = polytope(A,b);
p = 1;
res(end+1,1) = isIntersecting(P1,p);


% 2D, fully empty & unbounded
A = zeros(0,2); b = zeros(0,0);
P1 = polytope(A,b);
A = [1 0]; b = 1;
P2 = polytope(A,b);
res(end+1,1) = isIntersecting(P1,P2);
res(end+1,1) = isIntersecting(P2,P1);

% 2D, bounded polytopes in all quadrants
V1 = [1 1; 4 0; 3 3; 2 4]';
P1 = polytope(V1);
V2 = [-2 1; -4 2; -3 4; -1 2]';
P2 = polytope(V2);
V3 = [-1 -2; -4 -1; -3 -4; -2 -3]';
P3 = polytope(V3);
V4 = [1 -1; 2 -4; 6 -5; 5 -2]';
P4 = polytope(V4);
% no combination should intersect
res(end+1,1) = ~isIntersecting(P1,P2);
res(end+1,1) = ~isIntersecting(P1,P3);
res(end+1,1) = ~isIntersecting(P1,P4);
res(end+1,1) = ~isIntersecting(P2,P3);
res(end+1,1) = ~isIntersecting(P2,P4);
res(end+1,1) = ~isIntersecting(P3,P4);

% 2D, bounded & point cloud (some contained, some not)
A = [1 0; -1 1; -1 -1]; b = [1;1;1];
P = polytope(A,b);
V = [0.5 0; 0 1.5; -0.5 -1; 0 -1.5; 1 -1]';
res(end+1,1) = isIntersecting(P,V);

% 2D, bounded & bounded (both contain the origin)
V = [2 2; 3 -1; -1 0; 0 3; 1 3]';
P = polytope(V);
% multiply by -1 -> resulting polytope also contains the origin
V_ = -1 * V;
P_ = polytope(V_);
res(end+1,1) = isIntersecting(P,P_);

% 2D, unbounded & unbounded
A1 = [1 0; -1 0; 0 1]; b1 = [1; 1; 5];
P1 = polytope(A1,b1);
A2 = [1 0; 0 1; 0 -1]; b2 = [5; 1; 1];
P2 = polytope(A2, b2);
res(end+1,1) = isIntersecting(P1,P2);

% 2D, unbounded & unbounded, intersection in single point
A1 = [1 0; -1 0; 0 1]; b1 = [1; -1; 5];
P1 = polytope(A1,b1);
A2 = [1 0; 0 1; 0 -1]; b2 = [5; 1; -1];
P2 = polytope(A2, b2);
res(end+1,1) = isIntersecting(P1, P2);

% 2D, unbounded & unbounded
A1 = [1 0; -1 0; 0 1]; b1 = [1; -1; 5];
P1 = polytope(A1,b1);
A2 = [1 0; 0 1; 0 -1]; b2 = [-5; 1; -1];
P2 = polytope(A2, b2);
res(end+1,1) = ~isIntersecting(P1, P2);

% 2D, degenerate contained in unbounded
A1 = [1 0; -1 0; 0 1]; b1 = [1; 1; 1];
P1 = polytope(A1,b1);
A2 = [1 0; -1 0; 0 1; 0 -1]; b2 = [0.5; 0.5; 0.5; -0.5];
P2 = polytope(A2,b2);
res(end+1,1) = isIntersecting(P1,P2);


% 2D, intersections of polytope and constrained zonotope
Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
A = [1 1 1]; b = 1;
cZ = conZonotope(Z,A,b);
% unbounded polytope, no intersection
A = [1 0]; b = -3;
P = polytope(A,b);
res(end+1,1) = ~isIntersecting(P,cZ);
% unbounded polytope, intersection
A = [1 0]; b = 0;
P = polytope(A,b);
res(end+1,1) = isIntersecting(P,cZ);
% unbounded polytope, intersection in a single point
A = [0 1]; b = -1.5;
P = polytope(A,b);
res(end+1,1) = isIntersecting(P,cZ);
% bounded polytope, no intersection
A = [-1 -1; 1 0; 0 1]; b = [-4; 5; 5];
P = polytope(A,b);
res(end+1,1) = ~isIntersecting(P,cZ);
% bounded polytope, intersection
A = [-1 -1; 1 0; 0 1]; b = [-2; 5; 5];
P = polytope(A,b);
res(end+1,1) = isIntersecting(P,cZ);
% bounded polytope, intersection in a single point
A = [-1 -1; 1 0; 0 1]; b = [-3; 5; 5];
P = polytope(A,b);
res(end+1,1) = isIntersecting(P,cZ);
% intersection with empty polytope
A = [1 0; -1 0]; b = [-1; -1];
P = polytope(A,b);
res(end+1,1) = ~isIntersecting(P,cZ);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
