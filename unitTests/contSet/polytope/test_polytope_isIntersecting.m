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
P1 = polytope(1,1);
P2 = polytope(-1,5);
res(end+1,1) = isIntersecting(P1,P2);

% 1D, unbounded & degenerate
P1 = polytope(1,1);
P2 = polytope([],[],1,-5);
res(end+1,1) = isIntersecting(P1,P2);

% 1D, bounded & unbounded, intersection in a single point
P1 = polytope([1;-1],[1;1]);
P2 = polytope([],[],-1,1);
res(end+1,1) = isIntersecting(P1,P2);

% 1D, unbounded & point
P1 = polytope(1,1);
P2 = 5;
res(end+1,1) = ~isIntersecting(P1,P2);


% 2D, intersection with fully empty set
P1 = polytope(zeros(0,2),[]);
P2 = polytope([1 0],1);
res(end+1,1) = ~isIntersecting(P1,P2);

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
P = polytope([1 0],-3);
res(end+1,1) = ~isIntersecting(P,cZ);
% unbounded polytope, intersection
P = polytope([1 0],0);
res(end+1,1) = isIntersecting(P,cZ);
% unbounded polytope, intersection in a single point
P = polytope([0 1],-1.5);
res(end+1,1) = isIntersecting(P,cZ);
% bounded polytope, no intersection
P = polytope([-1 -1; 1 0; 0 1],[-4; 5; 5]);
res(end+1,1) = ~isIntersecting(P,cZ);
% bounded polytope, intersection
P = polytope([-1 -1; 1 0; 0 1],[-2; 5; 5]);
res(end+1,1) = isIntersecting(P,cZ);
% bounded polytope, intersection in a single point
P = polytope([-1 -1; 1 0; 0 1],[-3; 5; 5]);
res(end+1,1) = isIntersecting(P,cZ);
% intersection with empty polytope
P = polytope([1 0; -1 0],[-1; -1]);
res(end+1,1) = ~isIntersecting(P,cZ);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
