function res = test_polytope_minkDiff
% test_polytope_minkDiff - unit test function of minkDiff
%
% Syntax:
%    res = test_polytope_minkDiff
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
% Written:       07-September-2022
% Last update:   01-December-2022
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% 1D, bounded - bounded
% {1 <= x <= 1} - {-0.7 <= x <= 0.2} = {-0.3 <= x <= 0.8}
P1 = polytope([1; -1],[1; 1]);
P2 = polytope([1; -1],[0.2; 0.7]);
P = minkDiff(P1,P2);
P_true = polytope([1;-1],[0.8; 0.3]);
res(end+1,1) = isequal(P,P_true,1e-10);


% 2D, convert zonotope to polytopes
Z_m = zonotope([1; 1],[1 0 1; 0 1 1]);
P_m = polytope(Z_m);
Z_s = zonotope([-0.5; 1],[0.5 0; 0 0.5]);
P_s = polytope(Z_s);
% compute Minkowski difference
P_minkDiff = minkDiff(P_m,P_s);
% true result in zonotope representation and converted to a polytope
Z_true = zonotope([1.5; 0],[0.5 0 1; 0 0.5 1]);
P_true = polytope(Z_true);
% compare results
res(end+1,1) = P_minkDiff == P_true;
% visualization
% figure; hold on;
% plot(Z_m);
% plot(Z_s,[1,2],'k');
% plot(Z_true,[1,2],'Color',colorblind('r'));
% plot(Z_s + Z_true,[1,2],'g--');

% 2D, bounded - unbounded
P1 = polytope([1 0; -1 0; 0 1; 0 -1],[1;1;1;1]);
P2 = polytope([1 0; -1 0; 0 1],[1;1;1]);
P_minkDiff = minkDiff(P1,P2);
res(end+1,1) = representsa(P_minkDiff, 'emptySet');

% 2D, degenerate - unbounded
P1 = polytope([1 0; -1 0; 0 1; 0 -1],[1;1;0;0]);
P2 = polytope([1 0; -1 0; 0 1],[1;1;1]);
P_minkDiff = minkDiff(P1,P2);
res(end+1,1) = representsa(P_minkDiff, 'emptySet');

% 2D, unbounded - degenerate
P1 = polytope([1 0; -1 0; 0 1],[1;1;1]);
P2 = polytope([1 0; -1 0; 0 1; 0 -1],[1;1;0;0]);
P_minkDiff = minkDiff(P1,P2);
P_true = polytope([1 0; -1 0; 0 1],[0;0;1]);
res(end+1,1) = P_minkDiff == P_true;

% 2D, unbounded - unbounded
P1 = polytope([1 0; -1 0; 0 -1],[1;1;1]);
P2 = polytope([1 0; -1 0; 0 1],[1;1;1]);
P_minkDiff = minkDiff(P1,P2);
res(end+1,1) = representsa(P_minkDiff, 'emptySet');

% 2D, two identical polytopes, unbounded
P1 = polytope([1 0; -1 0; 0 1],[1;1;1]);
P_minkDiff = minkDiff(P1,P1);
P_true = polytope([1 0; -1 0; 0 1],[0;0;0]);
res(end+1,1) = P_minkDiff == P_true;

% 2D, two identical polytopes, bounded
P1 = polytope([-1 0; 0 -1],[0;0]);
P_minkDiff = minkDiff(P1,P1);
res(end+1,1) = P_minkDiff == P1;

% 2D, degenerate - degenerate (equality constraints)
P1 = polytope([1 0 0; -1 0 0; 0 1 0; 0 -1 0],[1;1;1;1],[0 0 1],0);
P2 = polytope([1 0 0; -1 0 0; 0 1 0; 0 -1 0],0.5*[1;1;1;1],[0 0 1],0);
P_minkDiff = minkDiff(P1,P2);
P_true = polytope([1 0 0; -1 0 0; 0 1 0; 0 -1 0],0.5*[1;1;1;1],[0 0 1],0);
res(end+1,1) = P_minkDiff == P_true;

% 2D, degenerate - vector (only inequality constraints)
% ...line from (-1,0) to (1,0)
P1 = polytope([1 0; -1 0],[1;1],[0 1],0);
P2 = [-1;1];
P_minkDiff = minkDiff(P1,P2);
P_true = polytope([1 0; -1 0],[2;0],[0 1],-1);
% ...line from (0,-1) to (2,-1)
res(end+1,1) = P_minkDiff == P_true;


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
