function res = test_GJKalgorithm()
% test_GJKalgorithm - unit test function for the GJK algorithm; only 2D
%    tests with corner cases for different set representations
%
% Syntax:
%    res = test_GJKalgorithm()
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
% Written:       24-April-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% zonotopes
Z1 = zonotope([1;-2],[1 0 1; -1 2 1]);
Z2 = zonotope([5; 3],[0 1; -1 1]); % no intersection
Z3 = zonotope([2; 4],[0 1; -1 1]); % vertex-edge intersection
Z4 = zonotope([3; 3],[0 1; -1 1]); % edge-edge intersection
Z5 = zonotope([1; 3],[0 1; -1 1]); % vertex-vertex intersection
Z6 = zonotope([1; 2],[0 1; -1 1]); % full-dimensional intersection

% check for intersection with GJK algorithm
res(end+1,1) = ~GJKalgorithm(Z1,Z2);
res(end+1,1) = GJKalgorithm(Z1,Z3);
res(end+1,1) = GJKalgorithm(Z1,Z4);
res(end+1,1) = GJKalgorithm(Z1,Z5);
res(end+1,1) = GJKalgorithm(Z1,Z6);


% intervals
I1 = interval([-3;1],[2;4]);
I2 = interval([0;5],[3;6]); % no intersection
I3 = interval([2;3],[3;5]); % edge-edge intersection
I4 = interval([2;4],[3;6]); % vertex-vertex intersection
I5 = interval([0;3],[1;5]); % full-dimensional intersection

res(end+1,1) = ~GJKalgorithm(I1,I2);
res(end+1,1) = GJKalgorithm(I1,I3);
res(end+1,1) = GJKalgorithm(I1,I4);
res(end+1,1) = GJKalgorithm(I1,I5);


% polytopes
P1 = polytope([1 0; -1 1; -1 -1],[1;1;1]);
P2 = polytope([0 1.5; -1 2; -1.5 0.5]'); % no intersection
P3 = polytope([0 1; -1 2; -1.5 0.5]'); % vertex-edge intersection
P4 = polytope([-1 0; -1 2; -1.5 0.5]'); % vertex-vertex intersection
P5 = polytope([-0.5 0; -1 2; -1.5 0.5]'); % full-dimensional intersection
P6 = polytope([0 1; -1 0; -1 2; -1.5 0.5]'); % edge-edge intersection

res(end+1,1) = ~GJKalgorithm(P1,P2);
res(end+1,1) = GJKalgorithm(P1,P3);
res(end+1,1) = GJKalgorithm(P1,P4);
res(end+1,1) = GJKalgorithm(P1,P5);
res(end+1,1) = GJKalgorithm(P1,P6);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
