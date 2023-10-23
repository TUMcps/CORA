function res = test_polytope_distance
% test_polytope_distance - unit test function of distance
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
% See also: -

% Authors:       Viktor Kotsev
% Written:       07-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% 2D, bounded; bounded
P1 = polytope([ 1 0;-1 0; 0 1; 0 -1],[6;-5;1;1]);
P2 = polytope([ 1 0;-1 0; 0 1; 0 -1],[1;1;1;1]);
% compute distance and compare to true value
dist = distance(P1,P2);
dist_true = 4;
res(end+1,1) = withinTol(dist,dist_true);

% 2D, bounded, degenerate; bounded, degenerate
P1 = polytope([1 0; 0 1; -1 0; 0 -1],[1;1;1;-1]);
P2 = polytope([1 0; 0 1; -1 0; 0 -1],[6;3;-4;-3]);
% compute distance and compare to true value
dist = distance(P1,P2);
dist_true = 3.605551275;
res(end+1,1) = withinTol(dist,dist_true);


% Unbounded case
% U1 = polytope([1 0; 0 1;-1 0],[1;1;1]);
% U2 = polytope([1 0;0 1;-1 0],[6;3 ;-4]);
% dist_3 =(U1,U2);
% a3 = 3.605551275;
%res(3) = withinTol(dist_3, a3);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
