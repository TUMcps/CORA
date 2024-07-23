function res = test_polytope_matPolytope
% test_polytope_matPolytope - unit test function of matPolytope
%
% Syntax:
%    res = test_polytope_matPolytope
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
% Written:       12-July-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% 2D, instantiation via vertices
V = [1 0; 0 1; -1 1; -1 -1; 0 -1]';
P = polytope(V);
matP = matPolytope(P);
res(end+1,1) = all(all(V == reshape(matP.V,[2,5])));

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
