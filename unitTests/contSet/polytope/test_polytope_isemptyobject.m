function res = test_polytope_isemptyobject
% test_polytope_isemptyobject - unit test function of object emptiness
%    check
%
% Syntax:
%    res = test_polytope_isemptyobject
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
% Written:       03-June-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% 2D, only inequalities
A = [-1 0; 2 4; 1 -2]; b = [-1; 14; -1];
P = polytope(A,b);
res(end+1,1) = ~isemptyobject(P);

% 3D, only equalities, single point
Ae = [1 0 1; 0 1 -1; 1 0 -1]; be = [1;4;2];
P = polytope([],[],Ae,be);
res(end+1,1) = ~isemptyobject(P);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
