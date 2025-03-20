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

% 2D, no constraints
A = zeros(0,2); b = zeros(0,0);
P = polytope(A,b);
assert(isemptyobject(P));

% 2D, only inequalities
A = [-1 0; 2 4; 1 -2]; b = [-1; 14; -1];
P = polytope(A,b);
assert(~isemptyobject(P));

% 3D, only equalities, single point
Ae = [1 0 1; 0 1 -1; 1 0 -1]; be = [1;4;2];
P = polytope([],[],Ae,be);
assert(~isemptyobject(P));

% 3D, no vertices
V = zeros(3,0);
P = polytope(V);
assert(isemptyobject(P));

% check special cases
P = polytope.Inf(2);
assert(~isemptyobject(P))

P = polytope.empty(2);
assert(~isemptyobject(P))

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
