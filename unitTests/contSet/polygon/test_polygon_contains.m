function res = test_polygon_contains()
% test_polygon_contains - unit test function for polygon/contains
%
% Syntax:
%    res = test_polygon_contains
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
% See also: polygon

% Authors:       Tobias Ladner
% Written:       11-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init polygon via vertices
V = [4 0; 3 2; 0 3; -3 0; -2 -3; 0 -3; 2 -2]';
pgon = polygon(V);

% check if vertices are contained
assert(all(contains(pgon,V)))

% check if center is contained
assert(all(contains(pgon,center(pgon))))

% check containment with itself
assert(contains(pgon,pgon))

% check containment with slightly enlarged polygon
pgon2 = enlarge(pgon,1.1);
assert(contains(pgon2,pgon))
assert(~contains(pgon,pgon2))

% check containment of zonotope
Z = zonotope([0;0],eye(2)*0.5);
assert(contains(pgon,Z))

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
