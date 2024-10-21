function res = test_polygon_vertices()
% test_polygon_vertices - unit test function for polygon/vertices
%
% Syntax:
%    res = test_polygon_vertices
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

% check empty
pgon = polygon.empty(2);
V = vertices(pgon);
assert(compareMatrices(V,zeros(2,0)))

% init polygon via vertices
V_true = [4 0; 3 2; 0 3; -3 0; -2 -3; 0 -3; 2 -2]';
pgon = polygon(V_true);

V = vertices(pgon);
assert(compareMatrices(V_true,V))

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
