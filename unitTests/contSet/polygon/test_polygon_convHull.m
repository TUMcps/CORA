function res = test_polygon_convHull()
% test_polygon_convHull - unit test function for polygon/convHull
%
% Syntax:
%    res = test_polygon_convHull
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
assert(isequal(pgon,convHull(pgon)));

% init polygon via vertices
V_conv = [4 0; 3 2; 0 3; -3 0; -2 -3; 0 -3; 2 -2]';
V = [4 0; 3 2; 0 3; -3 0; -2 -3; 0 -3; 1 -1; 2 -2]';
pgon = polygon(V);

% compute convHull
pgon_conv = convHull(pgon);
V = vertices(pgon_conv);
assert(compareMatrices(V_conv,V))

% init second polygon
% generate data
x = gallery('uniformdata',30,1,1);
y = gallery('uniformdata',30,1,10);
ind = boundary(x,y);
x = x(ind);
y = y(ind);
pgon2 = polygon(x,y)+[-2;1];

% compute convHull
pgon_conv = convHull(pgon,pgon2);

assert(isConvex(pgon_conv))
assert(contains(pgon_conv,pgon))
assert(contains(pgon_conv,pgon2,'exact',1e-6))

% init point polygon
p = [2;4];

% compute convHull
pgon_conv = convHull(pgon,p);

assert(isConvex(pgon_conv))
assert(contains(pgon_conv,pgon))
assert(contains(pgon_conv,p))

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
