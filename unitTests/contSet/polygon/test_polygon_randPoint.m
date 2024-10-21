function res = test_polygon_randPoint()
% test_polygon_randPoint - unit test function for polygon/randPoint
%
% Syntax:
%    res = test_polygon_randPoint()
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

tol = 1e-6;

% generate data
x = gallery('uniformdata',30,1,1);
y = gallery('uniformdata',30,1,10);
ind = boundary(x,y);
x = x(ind);
y = y(ind);

% get polygon
pgon = polygon(x,y);

% sample points
p = randPoint(pgon,100);
assert(all(contains(pgon,p,'exact',tol)));

% sample extreme points
p = randPoint(pgon,100,'extreme');
assert(all(contains(pgon,p,'exact',tol)));

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
