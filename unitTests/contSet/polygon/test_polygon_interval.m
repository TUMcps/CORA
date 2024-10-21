function res = test_polygon_interval()
% test_polygon_interval - unit test function for polygon/interval
%
% Syntax:
%    res = test_polygon_interval()
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
% Written:       12-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% test empty
pgon = polygon.empty(2);
I = interval(pgon);
assert(representsa(I,'emptySet'))

% test scalar
V = [1;2];
pgon = polygon(V);
I = interval(pgon);
assert(contains(I,V));

% generate data
x = gallery('uniformdata',30,1,1);
y = gallery('uniformdata',30,1,10);
ind = boundary(x,y);
x = x(ind);
y = y(ind);
V = [x y]';

% get polygon
pgon = polygon(V);
I = interval(pgon);
assert(all(contains(I,V)));

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
