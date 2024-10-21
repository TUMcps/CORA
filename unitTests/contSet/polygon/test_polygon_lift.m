function res = test_polygon_lift()
% test_polygon_lift - unit test function for polygon/lift
%
% Syntax:
%    res = test_polygon_lift()
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

% generate data
x = gallery('uniformdata',30,1,1);
y = gallery('uniformdata',30,1,10);
ind = boundary(x,y);
x = x(ind);
y = y(ind);

% get polygon
pgon = polygon(x,y);
pgon_lift = lift(pgon,2,1:2);
assert(isequal(pgon,pgon_lift));

assertThrowsAs(@lift,'CORA:notDefined',pgon,3)
assertThrowsAs(@lift,'CORA:notDefined',pgon,4)
assertThrowsAs(@lift,'CORA:wrongValue',pgon,1)

% check empty case
pgon = polygon();
pgon_lift = lift(pgon,2,1:2);
assert(isequal(pgon,pgon_lift));

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
