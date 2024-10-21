function res = test_polygon_projectHighDim()
% test_polygon_projectHighDim - unit test function for polygon/projectHighDim
%
% Syntax:
%    res = test_polygon_projectHighDim()
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
pgon_proj = projectHighDim(pgon,2,1:2);
assert(isequal(pgon,pgon_proj));

assertThrowsAs(@projectHighDim,'CORA:wrongValue',pgon,3)
assertThrowsAs(@projectHighDim,'CORA:wrongValue',pgon,4)
assertThrowsAs(@projectHighDim,'CORA:wrongValue',pgon,1)

% check empty case
pgon = polygon();
pgon_proj = projectHighDim(pgon,2,1:2);
assert(isequal(pgon,pgon_proj));

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
