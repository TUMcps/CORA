function res = test_polygon_supportFunc()
% test_polygon_supportFunc - unit test function for polygon/supportFunc
%
% Syntax:
%    res = test_polygon_supportFunc
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
dir = [1 1];
[val,x] = supportFunc(pgon,dir);
assert(val == -Inf);
assert(isempty(x));

% init polygon via vertices
V_true = [4 0; 3 2; 0 3; -3 0; -2 -3; 0 -3; 2 -2]';
pgon = polygon(V_true);

[val,x] = supportFunc(pgon,dir,'lower');
assert(val == -5)
assert(compareMatrices(x,[-2;-3]))

[val,x] = supportFunc(pgon,dir,'upper');
assert(val == 5)
assert(compareMatrices(x,[3;2]))

[val,x] = supportFunc(pgon,dir,'range');
assert(isequal(val,interval(-5,5)));
assert(compareMatrices(x,[-2 3;-3 2]))

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
