function res = test_polygon_getRegions()
% test_polygon_getRegions - unit test function for polygon/getRegions
%
% Syntax:
%    res = test_polygon_getRegions()
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
% Written:       13-February-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

I1 = interval([-3;-2],[2;3]);
I2 = interval([-1;-1],[1;1]);

pgon = polygon(I1);
regions = getRegions(pgon);
assert(isequal(regions{1},pgon))

pgon = polygon(I2-5) | polygon(I2+5);
regions = getRegions(pgon);
assert(isequal(regions{1},polygon(I2-5)))
assert(isequal(regions{2},polygon(I2+5)))

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
