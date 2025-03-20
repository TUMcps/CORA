function res = test_polygon_getHoles()
% test_polygon_getHoles - unit test function for polygon/getHoles
%
% Syntax:
%    res = test_polygon_getHoles()
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

% init
I1 = interval([-3;-2],[2;3]);
I2 = interval([-1;-1],[1;1]);

% no hole
pgon = polygon(I1);
holes = getHoles(pgon);
assert(isempty(holes))

pgon = polygon(I2-5) | polygon(I2+5);
holes = getHoles(pgon);
assert(isempty(holes))

% single hole
pgon = subtract(polygon(I1),polygon(I2*0.1+1));
holes = getHoles(pgon);
assert(isequal(holes{1},polygon(I2*0.1+1)))

% two holes
pgon = subtract(subtract(polygon(I1),polygon(I2*0.1+1)),polygon(I2*0.1-1));
holes = getHoles(pgon);
assert(isequal(holes{1},polygon(I2*0.1+1)))
assert(isequal(holes{2},polygon(I2*0.1-1)))

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
