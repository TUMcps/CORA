function res = test_polygon_subtract()
% test_polygon_subtract - unit test function for polygon/subtract
%
% Syntax:
%    res = test_polygon_subtract()
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

pgon = subtract(polygon(I1),polygon(I2*0.1+1));
assert(pgon.nrOfHoles == 1)
holes = getHoles(pgon);
assert(isequal(holes{1},polygon(I2*0.1+1)))

pgon = subtract(subtract(polygon(I1),polygon(I2*0.1+1)),polygon(I2*0.1-1));
assert(pgon.nrOfHoles == 2)
holes = getHoles(pgon);
assert(isequal(holes{1},polygon(I2*0.1+1)))
assert(isequal(holes{2},polygon(I2*0.1-1)))

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
