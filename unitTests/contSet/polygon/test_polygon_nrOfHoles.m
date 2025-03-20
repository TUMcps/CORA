function res = test_polygon_nrOfHoles()
% test_polygon_nrOfHoles - unit test function for polygon/nrOfHoles
%
% Syntax:
%    res = test_polygon_nrOfHoles()
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
assert(pgon.nrOfHoles == 0)

pgon = subtract(polygon(I1),polygon(I2*0.1+1));
assert(pgon.nrOfHoles == 1)

pgon = subtract(subtract(polygon(I1),polygon(I2*0.1+1)),polygon(I2*0.1-1));
assert(pgon.nrOfHoles == 2)

pgon = polygon(I2-5) | polygon(I2+5);
assert(pgon.nrOfHoles == 0)

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
