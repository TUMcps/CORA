function res = test_polygon_isBounded()
% test_polygon_isBounded - unit test function for polygon/isBounded
%
% Syntax:
%    res = test_polygon_isBounded()
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
% Written:       16-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty polygon
pgon = polygon.empty(2);
assert(isBounded(pgon));

% get normal polygon
V = [1 0 0 1; 0 0 1 1];
pgon = polygon(V);
assert(isBounded(pgon))

% test completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
