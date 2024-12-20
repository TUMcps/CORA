function res = test_polygon_uplus()
% test_polygon_uplus - unit test function for polygon/uplus
%
% Syntax:
%    res = test_polygon_uplus()
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
% Written:       28-June-2023
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
assert(isequal(pgon, +pgon));

% check empty case
pgon = polygon();
assert(isequal(pgon, +pgon));

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
