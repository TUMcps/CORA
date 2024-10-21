function res = test_polygon_uminus()
% test_polygon_uminus - unit test function for polygon/uminus
%
% Syntax:
%    res = test_polygon_uminus()
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
assert(isequal(-1 * pgon, -pgon));

% check empty case
pgon = polygon();
assert(isequal(-1 * pgon, -pgon));

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
