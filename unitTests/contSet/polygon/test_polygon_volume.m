function res = test_polygon_volume()
% test_polygon_volume - unit test function for polygon/volume
%
% Syntax:
%    res = test_polygon_volume()
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
assert(withinTol(volume(pgon),0.685901695925958));

% check empty case
pgon = polygon();
assert(withinTol(volume(pgon),0));

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
