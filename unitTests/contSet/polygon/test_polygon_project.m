function res = test_polygon_project()
% test_polygon_project - unit test function for polygon/project
%
% Syntax:
%    res = test_polygon_project()
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
pgon_proj = project(pgon,1:2);
assert(isequal(pgon,pgon_proj));

assert(isequal(project(pgon,1), interval(0.0222097785726014, 0.9538774735922314)))
assertThrowsAs(@project,'CORA:wrongValue',pgon,3)
assertThrowsAs(@project,'CORA:wrongValue',pgon,2:3)

% check empty case
pgon = polygon();
pgon_proj = project(pgon,1:2);
assert(isequal(pgon,pgon_proj));

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
