function res = test_polygon_polygon()
% test_polygon_polygon - unit test function for the polygon constructor
%
% Syntax:
%    res = test_polygon_polygon()
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
% Written:       25-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

% test constructor

pgon = polygon();
pgon = polygon([1 2 3],[2 7 3]);
pgon = polygon([1;2;3],[2;7;3]);
pgon = polygon([1 2 3; 2 7 3]);
resvec(end+1) = true;

pgon = polygon(pgon);
pgon = polygon(pgon.set);
resvec(end+1) = true;

% check if points are there
V = [1 7 5 9; 2 3 9 0];

pgon = polygon(V);
resvec(end+1) = isequal(V,pgon.set.Vertices');

pgon = polygon(V(1,:),V(2,:));
resvec(end+1) = isequal(V,pgon.set.Vertices');

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
