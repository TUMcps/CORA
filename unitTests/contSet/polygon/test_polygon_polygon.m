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

% test constructor

pgon = polygon();
pgon = polygon([1 2 3],[2 7 3]);
pgon = polygon([1;2;3],[2;7;3]);
pgon = polygon([1 2 3; 2 7 3]);

pgon = polygon(pgon);
pgon = polygon(pgon.set);

% check if points are there
V = [1 7 5 9; 2 3 9 0];

pgon = polygon(V);
assert(isequal(V,pgon.set.Vertices'));

pgon = polygon(V(1,:),V(2,:));
assert(isequal(V,vertices(pgon)));

% check degenerate polygons
V = [1;1];
pgon = polygon(V);
assert(contains(pgon,V))

V = [1:10;2:11];
pgon = polygon(V);
assert(all(contains(pgon,V)))

% unable to initialize infinite vertices
assertThrowsAs(@polygon, 'CORA:wrongInputInConstructor',Inf(2,1));

% nan is ok due to multiple regions
V1 = [0 0 1 1; 1 0 0 1];
V2 = V1 + 4;
pgon = polygon([V1 nan(2,1) V2]);

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
