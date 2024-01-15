function res = test_zonotope_vertices
% test_zonotope_vertices - unit test function of vertices
%
% Syntax:
%    res = test_zonotope_vertices
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
% See also: none

% Authors:       Matthias Althoff, Tobias Ladner
% Written:       26-July-2016
% Last update:   09-February-2017
%                01-June-2023 (TL, more cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

% empty set
Z = zonotope.empty(2);
V = vertices(Z);
resvec(end+1) = isempty(V) && isnumeric(V) && all(size(V) == [2,0]);

% simple zonotope
Z = zonotope([0 1 0 1; 0 1 1 0]);
V = vertices(Z);
V0 = [2 0 -2 -2 0 2; 2 2 0 -2 -2 0];
resvec(end+1) = compareMatrices(V,V0);

% zonotope
Z = zonotope([-4, -3, -2, -1; 1, 2, 3, 4]);
V = vertices(Z);
V0 = [-8, 0, 2, -4, -4, -10; ...
              2, 0, -8, -4, 6, 10];
resvec(end+1) = compareMatrices(V,V0);

% 1d zonotope
Z = zonotope([2 1]);
V = vertices(Z);
V0 = [1 3];
resvec(end+1) = compareMatrices(V,V0);

% degenerate case
Z = zonotope([1;2;3],[1 1; 0 1; 1 1]);
V = vertices(Z);
V0 = [ ...
 1.000, 3.000, 1.000, -1.000 ; ...
 3.000, 3.000, 1.000, 1.000 ; ...
 3.000, 5.000, 3.000, 1.000 ; ...
];
resvec(end+1) = compareMatrices(V,V0);

% another degenerate case
Z = zonotope([1;2;3],[1 1 2 -1; 0 1 2 0; 1 1 2 -1]);
V = vertices(Z);
V0 = [ ...
 2.000, -4.000, 6.000, 0.000 ; ...
 5.000, -1.000, 5.000, -1.000 ; ...
 4.000, -2.000, 8.000, 2.000 ; ...
];
resvec(end+1) = compareMatrices(V,V0);

% case with zero dimension
Z = zonotope([1;2;0],[1 1; 0 1; 0 0]);
V = vertices(Z);
V0 = [ ...
 3.000, 1.000, 1.000, -1.000 ; ...
 3.000, 3.000, 1.000, 1.000 ; ...
 0.000, 0.000, 0.000, 0.000 ; ...
];
resvec(end+1) = compareMatrices(V,V0);

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
