function res = testLong_polyZonotope_polygon
% testLong_polyZonotope_polygon - unit test function for the
%    conversion of a polynomial zonotope to as polygon
%
% Syntax:
%    res = testLong_polyZonotope_polygon
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
% See also: -

% Authors:       Tobias Ladner
% Written:       11-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = true;

% empty case
pgon = polygon(polyZonotope.empty(1));
resvec(end+1) = representsa_(pgon,'emptySet',eps);

% point case
pgon = polygon(polyZonotope([1; 2]));
resvec(end+1) = ~representsa_(pgon,'emptySet',eps);

% line case
pgon = polygon(polyZonotope([1; 2], [1; 1]));
resvec(end+1) = ~representsa_(pgon,'emptySet',eps);
pgon = polygon(polyZonotope([1; 2], [2, 4; 1, 2]));
resvec(end+1) = ~representsa_(pgon,'emptySet',eps);

% complex cases
pgon = polygon(polyZonotope([1; 2], [1, -1; 1, 3]));
resvec(end+1) = compareMatrices( ...
    pgon.set.Vertices, ...
    [1.000, -2.000; -1.000, 4.000; 1.000, 6.000; 3.000, 0.000], ...
    1e-8, "equal", true);

pgon = polygon(polyZonotope([1; 2], [1, 0; 0, 1; 1, 1]', [], [1, 0, 1; 0, 1, 1]));
resvec(end+1) = compareMatrices( ...
    [0.7421875000000000, 1.7421875000000000; ...
    -0.0078125000000000, 1.9921875000000000; ...
    -0.6328125000000000, 1.9921875000000000; ...
    -0.4921875000000000, 2.2578125000000000; ...
    0.1250000000000000, 2.5625000000000000; ...
    0.7500000000000000, 2.8750000000000000; ...
    1.5078125000000000, 3.2578125000000000; ...
    2.7578125000000000, 3.5078125000000000; ...
    1.8828125000000000, 1.7578125000000000; ...
    1.5000000000000000, 1.0000000000000000; ...
    1.1250000000000000, 0.2500000000000000; ...
    0.9993489583333000, 0.3847656250000000; ...
    0.9980468750000000, 0.9042968750000000; ...
    0.9296875000000000, 1.4296875000000000; ...
    0.7773437500000000, 1.7148437500000000; ...
    ]', ...
    pgon.set.Vertices', ...
    1e-6, "subset", true);

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
