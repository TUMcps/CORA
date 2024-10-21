function res = test_zonotope_compact
% test_zonotope_compact - unit test function of redundancy removal in
%    zonotope representation
%    this encompasses checking the function nonzeroFilter
%
% Syntax:
%    res = test_zonotope_compact
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

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       26-July-2016
% Last update:   09-August-2020 (MW, enhance randomness)
% Last revision: 29-July-2023 (MW, rename '...compact')

% ------------------------------ BEGIN CODE -------------------------------

% 1D, zero removal/all
Z = zonotope(4, [2 -3 1 0 2]);
Z_compact = compact(Z,'zeros');
G_true = [2 -3 1 2];
assert(compareMatrices(Z_compact.G,G_true));
Z_compact = compact(Z,'all');
G_true = 8;
assert(compareMatrices(Z_compact.G,G_true));

% 2D, zero removal
Z = zonotope([1;5],[2 0 4; 6 0 0]);
Z_compact = compact(Z,'zeros');
G_true = [2, 4; 6, 0];
assert(compareMatrices(Z_compact.G,G_true));

% 2D, aligned generators differ by scaling, sign
Z = zonotope(zeros(2,1), [4 2 2 3 1 -4; 2 3 1 0 2 -2]);
Z_compact = compact(Z,'all');
G_true = [10 2 3 1; 5 3 0 2];
assert(compareMatrices(Z_compact.G,G_true));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
