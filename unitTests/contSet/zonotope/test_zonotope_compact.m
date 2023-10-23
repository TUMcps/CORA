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
% See also: -

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       26-July-2016
% Last update:   09-August-2020 (MW, enhance randomness)
% Last revision: 29-July-2023 (MW, rename '...compact')

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% create zonotope
c = [1;5];
G = [2,0,4; 6 0 0];
Z = zonotope(c,G);

% obtain zonotope without zeros
Z_ = compact(Z,'zeros');
G_ = Z_.G;

% true result
true_mat = [2, 4; 6, 0];

% check result
res(end+1,1) = compareMatrices(G_,true_mat);


% create a zonotope
Z_cent = zeros(2,1);
% aligned generators differ by scaling, sign
Z_gens = [4 2 2 3 1 -4;
          2 3 1 0 2 -2];
Z = zonotope(Z_cent, Z_gens);

% delete aligned generators
Z_del = compact(Z,'aligned');

% true matrix
Z_gens_true = [10 2 3 1;
               5  3 0 2];

% check results
res(end+1,1) = compareMatrices(Z_gens_true,Z_del.G);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
