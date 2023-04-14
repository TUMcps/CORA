function res = test_zonotope_deleteAligned
% test_zonotope_deleteAligned - unit test function of deleteAligned
%
% Syntax:  
%    res = test_zonotope_deleteAligned
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

% Author:       Mark Wetzlinger
% Written:      26-August-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% create a zonotope
Z_cent = zeros(2,1);
% aligned generators differ by scaling, sign
Z_gens = [4 2 2 3 1 -4;
          2 3 1 0 2 -2];
Z = zonotope(Z_cent, Z_gens);

% delete aligned generators
Z_del = deleteAligned(Z);

% true matrix
Z_gens_true = [10 2 3 1;
               5  3 0 2];

% check results
res = compareMatrices(Z_gens_true,generators(Z_del));

%------------- END OF CODE --------------