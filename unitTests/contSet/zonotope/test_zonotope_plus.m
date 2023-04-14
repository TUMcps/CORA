function res = test_zonotope_plus
% test_zonotope_plus - unit test function of plus
%
% Syntax:  
%    res = test_zonotope_plus
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

% Author:       Matthias Althoff, Mark Wetzlinger
% Written:      26-July-2016
% Last update:  09-August-2020
% Last revision:---

%------------- BEGIN CODE --------------

% create zonotopes
Z1 = zonotope([-4, -3, -2, -1; 1, 2, 3, 4]);
Z2 = zonotope([1 10; -1 -10]);

% obtain results
Z_ = Z1+Z2;

% obtain center and generator matrix
c_ = center(Z_);
G_ = generators(Z_);

% true result
true_c = [-3; 0];
true_G = [-3, -2, -1, 10; ...
            2, 3, 4, -10];

% check result
res_val = compareMatrices(c_,true_c) && compareMatrices(G_,true_G);

% empty set
res_e = isempty(Z1+zonotope());

% add results
res = res_val && res_e;

%------------- END OF CODE --------------
