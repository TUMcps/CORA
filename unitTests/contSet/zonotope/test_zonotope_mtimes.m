function res = test_zonotope_mtimes
% test_zonotope_mtimes - unit test function of mtimes
%
% Syntax:  
%    res = test_zonotope_mtimes
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

% Author:       Matthias Althoff
% Written:      26-July-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% instantiate zonotopes
Z_empty = zonotope();
Z = zonotope([-4, -3, -2, -1; 1, 2, 3, 4]);

% create parallelotopes
M = [-1 2; 3 -4];

% compute linear maps
Z_empty_ = M*Z_empty;
Z_ = M*Z;

% obtain center and generator matrix
c_ = center(Z_);
G_ = generators(Z_);

% true result
true_c = [6; -16];
true_G = [7, 8, 9; ...
         -17, -18, -19];

% check result
res = isempty(Z_empty_) && compareMatrices(c_,true_c) ...
    && compareMatrices(G_,true_G);

%------------- END OF CODE --------------
