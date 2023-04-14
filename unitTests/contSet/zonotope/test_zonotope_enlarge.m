function res = test_zonotope_enlarge
% test_zonotope_enlarge - unit test function of enlarge
%
% Syntax:  
%    res = test_zonotope_enlarge
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

% create zonotope
c = [-4; 1];
G = [-3, -2, -1; 2, 3, 4];
Z = zonotope(c,G);

% compute enlarged zonotope
Z_ = enlarge(Z,[2;1.5]);

% obtain center and generator matrix
c_ = center(Z_);
G_ = generators(Z_);

% true result
true_c = [-4; 1];
true_G = [-6, -4, -2; ...
            3, 4.5, 6];

% check result
res = compareMatrices(c_,true_c) && compareMatrices(G_,true_G);

%------------- END OF CODE --------------
