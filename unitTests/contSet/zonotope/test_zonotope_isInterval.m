function res = test_zonotope_isInterval
% test_zonotope_isInterval - unit test function of isInterval
%
% Syntax:  
%    res = test_zonotope_isInterval
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
% Written:      09-August-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% create zonotopes
c1 = [0; 0];
G1 = [2 0; 0 1];
Z1 = zonotope(c1,G1);

c2 = [1; 0];
G2 = [2 1; -1 4];
Z2 = zonotope(c2,G2);

% check result
res = isInterval(Z1) && ~isInterval(Z2);

%------------- END OF CODE --------------
