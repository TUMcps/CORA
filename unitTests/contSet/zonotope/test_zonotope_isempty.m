function res = test_zonotope_isempty
% test_zonotope_isempty - unit test function of isempty
%
% Syntax:  
%    res = test_zonotope_isempty
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
% Written:      17-September-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% create zonotope
Z1 = zonotope([1, 2, 4;
               5, 6, 0;
              -1, 4, 8]);
Z2 = zonotope();

% check result
res = ~isempty(Z1) && isempty(Z2);

%------------- END OF CODE --------------
