function res = test_emptySet_isempty
% test_emptySet_isempty - unit test function of isempty
%
% Syntax:  
%    res = test_emptySet_isempty
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
% Written:      05-April-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% init empty set
n = 2;
O = emptySet(n);

% check result
res = isempty(O);

%------------- END OF CODE --------------