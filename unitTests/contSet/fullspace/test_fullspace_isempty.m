function res = test_fullspace_isempty
% test_fullspace_isempty - unit test function of isempty
%
% Syntax:  
%    res = test_fullspace_isempty
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

% init fullspace
n = 2;
fs = fullspace(n);

% check emptiness
res = ~isempty(fs);

%------------- END OF CODE --------------