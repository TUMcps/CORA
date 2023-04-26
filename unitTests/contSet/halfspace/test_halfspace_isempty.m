function res = test_halfspace_isempty
% test_halfspace_isempty - unit test function of isempty
%
% Syntax:  
%    res = test_halfspace_isempty
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
% Last update:  03-May-2020 (add empty case)
% Last revision:---

%------------- BEGIN CODE --------------

% instantiate halfspaces
hempty = halfspace();
h = halfspace([2;3;-1],3);

% combine tests
res = isempty(hempty) && ~isempty(h);

%------------- END OF CODE --------------