function res = test_halfspace_plus
% test_halfspace_plus - unit test function of plus
%
% Syntax:  
%    res = test_halfspace_plus
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
% Written:      16-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% 1. empty case
h = halfspace();
v = rand(2,1);
res_empty = isempty(h + v);

% 2. dimension mismatch
res_dim = true;
h = halfspace(randn(2,1),1);
v = rand(3,1);
try
    h + v; % should throw error here
    res_dim = false;
end

% combine tests
res = res_empty && res_dim;

%------------- END OF CODE --------------