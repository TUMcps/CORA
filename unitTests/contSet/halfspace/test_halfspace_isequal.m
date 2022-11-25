function res = test_halfspace_isequal
% test_halfspace_isequal - unit test function of isequal
%
% Syntax:  
%    res = test_halfspace_isequal
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Mark Wetzlinger
% Written:      17-Sep-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% instantiate halfspaces
h1 = halfspace([2;3;-1],3);
h2 = halfspace([1;3;-1],3);
h3 = h1;

% combine tests
res = ~isequal(h1,h2) && isequal(h1,h3);

if res
    disp('test_isequal successful');
else
    disp('test_isequal failed');
end

%------------- END OF CODE --------------