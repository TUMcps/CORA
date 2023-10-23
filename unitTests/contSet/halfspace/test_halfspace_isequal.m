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
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Mark Wetzlinger
% Written:       17-September-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% instantiate halfspaces
h1 = halfspace([2;3;-1],3);
h2 = halfspace([1;3;-1],3);
h3 = h1;

% combine tests
res = ~isequal(h1,h2) && isequal(h1,h3);

% ------------------------------ END OF CODE ------------------------------
