function res = test_fullspace_dim
% test_fullspace_dim - unit test function of dim
%
% Syntax:
%    res = test_fullspace_dim
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
% Written:       05-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% n-dimensional fullspace
n = 2;
fs = fullspace(n);
assert(dim(fs) == n);

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
