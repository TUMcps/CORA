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

res = true;

% n-dimensional fullspace
n = 2;
fs = fullspace(n);
res(end+1,1) = dim(fs) == n;

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
