function res = test_halfspace_dim
% test_halfspace_dim - unit test function of dim
%
% Syntax:
%    res = test_halfspace_dim
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
% Written:       27-September-2019
% Last update:   16-March-2021 (MW, add empty case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% empty case
n = 2;
hs = halfspace.empty(n);
res(end+1,1) = dim(hs) == n;

% fullspace case
n = 3;
hs = halfspace.Inf(n);
res(end+1,1) = dim(hs) == n;

% 2D halfspace
hs = halfspace([1;1],0.5);
res(end+1,1) = dim(hs) == 2;

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
