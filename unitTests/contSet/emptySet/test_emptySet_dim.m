function res = test_emptySet_dim
% test_emptySet_dim - unit test function of dim
%
% Syntax:
%    res = test_emptySet_dim
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

% init empty set
n = 2;
O = emptySet(n);

% check dimension
res = dim(O) == n;

% ------------------------------ END OF CODE ------------------------------
