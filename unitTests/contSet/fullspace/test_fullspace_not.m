function res = test_fullspace_not
% test_fullspace_not - unit test function of '~'-operator
%
% Syntax:
%    res = test_fullspace_not
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
% Written:       07-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init fullspace
n = 2;
fs = fullspace(n);

% compute complement
O = ~fs;

% true solution
O_ = emptySet(n);

% check solution
res = O == O_;

% different syntax
O = not(fs);
res(end+1,1) = O == O_;

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
