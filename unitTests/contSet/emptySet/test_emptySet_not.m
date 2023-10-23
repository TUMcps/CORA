function res = test_emptySet_not
% test_emptySet_not - unit test function of '~'-operator
%
% Syntax:
%    res = test_emptySet_not
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

% init empty set
n = 2;
O = emptySet(n);

% compute complement
fs = ~O;

% true solution
fs_ = fullspace(n);

% check solution
res = fs == fs_;

% different syntax
fs = not(O);
res(end+1,1) = fs == fs_;

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
