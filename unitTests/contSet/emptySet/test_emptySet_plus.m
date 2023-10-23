function res = test_emptySet_plus
% test_emptySet_plus - unit test function of plus
%
% Syntax:
%    res = test_emptySet_plus
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

% addition with empty vector
p = double.empty(n,0);
O_ = O + p;
res = isequal(O,O_);
% different order
O_ = p + O;
res(end+1,1) = isequal(O,O_);

% addition with another empty set
O2 = emptySet(n);
O_ = O + O2;
res(end+1,1) = isequal(O,O_);

% init zonotope
Z = zonotope(zeros(n,1),eye(n));
O_ = O + Z;
res(end+1,1) = isequal(O,O_);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
