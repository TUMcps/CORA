function res = test_emptySet_and
% test_emptySet_and - unit test function of and
%
% Syntax:
%    res = test_emptySet_and
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

% intersection with itself
O_ = O & O;
res = isequal(O,O_);

% init zonotope
Z = zonotope(zeros(n,1),eye(n));

% intersection with zonotope
O_ = O & Z;
res(end+1,1) = isequal(O,O_);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
