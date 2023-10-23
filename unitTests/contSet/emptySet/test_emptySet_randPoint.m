function res = test_emptySet_randPoint
% test_emptySet_randPoint - unit test function of randPoint
%
% Syntax:
%    res = test_emptySet_randPoint
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

% sample random point
p = randPoint(O);
res = all(size(p) == [n,0]);

% ------------------------------ END OF CODE ------------------------------
