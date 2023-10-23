function res = test_emptySet_isIntersecting
% test_emptySet_isIntersecting - unit test function of isIntersecting
%
% Syntax:
%    res = test_emptySet_isIntersecting
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
res = ~isIntersecting(O,O);

% init zonotope
Z = zonotope(zeros(2,1));

% intersection with zonotope
res(end+1,1) = ~isIntersecting(O,Z);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
