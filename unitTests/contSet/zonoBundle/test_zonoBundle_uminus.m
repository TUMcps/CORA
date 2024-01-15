function res = test_zonoBundle_uminus
% test_zonoBundle_uminus - unit test function of uminus
%
% Syntax:
%    res = test_zonoBundle_uminus
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

% Authors:       Tobias Ladner
% Written:       06-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = true(0);

% init
Z1 = zonotope([1 3 0; 1 0 2]);
Z2 = zonotope([0 2 2; 0 2 -2]);
zB = zonoBundle({Z1,Z2});

% negate
nzB = -zB;

% compare with -1 * zB
resvec(end+1) = isequal(nzB, -1*zB);

% test empty case
resvec(end+1) = isemptyobject(-zonoBundle.empty(2));

% add results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
