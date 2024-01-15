function res = test_zonoBundle_uplus
% test_zonoBundle_uplus - unit test function of uplus
%
% Syntax:
%    res = test_zonoBundle_uplus
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

% plus
pzB = +zB;

% compare with zB
resvec(end+1) = isequal(pzB, zB);

% test empty case
resvec(end+1) = isemptyobject(+zonoBundle.empty(2));

% add results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
