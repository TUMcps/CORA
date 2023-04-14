function res = test_zonotope_uplus
% test_zonotope_uplus - unit test function of uminus
%
% Syntax:  
%    res = test_zonotope_uplus
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

% Author:       Tobias Ladner
% Written:      06-April-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

resvec = true(0);

% init
cG = [0 2 0 2; 0 0 2 2];
Z = zonotope(cG);

% plus
pZ = +Z;
resvec(end+1) = all(pZ.Z == cG, 'all');

% compare with Z
resvec(end+1) = isequal(pZ, Z);

% test empty case
resvec(end+1) = isemptyobject(+zonotope());

% add results
res = all(resvec);

%------------- END OF CODE --------------