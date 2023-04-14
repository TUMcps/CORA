function res = test_zonotope_uminus
% test_zonotope_uminus - unit test function of uminus
%
% Syntax:  
%    res = test_zonotope_uminus
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

% negate
nZ = -Z;
resvec(end+1) = all(nZ.Z == -cG, 'all');

% compare with -1 * Z
resvec(end+1) = isequal(nZ, -1*Z);

% test empty case
resvec(end+1) = isemptyobject(-zonotope());

% add results
res = all(resvec);

%------------- END OF CODE --------------