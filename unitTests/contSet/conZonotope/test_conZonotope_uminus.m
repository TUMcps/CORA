function res = test_conZonotope_uminus
% test_conZonotope_uminus - unit test function of uminus
%
% Syntax:  
%    res = test_conZonotope_uminus
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
Z = [0 2 0 2; 0 0 2 2];
A = [2 -4 2];
b = 0;
cZ = conZonotope(Z,A,b);

% negate
ncZ = -cZ;
resvec(end+1) = all(ncZ.Z == -Z, 'all');
resvec(end+1) = all(ncZ.A == A, 'all');
resvec(end+1) = all(ncZ.b == b, 'all');

% compare with -1 * cPZ
resvec(end+1) = isequal(ncZ, -1*cZ);

% test empty case
resvec(end+1) = isemptyobject(-conZonotope());

% add results
res = all(resvec);

%------------- END OF CODE --------------