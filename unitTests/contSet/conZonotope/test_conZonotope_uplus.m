function res = test_conZonotope_uplus
% test_conZonotope_uplus - unit test function of uplus
%
% Syntax:  
%    res = test_conZonotope_uplus
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

% plus
pcZ = +cZ;
resvec(end+1) = all(pcZ.Z == Z, 'all');
resvec(end+1) = all(pcZ.A == A, 'all');
resvec(end+1) = all(pcZ.b == b, 'all');

% compare with cPZ
resvec(end+1) = isequal(pcZ, cZ);

% test empty case
resvec(end+1) = isemptyobject(+conZonotope());

% add results
res = all(resvec);

%------------- END OF CODE --------------