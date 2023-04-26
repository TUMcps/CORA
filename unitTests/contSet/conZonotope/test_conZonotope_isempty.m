function res = test_conZonotope_isempty
% test_conZonotope_isempty - unit test function of isempty
%
% Syntax:  
%    res = test_conZonotope_isempty
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

% Author:       Mark Wetzlinger
% Written:      21-April-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% check empty conZonotope object
cZ = conZonotope();
res = isempty(cZ);

% constrained zonotope
Z = [0 3 0 1;0 0 2 1];
A = [1 0 1]; b = 1;
cZ = conZonotope(Z,A,b);
res(end+1,1) = ~isempty(cZ);

% combine results
res = all(res);

%------------- END OF CODE --------------
