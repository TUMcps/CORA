function res = test_conZonotope_isZero
% test_conZonotope_isZero - unit test function of isZero
%
% Syntax:  
%    res = test_conZonotope_isZero
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
% See also: none

% Author:       Mark Wetzlinger
% Written:      16-March-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% empty case
res(1) = ~isZero(conZonotope());

% check different cases...

% true isZero cases
cZ = conZonotope(zeros(3,1));
res(end+1) = isZero(cZ);

cZ = conZonotope(zeros(2,1),zeros(2));
res(end+1) = isZero(cZ);

% shifted center
cZ = conZonotope(ones(2,1),zeros(2,1));
res(end+1) = ~isZero(cZ);

% zero-centered, but with generator
cZ = conZonotope(zeros(3,1),[1; 0; -1]);
res(end+1) = ~isZero(cZ);

% below tolerance
cZ = conZonotope([0.5;0.5],[0.02 -0.1; 0.2 0.05]);
tol = 1;
res(end+1) = isZero(cZ,tol);

% combine tests
res = all(res);

%------------- END OF CODE --------------