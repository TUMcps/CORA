function res = test_zonotope_isZero
% test_zonotope_isZero - unit test function of isZero
%
% Syntax:  
%    res = test_zonotope_isZero
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
% Written:      17-March-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% empty zonotope
Z = zonotope();
res(1) = ~isZero(Z);

% only origin
Z = zonotope(zeros(3,1));
res(end+1) = isZero(Z);

% shifted center
Z = zonotope(0.01*ones(4,1));
res(end+1) = ~isZero(Z);

% ...add tolerance
tol = 0.02;
res(end+1) = isZero(Z,tol);

% include generator matrix
Z = zonotope(ones(2,1),0.1*eye(2));
tol = 2;
res(end+1) = isZero(Z,tol);

% combine results
res = all(res);

%------------- END OF CODE --------------