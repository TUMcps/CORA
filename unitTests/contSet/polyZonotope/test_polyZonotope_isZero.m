function res = test_polyZonotope_isZero
% test_polyZonotope_isZero - unit test function of isZero
%
% Syntax:  
%    res = test_polyZonotope_isZero
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

% empty polyZonotope
pZ = polyZonotope();
res(1) = ~isZero(pZ);

% only origin
pZ = polyZonotope(zeros(3,1));
res(end+1) = isZero(pZ);

% shifted center
pZ = polyZonotope(0.01*ones(4,1));
res(end+1) = ~isZero(pZ);

% ...add tolerance
tol = 0.02;
res(end+1) = isZero(pZ,tol);

% include dependent generator matrix
pZ = polyZonotope(ones(2,1),0.1*eye(2));
tol = 2;
res(end+1) = isZero(pZ,tol);

% dependent and independent generators
pZ = polyZonotope(zeros(2,1),[0.01 -0.02; 0.03 0.01],[0.05; -0.02],[2 0; 0 1]);
tol = 0.1;
res(end+1) = isZero(pZ,tol);


% combine results
res = all(res);

%------------- END OF CODE --------------