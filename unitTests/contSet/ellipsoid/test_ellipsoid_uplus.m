function res = test_ellipsoid_uplus
% test_ellipsoid_uplus - unit test function of uplus
%
% Syntax:
%    res = test_ellipsoid_uplus
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
Q = [2.7 -0.2;-0.2 2.4];
q = [1;2];
E = ellipsoid(Q, q);

pE = +E;

% compare with E
resvec(end+1) = isequal(pE, E);

% test empty case
resvec(end+1) = isemptyobject(+ellipsoid.empty(2));

% add results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
