function res = test_polyZonotope_isequal
% test_polyZonotope_isequal - unit test function of isequal
%
% Syntax:
%    res = test_polyZonotope_isequal
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

% Authors:       Mark Wetzlinger
% Written:       01-May-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% create polyZonotopes
c = [2; -3];
G1 = [2, 3;
     -1, 0];
GI1 = [1, 2, 4;
       5, 6, 0];
E1 = [2, 0;
      0, 1];
pZ1 = polyZonotope(c,G1,GI1,E1);
G2 = G1 + ones(2,2);
pZ2 = polyZonotope(c,G2,GI1,E1);
GI3 = [1, 2, 0, 4;
       5, 6, 0, 0];
pZ3 = polyZonotope(c,G1,GI3,E1);

% check result
res(end+1,1) = isequal(pZ1,pZ3);
res(end+1,1) = ~isequal(pZ1,pZ2);

% test random
pZ_random = polyZonotope.generateRandom();
res(end+1,1) = isequal(pZ_random, pZ_random);

% test non equal polyZonotopes
res(end+1,1) = ~isequal(pZ1, 1+pZ1);
% ...but should work with increased tolerance
tol = 0.001;
res(end+1,1) = isequal(pZ1, tol/2+pZ1, tol);

% slightly shrink polyZonotope
res(end+1,1) = isequal(pZ1, (1-tol/2) * pZ1, tol);
% slightly enlarge polyZonotope
res(end+1,1) = isequal(pZ1, (1+tol/2) * pZ1, tol);

% different sign in independent generator matrix
pZ4 = polyZonotope(c,G1,-GI1,E1);
res(end+1,1) = isequal(pZ1,pZ4);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
