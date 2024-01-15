function res = test_ellipsoid_ne
% test_ellipsoid_ne - unit test function of '~=' operator
%
% Syntax:
%    res = test_ellipsoid_ne
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

% Authors:       Mark Wetzlinger
% Written:       28-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty ellipsoid
E = ellipsoid.empty(2);

res = ~(E ~= E);

% 2D ellipsoids
E1 = ellipsoid([2 0; 0 1],[1;-1]);
E2 = ellipsoid([2 0; 0 1],[1;0]);
E3 = ellipsoid([2 0; 0 1],[1;-1],1e-10);

% compare
res(end+1,1) = ~(E1 ~= E1);
res(end+1,1) = E1 ~= E2;
res(end+1,1) = ~ne(E1,E3);

% different dimensions
E4 = ellipsoid([5 -8 -8; -8 16 16; -8 16 20]);

res(end+1,1) = E1 ~= E4;

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
