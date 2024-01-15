function res = test_capsule_volume
% test_capsule_volume - unit test function of volume
%
% Syntax:
%    res = test_capsule_volume
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
% Written:       27-July-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% empty case
C = capsule.empty(2);
res(end+1,1) = volume(C) == 0;

% instantiate capsule as ball
n = 3;
c = [3; 2; -1];
g = [0; 0; 0];
r = 0.5;
C = capsule(c,g,r);

% calculate volume
vol = volume(C);
% true volume (n-dim sphere)
vol_true = (pi^(n/2) / gamma(1+n/2)) * r^n;

% compare results
tol = 1e-9;
res(end+1,1) = withinTol(vol_true,vol,tol);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
