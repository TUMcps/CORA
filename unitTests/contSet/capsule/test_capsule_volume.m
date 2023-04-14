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
% See also: -

% Author:       Mark Wetzlinger
% Written:      27-July-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

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
res = withinTol(vol_true,vol,tol);
res = res && volume(capsule()) == 0;

%------------- END OF CODE --------------