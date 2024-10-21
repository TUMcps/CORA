function res = test_zonotope_volume
% test_zonotope_volume - unit test function of volume
%
% Syntax:
%    res = test_zonotope_volume
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

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       26-July-2016
% Last update:   01-May-2020 (MW, add second case)
%                09-September-2020 (MA, approximate computation added)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty set
Z_empty = zonotope.empty(2);
assert(volume(Z_empty) == 0);

% 2D zonotope
Z1 = zonotope([-4, -3, -2, -1; 1, 2, 3, 4]);
vol = volume(Z1);
true_vol = 80;
assert(withinTol(vol,true_vol));

% compare to interval
I1 = interval(Z1);
volInt = volume(I1);
% convert back to zonotope
IZ1 = zonotope(I1);
volIntZon = volume(IZ1); % has to be equal to interval volume

assert(vol < volInt);
assert(withinTol(volIntZon,volInt));

% approximate computation
% order reduction
volApprox_red = volume(Z1, 'reduce', 1);
true_vol_approx_red = 122.8162136821466106;
assert(withinTol(volApprox_red,true_vol_approx_red));

% Alamo technique
volApprox_red = volume(Z1, 'alamo');
true_vol_approx_red = 48.9897948556635612;
assert(withinTol(volApprox_red,true_vol_approx_red));


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
