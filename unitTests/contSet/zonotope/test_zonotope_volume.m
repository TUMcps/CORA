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
% See also: -

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       26-July-2016
% Last update:   01-May-2020 (MW, add second case)
%                09-September-2020 (MA, approximate computation added)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty set
res_e = (volume(zonotope()) == 0);

% create zonotope
Z1 = zonotope([-4, -3, -2, -1; 1, 2, 3, 4]);

%% obtain result
vol = volume(Z1);

% true result 1
true_vol = 80;

res_int(1) = withinTol(vol,true_vol);

%% compare to interval
I1 = interval(Z1);
volInt = volume(I1);
% convert back to zonotope
IZ1 = zonotope(I1);
volIntZon = volume(IZ1); % has to be equal to interval volume

res_int(2) = vol < volInt;
res_int(3) = withinTol(volIntZon,volInt);

%% approximate computation
% order reduction
volApprox_red = volume(Z1, 'reduce', 1);
true_vol_approx_red = 122.8162136821466106;
res_int(4) = withinTol(volApprox_red,true_vol_approx_red);

% Alamo technique
volApprox_red = volume(Z1, 'alamo');
true_vol_approx_red = 48.9897948556635612;
res_int(5) = withinTol(volApprox_red,true_vol_approx_red);

%% final result
res = all(res_int) && res_e;

% ------------------------------ END OF CODE ------------------------------
