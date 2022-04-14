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
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Matthias Althoff, Mark Wetzlinger
% Written:      26-July-2016
% Last update:  01-May-2020 (MW, add second case)
%               09-September-2020 (MA, approximate computation added)
% Last revision:---

%------------- BEGIN CODE --------------

% create zonotope
Z1 = zonotope([-4, -3, -2, -1; 1, 2, 3, 4]);

%% obtain result
vol = volume(Z1);

% true result 1
true_vol = 80;

res_int(1) = vol == true_vol;

%% compare to interval
I1 = interval(Z1);
volInt = volume(I1);
% convert back to zonotope
IZ1 = zonotope(I1);
volIntZon = volume(IZ1); % has to be equal to interval volume

res_int(2) = vol < volInt;
res_int(3) = volIntZon == volInt;

%% approximate computation
% order reduction
volApprox_red = volume(Z1, 'reduce', 1);
true_vol_approx_red = 122.8162136821466106;
res_int(4) = (abs(volApprox_red - true_vol_approx_red)<1e-10);

% Alamo technique
volApprox_red = volume(Z1, 'alamo');
true_vol_approx_red = 48.9897948556635612;
res_int(5) = (abs(volApprox_red - true_vol_approx_red)<1e-10);

%% final result
res = all(res_int);


if res
    disp('test_zonotope_volume successful');
else
    disp('test_zonotope_volume failed');
end

%------------- END OF CODE --------------
