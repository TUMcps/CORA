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
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Mark Wetzlinger
% Written:      28-August-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% instantiate capsule = unit circle / ball centered at origin
unit_radius = 1;
C_circle = capsule([0; 0], [0; 0], unit_radius);
C_ball   = capsule([0; 0; 0], [0; 0; 0], unit_radius);

% calculate volume
vol_circle = volume(C_circle);
vol_ball   = volume(C_ball);

% true volume
vol_circle_true = pi * unit_radius^2;
vol_ball_true   = 4 / 3 * pi * unit_radius^3;

% compare results
tol = 1e-9;
res = abs(vol_circle_true - vol_circle) < tol && ...
    abs(vol_ball_true - vol_ball) < tol;

if res
    disp('test_volume successful');
else
    disp('test_volume failed');
end

%------------- END OF CODE --------------