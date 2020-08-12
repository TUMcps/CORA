function res = test_capsule_radius
% test_capsule_radius - unit test function of radius
%
% Syntax:  
%    res = test_capsule_radius
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

% instantiate capsules
C = capsule([1;1],[1;0],0.5);
Cball = capsule([0;0],[0;0],1); % is a unit ball

% compute enclosing radius
C_rad = radius(C);
Cball_rad = radius(Cball);
% true solution
C_rad_true = 1.5;
Cball_rad_true = 1;

% compare results
tol = 1e-9;
res = abs(C_rad - C_rad_true) < tol && ...
    abs(Cball_rad - Cball_rad_true) < tol;

if res
    disp('test_radius successful');
else
    disp('test_radius failed');
end

%------------- END OF CODE --------------