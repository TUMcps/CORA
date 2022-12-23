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

% instantiate capsule
C = capsule([0;0],[1;0],0.5);

% compute enclosing radius
C_rad = radius(C);
% true solution
C_rad_true = 1.5;

% compare results
res = withinTol(C_rad,C_rad_true);

%------------- END OF CODE --------------