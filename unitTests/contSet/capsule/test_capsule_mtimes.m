function res = test_capsule_mtimes
% test_capsule_mtimes - unit test function of mtimes
%
% Syntax:  
%    res = test_capsule_mtimes
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
C = capsule([1; 0], [1; -1], 2);

% linear map: 90° = pi/2
angle = pi/2;
R=[cos(angle) -sin(angle); sin(angle) cos(angle)];
C_mtimes = R * C;

% true solution
C_mtimes_true = capsule([0; 1], [1; 1], 2);

% compare results
tol = 1e-9;
res(1) = all(abs(center(C_mtimes) - center(C_mtimes_true)) < tol);
res(2) = all(abs(C_mtimes.g - C_mtimes_true.g) < tol);
res(3) = abs(radius(C_mtimes) - radius(C_mtimes_true)) < tol;

% add results
res = all(res);

if res
    disp('test_mtimes successful');
else
    disp('test_mtimes failed');
end

%------------- END OF CODE --------------