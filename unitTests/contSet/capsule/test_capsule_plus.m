function res = test_capsule_plus
% test_capsule_plus - unit test function of plus
%
% Syntax:  
%    res = test_capsule_plus
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

% instantiate random capsule
C = capsule([1; 0], [-4; 3], 1);

% Minkowski addition: vector
vect = [1; 1];
C_vect = C + vect;

% Minkowski addition: capsule
Cadd = capsule([0; 1], [sqrt(2); sqrt(2)], 1);
C_caps = C + Cadd;

% true solution
C_vect_true = capsule([2; 1], [-4; 3], 1);
C_caps_true = capsule([1; 1], [-4; 3], 4);

% compare solutions
tol = 1e-9;
res_vect(1) = all(abs(center(C_vect) - center(C_vect_true)) < tol);
res_vect(2) = all(abs(C_vect.g - C_vect_true.g) < tol);
res_vect(3) = abs(radius(C_vect) - radius(C_vect_true)) < tol;
res_caps(1) = all(abs(center(C_caps) - center(C_caps_true)) < tol);
res_caps(2) = all(abs(C_caps.g - C_caps_true.g) < tol);
res_caps(3) = abs(radius(C_caps) - radius(C_caps_true)) < tol;

% add results
res = all([res_vect, res_caps]);

if res
    disp('test_plus successful');
else
    disp('test_plus failed');
end

%------------- END OF CODE --------------