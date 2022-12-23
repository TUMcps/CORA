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
C_add = capsule([0; 1], [sqrt(2); sqrt(2)], 1);
C_caps = C + C_add;

% true solution
C_vect_true = capsule([2; 1], [-4; 3], 1);
C_caps_true = capsule([1; 1], [-4; 3], 4);

% compare solutions
res_vect(1) = compareMatrices(C_vect.c,C_vect_true.c);
res_vect(2) = compareMatrices(C_vect.g,C_vect_true.g);
res_vect(3) = withinTol(C_vect.r,C_vect_true.r);
res_caps(1) = compareMatrices(C_caps.c,C_caps_true.c);
res_caps(2) = compareMatrices(C_caps.g,C_caps_true.g);
res_caps(3) = withinTol(C_caps.r,C_caps_true.r);

% empty set
res_e = isempty(C_add + capsule());

% add results
res = all([res_vect, res_caps, res_e]);

%------------- END OF CODE --------------