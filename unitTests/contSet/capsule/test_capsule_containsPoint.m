function res = test_capsule_containsPoint
% test_capsule_containsPoint - unit test function of containsPoint
%
% Syntax:  
%    res = test_capsule_containsPoint
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
dim = 3;
C = capsule([1; 1; 0], [0.5; -1; 1], 0.5);

% Single points -----------------------------------------------------------
% create a point inside (center)
p_inside = center(C);
% ... and outside
p_outside = 10*ones(dim,1);

% check if correct results for containment
res_inside  = containsPoint(C, p_inside);
res_outside = containsPoint(C, p_outside);

% compare results
res_single = res_inside && ~res_outside;
% -------------------------------------------------------------------------

% Array of points ---------------------------------------------------------
% all outside...
num = 10;
p_array = 10*(ones(dim,num)+rand(dim,num));
res_array = any(containsPoint(C, p_array));
% -------------------------------------------------------------------------

res = res_single && ~res_array;

if res
    disp('test_containsPoint successful');
else
    disp('test_containsPoint failed');
end

%------------- END OF CODE --------------