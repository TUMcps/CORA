function res = test_capsule_center
% test_capsule_center - unit test function of center
%
% Syntax:  
%    res = test_capsule_center
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
C = capsule([1; 1; 0], [0.5; -1; 1], 0.5);

% read center
c = center(C);

% true solution
c_true = [1;1;0];

% compare results
res = all(c == c_true);

if res
    disp('test_center successful');
else
    disp('test_center failed');
end

%------------- END OF CODE --------------