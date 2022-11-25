function res = test_capsule_enlarge
% test_capsule_enlarge - unit test function of enlarge
%
% Syntax:  
%    res = test_capsule_enlarge
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
% Written:      15-Sep-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% instantiate capsule
C = capsule([1;1],[1;1],0.5);
factor = 2;

% compute enlarged capsule
C_enlarged = enlarge(C,factor);

% true solution
C_true = capsule([1;1],[2;2],1);

% compare results
tol = 1e-9;
res = isequal(C_true,C_enlarged,tol);

if res
    disp('test_enlarge successful');
else
    disp('test_enlarge failed');
end

%------------- END OF CODE --------------