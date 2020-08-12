function res = test_capsule_isequal
% test_capsule_isequal - unit test function of isequal
%
% Syntax:  
%    res = test_capsule_isequal
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
% Written:      17-Sep-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% instantiate capsule
C1 = capsule([1;1],[1;1],0.5);
C2 = capsule([1;1],[1;0],0.5);
C3 = C1;

% compare results
tol = 1e-9;
res = ~isequal(C1,C2,tol) && isequal(C1,C3,tol);

if res
    disp('test_isequal successful');
else
    disp('test_isequal failed');
end

%------------- END OF CODE --------------