function res = test_capsule_isempty
% test_capsule_isempty - unit test function of isempty
%
% Syntax:  
%    res = test_capsule_isempty
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

% instantiate capsules
C1 = capsule();
C2 = capsule([1;1],[],0.5);
C3 = capsule([1;1],[0;1],0.5);
C4 = capsule([1;1],[0;1],0);

% compare results
res = isempty(C1) && isempty(C2) && ~isempty(C3) && ~isempty(C4);

if res
    disp('test_isempty successful');
else
    disp('test_isempty failed');
end

%------------- END OF CODE --------------