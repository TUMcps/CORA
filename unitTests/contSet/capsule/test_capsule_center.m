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
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Mark Wetzlinger
% Written:      27-July-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

% empty capsule
C_empty = capsule();
if ~isempty(center(C_empty))
    res = false;
end
    
% init capsule
c_true = [2; 0; -1];
C = capsule(c_true, [1; -1; 2], 0.5);

% read center
c = center(C);

% check result
if ~withinTol(c,c_true)
    res = false;
end

%------------- END OF CODE --------------