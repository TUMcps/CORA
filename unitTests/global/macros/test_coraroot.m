function res = test_coraroot
% test_coraroot - tests the macro CORAROOT
%
% Syntax:
%    res = test_coraroot()
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
% See also: CORAROOT

% Authors:       Tobias Ladner
% Written:       22-January-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
 
% check if can be executed
root = CORAROOT;

% check if is folder
assert(isfolder(root))

% check if contains expected subfolders
subfolders = {dir(CORAROOT).name};
assert(ismember('contSet',subfolders))
assert(ismember('contDynamics',subfolders))
assert(ismember('nn',subfolders))
assert(ismember('unitTests',subfolders))


% test completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
