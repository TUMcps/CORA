function res = test_unitTestsStatus
% test_unitTestsStatus - checks if all results stored in unitTestsStatus.mat
%    are transferred to ./unitTests/ci/results
%
% Syntax:
%    res = test_unitTestsStatus
%
% Inputs:
%    no
%
% Outputs:
%    res - true/false 
%
% Example: 
%

% Authors:       Tobias Ladner
% Written:       21-May-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = isempty(which('unitTestsStatus.mat'));

end

% ------------------------------ END OF CODE ------------------------------
