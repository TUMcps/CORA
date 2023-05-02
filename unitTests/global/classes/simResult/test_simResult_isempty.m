function res = test_simResult_isempty
% test_simResult_isempty - unit test function for isempty
%
% Syntax:  
%    res = test_simResult_isempty()
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
% See also: none

% Author:       Mark Wetzlinger
% Written:      01-May-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% empty simResult
simRes = simResult();
res = isempty(simRes);

% simResult with trajectory
t = {[0; 0.02; 0.05]};
x = {[1 1; 0.9 1.1; 0.8 1.2]};
simRes = simResult(x,t);
res(end+1,1) = ~isempty(simRes);

% combine results
res = all(res);

%------------- END OF CODE --------------