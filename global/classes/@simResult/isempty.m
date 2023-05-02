function res = isempty(simRes)
% isempty - checks if a simResult object is empty
%
% Syntax:  
%    res = isempty(simRes)
%
% Inputs:
%    simRes - simResult object
%
% Outputs:
%    res - true/false
%
% Example:
%    simRes = simResult();
%    isempty(simRes)
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

% check time
res = isempty(simRes.t);

%------------- END OF CODE --------------