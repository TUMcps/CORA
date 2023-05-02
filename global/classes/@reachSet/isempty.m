function res = isempty(R)
% isempty - checks if a reachSet object is empty
%
% Syntax:  
%    res = isempty(R)
%
% Inputs:
%    R - reachSet object
%
% Outputs:
%    res - true/false
%
% Example:
%    R = reachSet();
%    isempty(R)
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

res = all(arrayfun(@(x) isempty(x.timePoint),R,'UniformOutput',true));

%------------- END OF CODE --------------
