function res = isemptyobject(R)
% isemptyobject - checks if a reachSet object is empty
%
% Syntax:
%    res = isemptyobject(R)
%
% Inputs:
%    R - reachSet object
%
% Outputs:
%    res - true/false
%
% Example:
%    R = reachSet();
%    isemptyobject(R)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       01-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = all(arrayfun(@(x) isempty(x.timePoint),R,'UniformOutput',true));

% ------------------------------ END OF CODE ------------------------------
