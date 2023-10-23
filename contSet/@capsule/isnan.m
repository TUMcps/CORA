function res = isnan(C)
% isnan - checks if any value in the capsule is NaN
%
% Syntax:
%    res = isnan(C)
%
% Inputs:
%    C - capsule object
%
% Outputs:
%    res - false
%
% Example: 
%    C = capsule([1; 1; 0], [0.5; -1; 1], 0.5);
%    res = isnan(C)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       01-June-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% NaN values are not possible by constructor
res = false;

% ------------------------------ END OF CODE ------------------------------
