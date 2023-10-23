function res = isnan(hs)
% isnan - checks if any value in the halfspace is NaN
%
% Syntax:
%    res = isnan(hs)
%
% Inputs:
%    hs - halfspace object
%
% Outputs:
%    res - false
%
% Example: 
%    hs = halfspace([-1;-1],0);
%    res = isnan(hs)
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
