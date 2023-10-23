function res = isnan(I)
% isnan - checks if any value in the interval is NaN
%
% Syntax:
%    res = isnan(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    res - false
%
% Example: 
%    I = interval([-2;-1],[2;1]);
%    res = isnan(I)
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
