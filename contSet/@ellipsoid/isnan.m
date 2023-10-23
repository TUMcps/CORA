function res = isnan(E)
% isnan - checks if any value in the ellipsoid is NaN
%
% Syntax:
%    res = isnan(E)
%
% Inputs:
%    E - ellipsoid object
%
% Outputs:
%    res - false
%
% Example: 
%    E = ellipsoid([5 7;7 13],[1;2]);
%    res = isnan(E)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       01-June-2022
% Last update:   04-July-2022 (VG, support class arrays)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = false(size(E));

% ------------------------------ END OF CODE ------------------------------
