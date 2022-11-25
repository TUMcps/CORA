function res = isempty(E)
% isempty - checks if ellipsoid is empty
%
% Syntax:  
%    res = isempty(E)
%
% Inputs:
%    E - ellipsoid
%
% Outputs:
%    res - boolean whether ellipsoid is empty or not
%
% Example: 
%    E = ellipsoid([1 0; 0 1], [0.5; -1]);
%    res = isempty(E)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      16-Sep-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = isempty(E.Q) && isempty(E.q);

%------------- END OF CODE --------------