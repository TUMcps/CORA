function res = isempty(E)
% isempty - checks if an ellipsoid is the empty set
%
% Syntax:  
%    res = isempty(E)
%
% Inputs:
%    E - ellipsoid
%
% Outputs:
%    res - true/false
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

% Author:       Mark Wetzlinger, Victor Gassmann
% Written:      16-Sep-2019
% Last update:  04-July-2022 (VG: added class-array support)
% Last revision:---

%------------- BEGIN CODE --------------

res = false(size(E));
for i=1:numel(E)
    res(i) = isempty(E(i).Q) && isempty(E(i).q);
end

%------------- END OF CODE --------------