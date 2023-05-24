function res = isempty(loc)
% isempty - checks if a location object is empty
%
% Syntax:  
%    res = isequal(loc)
%
% Inputs:
%    loc - location object
%
% Outputs:
%    res - true/false
%
% Example: 
%    res = isempty(location())
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      16-May-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

[r,c] = size(loc);
res = false(r,c);

% loop over all locations
for i=1:r
    for j=1:c
        res(i,j) = isempty(loc(i,j).contDynamics);
    end
end

%------------- END OF CODE --------------