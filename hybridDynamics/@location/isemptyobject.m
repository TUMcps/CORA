function res = isemptyobject(loc)
% isemptyobject - checks if a location object is empty
%
% Syntax:
%    res = isemptyobject(loc)
%
% Inputs:
%    loc - location object
%
% Outputs:
%    res - true/false
%
% Example: 
%    res = isemptyobject(location())
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       16-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

[r,c] = size(loc);
res = false(r,c);

% loop over all locations
for i=1:r
    for j=1:c
        res(i,j) = isemptyobject(loc(i,j).contDynamics);
    end
end

% ------------------------------ END OF CODE ------------------------------
