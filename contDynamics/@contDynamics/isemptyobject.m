function res = isemptyobject(sys)
% isemptyobject - checks if a contDynamics object is empty
%
% Syntax:
%    res = isemptyobject(sys)
%
% Inputs:
%    sys - contDynamics object
%
% Outputs:
%    res - true/false
%
% Example: 
%    ---
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

[r,c] = size(sys);
res = false(r,c);

% loop over all contDynamics
for i=1:r
    for j=1:c
        res(r,c) = sys(r,c).dim == 0;
    end
end

% ------------------------------ END OF CODE ------------------------------
