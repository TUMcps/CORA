function res = isInterval(E)
% isInterval - checks if an ellipsoid can be equivalently represented by an
%    interval object
%
% Syntax:  
%    res = isInterval(E)
%
% Inputs:
%    E - ellipsoid object
%
% Outputs:
%    res - true/false
%
% Example:
%    E1 = ellipsoid(0.02);
%    isInterval(E1)
% 
%    E2 = ellipsoid(eye(2),[1;0]);
%    isInterval(E2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      17-March-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% one-dimensional capsule are always intervals
if dim(E) == 1
    res = true; return
end

% empty sets can also be represented by intervals
if isempty(E)
    res = true; return
end

% all other cases cannot be intervals
res = false;

%------------- END OF CODE --------------