function res = isInterval(C)
% isInterval - checks if a capsule can be equivalently represented by an
%    interval object
%
% Syntax:  
%    res = isInterval(C)
%
% Inputs:
%    C - capsule object
%
% Outputs:
%    res - true/false
%
% Example:
%    C1 = capsule(ones(2,1));
%    isInterval(C1)
%
%    C2 = capsule(ones(2,1),zeros(2,1),0.5);
%    isInterval(C2)
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

% assume true
res = true;

% one-dimensional capsule are always intervals
if dim(C) == 1
    return
end

% empty sets can also be represented by intervals
if isempty(C)
    return
end

% case >=2D: generator must be all-zero and radius = 0
if any(C.g) || C.r > 0
    res = false;
end

%------------- END OF CODE --------------