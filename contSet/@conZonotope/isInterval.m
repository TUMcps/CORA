function res = isInterval(cZ)
% isInterval - checks if a constrained zonotope can be equivalently
%    represented by an interval object
%
% Syntax:  
%    res = isInterval(cZ)
%
% Inputs:
%    cZ - conZonotope object
%
% Outputs:
%    res - true/false
%
% Example:
%    cZ1 = conZonotope(ones(2,1),eye(2));
%    isInterval(cZ1)
%
%    cZ2 = conZonotope(ones(2,1),[1 0.4; 0 -1],[1 -0.5],0);
%    isInterval(cZ2)
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
res = false;

% one-dimensional capsule are always intervals
if dim(cZ) == 1
    res = true; return
end

% empty sets can also be represented by intervals
if isemptyobject(cZ)
    res = true; return
end
% note: actual isempty not evaluated for speed reasons

if isempty(cZ.A)
    % check zonotope method
    res = isInterval(zonotope(cZ.Z)); return
end

% ...

%------------- END OF CODE --------------