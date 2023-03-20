function res = isInterval(zB)
% isInterval - checks if a zonotope bundle can be equivalently represented
%    by an interval object
%
% Syntax:  
%    res = isInterval(zB)
%
% Inputs:
%    zB - zonoBundle object
%
% Outputs:
%    res - true/false
%
% Example:
%    list{1} = zonotope([1;0],eye(2));
%    list{2] = zonotope([1;0.5],eye(2));
%    zB = zonoBundle(list);
%    isInterval(zB)
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

% one-dimensional zonoBundles are always intervals
if dim(zB) == 1
    res = true; return
end

% empty sets can also be represented by intervals
if isempty(zB)
    res = true; return
end

% all other cases: check individual zonotopes
% (note: there are cases where the intersection is still an interval)
res = true;
for i=1:zB.parallelSets
    if ~isInterval(zB.Z{1})
        res = false; return
    end
end

%------------- END OF CODE --------------