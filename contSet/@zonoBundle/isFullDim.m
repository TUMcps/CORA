function res = isFullDim(obj)
% isFullDim - check if a zonotope bundle is full-dimensional
%
% Syntax:  
%    res = isFullDim(obj)
%
% Inputs:
%    obj - zonoBundle object
%
% Outputs:
%    res - 1 if zonoBundle is full-dimensional, 0 else
%
% Example:
%    int1 = interval([0;0],[2;2]);
%    int2 = interval([0;1],[2;3]);
%    int3 = interval([0;2],[2;4]);
%
%    zB1 = zonoBundle({zonotope(int1),zonotope(int2)});
%    zB2 = zonoBundle({zonotope(int1),zonotope(int3)});
%
%    isFullDim(zB1)
%    isFullDim(zB2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/isFullDim

% Author:       Niklas Kochdumper
% Written:      02-January-2020 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % loop over all zonotopes
    for i = 1:obj.parallelSets
       if ~isFullDim(obj.Z{i})
          res = 0;
          return;
       end
    end

    % convert to polytope and call polytope/isFullDim
    res = isFullDim(mptPolytope(obj));

%------------- END OF CODE --------------