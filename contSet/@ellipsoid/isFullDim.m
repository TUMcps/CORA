function res = isFullDim(obj)
% isFullDim - check if an ellipsoid is full-dimensional
%
% Syntax:  
%    res = isFullDim(obj)
%
% Inputs:
%    obj - ellipsoid object
%
% Outputs:
%    res - 1 if ellipsoid is full-dimensional, 0 else
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

    res = ~obj.isdegenerate;

%------------- END OF CODE --------------