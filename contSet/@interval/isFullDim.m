function res = isFullDim(obj)
% isFullDim - check if an interval is full-dimensional
%
% Syntax:  
%    res = isFullDim(obj)
%
% Inputs:
%    obj - interval object
%
% Outputs:
%    res - 1 if interval is full-dimensional, 0 else
%
% Example:
%    int1 = interval([-1;-2],[1;2]);
%    int2 = interval([-1;-2],[1;-2]);
%
%    isFullDim(int1)
%    isFullDim(int2)
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

r = rad(obj);   
res = ~any(r == 0);

%------------- END OF CODE --------------