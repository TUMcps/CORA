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

% Author:       Niklas Kochdumper, Mark Wetzlinger
% Written:      02-January-2020 
% Last update:  12-March-2021 (MW, empty interval)
% Last revision:---

%------------- BEGIN CODE --------------

if isempty(obj)
    res = false;
else
    r = rad(obj);   
    res = ~any(r == 0);
end

%------------- END OF CODE --------------