function d = dim(obj)
% dim - Returns the defined dimension of a conZonotope
%
% Syntax:  
%    d = dim(obj)
%
% Inputs:
%    obj - conZonotope object
%
% Outputs:
%    d - dimension of obj
%
% Example: 
%    Z = [0 1 0 1;0 1 2 -1];
%    A = [-2 1 -1];
%    b = 2;
%    cZono = conZonotope(Z,A,b);
% 
%    d = dim(cZono)
%
% Other m-files required: center.m
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      15-Sep-2019 
% Last update:  14-March-2021 (MW, different approach)
% Last revision:---

%------------- BEGIN CODE --------------

if ~isempty(obj.Z)
    d = size(obj.Z,1);
else
    d = 0;
end

%------------- END OF CODE --------------