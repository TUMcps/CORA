function c = center(obj)
% center - Returns the center of an matZonotope
%
% Syntax:  
%    c = center(obj)
%
% Inputs:
%    obj - matZonotope object
%
% Outputs:
%    c - center of the matZonotope obj
%
% Example:
%    M = matZonotope(eye(2),{eye(2),2*eye(2)});
%    c = center(M)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Victor Gassmann
% Written:      23-July-2020 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

c = obj.center;

%------------- END OF CODE --------------