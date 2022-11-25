function [c] = center(Z)
% center - Returns the center of a zonotope bundle; the center is defined
% as the center of the last zonotope; alternative: mean of all centers
%
% Syntax:  
%    [c] = center(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    c - center of the zonotope Z
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      03-February-2011
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

c=center(Z.Z{end});

%------------- END OF CODE --------------