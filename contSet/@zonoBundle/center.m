function c = center(zB)
% center - Returns the center of the last zonotope of a zonotope bundle
%
% Syntax:  
%    c = center(zB)
%
% Inputs:
%    zB - zonoBundle object
%
% Outputs:
%    c - center
%
% Example:
%    Z1 = zonotope(zeros(2,1),[1 0.5; -0.2 1]);
%    Z2 = zonotope(ones(2,1),[1 -0.5; 0.2 1]);
%    zB = zonoBundle({Z1,Z2});
%    c = center(zB)
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

c = center(zB.Z{end});

%------------- END OF CODE --------------