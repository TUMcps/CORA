function r = radius(C)
% radius - Returns the radius of enclosing hyperball of a capsule
%
% Syntax:  
%    r = radius(C)
%
% Inputs:
%    C - capsule
%
% Outputs:
%    r - radius of enclosing hyperball of the capsule C
%
% Example: 
%    C = capsule([1; 1; 0], [0.5; -1; 1], 0.5);
%    r = radius(C)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      28-Aug-2019 
% Last update:  16-Sep-2019
% Last revision:---

%------------- BEGIN CODE --------------

r = norm(C.g,2) + C.r;

%------------- END OF CODE --------------