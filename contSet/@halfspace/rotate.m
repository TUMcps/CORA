function h = rotate(h, newDir, rotPoint)
% rotate - rotates a halfspace around a rotation point rotPoint such that
% the new normal vector is aligned with newDir
%
% Syntax:  
%    h = rotate(h, newDir, rotPoint)
%
% Inputs:
%    h - halfspace object
%    newDir - vector pointing in the new direction
%    rotPoint - rotation point
%
% Outputs:
%    h - halfspace object
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Matthias Althoff
% Written:      28-August-2013
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%obtain rotation matrix
rotMat = rotationMatrix(h, newDir);

%translate and rotate halfspace
h = rotMat*(h + (-rotPoint)) + rotPoint;


%------------- END OF CODE --------------