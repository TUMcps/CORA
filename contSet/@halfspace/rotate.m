function hs = rotate(hs, newDir, rotPoint)
% rotate - rotates a halfspace around a rotation point rotPoint such that
%    the new normal vector is aligned with newDir
%
% Syntax:
%    hs = rotate(hs, newDir, rotPoint)
%
% Inputs:
%    hs - halfspace object
%    newDir - vector pointing in the new direction
%    rotPoint - rotation point
%
% Outputs:
%    hs - halfspace object
%
% Example: 
%    hs = halfspace([1 -1],2);
%    newDir = [1;1];
%    rotPoint = [0;0];
% 
%    rotate(hs,newDir,rotPoint)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Authors:       Matthias Althoff
% Written:       28-August-2013
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%obtain rotation matrix
rotMat = rotationMatrix(hs, newDir);

%translate and rotate halfspace
hs = rotMat*(hs + (-rotPoint)) + rotPoint;

% ------------------------------ END OF CODE ------------------------------
