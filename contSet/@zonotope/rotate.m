function Z = rotate(Z,dims,angle)
% rotates - rotates a zonotope projected on two coordinates with the
%    specified angle
%
% Syntax:
%    Z = rotate(Z,dims,angle)
%
% Inputs:
%    Z - zonotope object
%    dims - projected dimensions
%    angle - rotation angle (in rad)
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Matthias Althoff
% Written:       07-October-2008
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input arguments
inputArgsCheck({{Z,'att','zonotope'};
                {dims,'att','numeric',{'nonnan','vector','nonnegative'}};
                {angle,'att','numeric',{'nonnan','scalar'}}});

%rotation matrix
R = [cos(angle) -sin(angle); sin(angle) cos(angle)];

%rotate points
Z.c(dims,:) = R*Z.c(dims,:);
Z.G(dims,:) = R*Z.G(dims,:);

% ------------------------------ END OF CODE ------------------------------
