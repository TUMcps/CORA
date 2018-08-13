function S = project(h, S)
% project - projects a set onto the halfspace
%
% Syntax:  
%    S = project(h, S)
%
% Inputs:
%    h - halfspace object
%    S - set to be projected
%
% Outputs:
%    S - projected set
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
% Written:      04-September-2013
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%obtain dimension
dim = length(h.c);

%obtain linear map
M = eye(dim) - h.c*h.c.'/(h.c.'*h.c);
n = (h.d*h.c)/(h.c.'*h.c);

%compute affine map
S = M*S + n;


% %create unit vector
% uVec = [1;zeros(dim-1,1)];
% 
% %obtain rotation matrix
% rotMat = rotationMatrix(h, uVec);
% 
% %rotate set and halfspace
% S = rotMat*S;
% h = rotMat*h;
% 
% %project onto first coordinate
% M = diag([0;ones(dim-1,1)]);
% dist = h.d/norm(h.c);
% S = M*S + [dist;zeros(dim-1,1)];
% 
% %rotate back
% S = rotMat.'*S;


%------------- END OF CODE --------------