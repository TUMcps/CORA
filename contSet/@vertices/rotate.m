function [V] = rotate(V,dims,angle)
% rotates - rotates a set of points projected on two coordinates with the
% specified angle
%
% Syntax:  
%    [V] = rotate(V,dims,angle)
%
% Inputs:
%    V - vertices object
%    dims - projected dimensions
%    angle - rotation angle
%
% Outputs:
%    V - vertices object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author: Matthias Althoff
% Written: 06-October-2008
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%rotation matrix
R=[cos(angle) -sin(angle); sin(angle) cos(angle)];

%rotate points
V.V(dims,:)=R*V.V(dims,:);

%------------- END OF CODE --------------