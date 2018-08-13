function plot3d(V)
% plot3d - Plots the convex hull of vertices in 3D
%
% Syntax:  
%    plot(V)
%
% Inputs:
%    V - vertices object 
%
% Outputs:
%    none
%
% Example: 
%    V=vertices(rand(3,6));
%    plot3d(V)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      06-October-2007
% Last update:  31-July-2018
% Last revision:---

%------------- BEGIN CODE --------------

%transpose matrix of vertices
vTrans = V.V';

%Plot if there are more than three points
if length(vTrans(:,1))>3
    %determine vertex indices for convex hull vertices from potential
    %vertices 
    vertexIndices=convhulln(vTrans);
else
    vertexIndices=[1 2 3];
end

%plot patches
patch('Faces',vertexIndices,'Vertices',vTrans,'FaceColor','red','FaceAlpha', 0.75);

%------------- END OF CODE --------------