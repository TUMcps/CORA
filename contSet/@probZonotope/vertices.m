function [V] = vertices(Z)
% vertices - Returns potential vertices of a zonotope
% WARNING: Do not use this function for high order zonotopes -
% computational complexity grows exponential!
%
% Syntax:  
%    [V] = vertices(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    V - vertices object
%
% Example: 
%    Z=zonotope(rand(2,5));
%    V=vertices(Z)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: intervallhull,  polytope

% Author: Matthias Althoff
% Written: 14-September-2006 
% Last update: 22-March-2007
% Last revision: ---

%------------- BEGIN CODE --------------

%get matrix from object
Z=Z.Z;

%first vertex is the center of the zonotope
vertexArray=Z(:,1);

%Generate further potential vertices in the loop
for iVertex=1:length(Z(1,2:end))
    translation=Z(:,iVertex+1)*ones(1,length(vertexArray(1,:)));
    V=[vertexArray+translation,vertexArray-translation];
    %remove inner points
    if iVertex>length(Z(:,1))
        K = convhulln(V');
        indices = unique(K);
        vertexArray=V(:,indices);
    else
        vertexArray=V;
    end
end

%create vertices object
V=vertices(vertexArray);

%------------- END OF CODE --------------