function V = vertices_(probZ,varargin)
% vertices_ - Returns potential vertices of a probabilistic zonotope
%    WARNING: Do not use this function for high order zonotopes -
%             computational complexity grows exponentially!
%
% Syntax:
%    V = vertices_(probZ)
%
% Inputs:
%    probZ - probZonotope object
%
% Outputs:
%    V - matrix storing the vertices V = [v1,...,vp]
%
% Example: 
%    Z1 = [10 1 -2; 0 1 1];
%    Z2 = [0.6 1.2; 0.6 -1.2];
%    probZ = probZonotope(Z1,Z2,2);
%    V = vertices(probZ);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/vertices, interval,  polytope

% Authors:       Matthias Althoff
% Written:       14-September-2006 
% Last update:   22-March-2007
% Last revision: 27-March-2023 (MW, rename vertices_)

% ------------------------------ BEGIN CODE -------------------------------

%get matrix from object
probZmat=probZ.Z;

%first vertex is the center of the zonotope
vertexArray=probZmat(:,1);

%Generate further potential vertices in the loop
for iVertex=1:length(probZmat(1,2:end))
    translation=probZmat(:,iVertex+1)*ones(1,length(vertexArray(1,:)));
    V=[vertexArray+translation,vertexArray-translation];
    %remove inner points
    if iVertex>length(probZmat(:,1))
        K = convhulln(V');
        indices = unique(K);
        vertexArray=V(:,indices);
    else
        vertexArray=V;
    end
end

%create vertices object
V = vertexArray;

% ------------------------------ END OF CODE ------------------------------
