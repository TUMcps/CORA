function p = polygon(V)
% polygon - computes ordered lists of vertices, defining a polygon; only
% the first two coordinates are considered;
%
% Syntax:  
%    polygon(V)
%
% Inputs:
%    V - vertices object 
%
% Outputs:
%    none
%
% Example: 
%    V = vertices(rand(2,6));
%    p = polygon(V)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      18-August-2016 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%check if vertex array is not empty
if ~isempty(V.V)

    %Obtain x- and y- coordinates of potential vertices
    xPotential=V.V(1,:);
    yPotential=V.V(2,:);

    %there are more than two points
    if length(xPotential)>2
        %determine vertex indices for convex hull vertices from potential
        %vertices 
        try
            vertexIndices=convhull(xPotential,yPotential);
        catch
            disp('plot failed');
            vertexIndices=[];
        end
    else
        vertexIndices=[1 2];
    end
    
    if ~isempty(vertexIndices)
        %Select convex hull vertices
        p(1,:) = xPotential(vertexIndices);
        p(2,:) = yPotential(vertexIndices);
    end
end


%------------- END OF CODE --------------