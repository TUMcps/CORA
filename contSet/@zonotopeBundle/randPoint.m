function [p] = randPoint(obj)
% randPoint - generates a random point within a zonotope bundle
%
% Syntax:  
%    [p] = randPoint(obj)
%
% Inputs:
%    obj - zonotope bundle object
%
% Outputs:
%    p - random point in R^n
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      18-August-2016 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% obtain vertices
V = vertices(obj);

% obtain number of vertices
Vmat = get(V,'V');
nrOfVertices = length(Vmat(1,:));

% random canvex combination
alpha = rand(nrOfVertices,1);
alphaNorm = alpha/sum(alpha);

% random point
p = 0*Vmat(:,1);
for i = 1:nrOfVertices
    p = p + alphaNorm(i)*Vmat(:,i);
end


%------------- END OF CODE --------------