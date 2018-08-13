function [Zrem] = cluster(obj,options,Z)
% cluster - clusters and unifies zonotopes
%
% Syntax:  
%    [Zrem] = cluster(Zdummy,Z,dir)
%
% Inputs:
%    Z - cell array of zonotope objects
%    dir - direction of the zonotopes/vector field
%
% Outputs:
%    Zrem - cell array of remaining zonotope objects
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author: Matthias Althoff
% Written: 09-October-2008
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%initialize Zrem
Zrem=[];

%obtain centers of the parallelotopes
for i=1:length(Z)
    c(i,:)=center(Z{i});
end

clusters=5;
if clusters>length(Z)
    clusters=length(Z);
end

[IDX,C] = kmeans(c,clusters);

for i=1:clusters
    %init vertices
    V{i}=[];
end

%unify vertices
for i=1:length(Z)
    %obtain cluster number
    cluster=IDX(i);
    %add vertices
    Vadd=extreme(polytope(Z{i}))';
    V{cluster}=[V{cluster}, Vadd];
end

%enclose vertices by zonotopes
for i=1:clusters
    
    %obtain actual vertices
    Vact=vertices(V{i});
    
    %obtain direction of the vector field
    x0=center(interval(Vact));
    fcnHandle=getfcn(obj.contDynamics,options);
    direction=fcnHandle(0,x0); 
        
    Zrem{i}=enclosingZonotope3(Vact,direction);
end

%------------- END OF CODE --------------