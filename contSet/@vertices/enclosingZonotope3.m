function [Zenclose] = enclosingZonotope3(V,direction)
% enclosingZonotope3 - Encloses a set of points by a hyperrectangle
%
% Syntax:  
%    [Zenclose] = enclosingZonotope3(V,direction)
%
% Inputs:
%    V - vertices object
%    direction - direction of the vector field
%
% Outputs:
%    Zenclose - enclosing zonotope
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author: Matthias Althoff
% Written: 07-October-2008
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%get dimension
dim=length(direction);
Vorig=V;

%convert direction to a vertices object
dirVec=vertices(direction);

%obtain rotation angles
for i=1:(dim-1)
    %project
    tmpDir=get(dirVec,'V');
    dirProj=tmpDir(dim-i:dim-i+1);
    
    %get angle to projected x1-axis
    angle(dim-i)=sign(dirProj(2))*acos(dirProj(1)/norm(dirProj));
    
    %rotate dirVec
    dirVec = rotate(dirVec,dim-i:dim-i+1,-angle(dim-i));
end

%rotate vertices 
for i=1:(dim-1)
    %rotate vertices for the specified dimensions
    V = rotate(V,dim-i:dim-i+1,-angle(dim-i));
end

%obtain enclosing axis aligned box
Min=min(get(V,'V'),[],2);
Max=max(get(V,'V'),[],2);
IH=intervalhull([Min,Max]);
Zenclose=zonotope(IH);

%rotate zonotope to the final position
for i=1:(dim-1)
    %rotate vertices for the specified dimensions
    Zenclose = rotate(Zenclose,i:i+1,angle(i));
end


% %check result by plotting a direction "arrow"
% o=vertices([1;zeros(dim-1,1)]);
% %rotate zonotope to the final position
% for i=1:(dim-1)
%     %rotate vertices for the specified dimensions
%     o = rotate(o,i:i+1,angle(i));
% end
% 
% o=get(o,'V');
% p=center(intervalhull(Vorig));
% Z=zonotope([p,o]);


    
%------------- END OF CODE --------------