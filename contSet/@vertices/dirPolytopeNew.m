function [Zenclose] = dirPolytopeNew(V,direction)
% polytope - Converts a zonotope to a polytope representation
%
% Syntax:  
%    [P] = polytope(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    P - polytope object
%
% Example: 
%    Z=zonotope(rand(2,5));
%    P=polytope(Z);
%    plot(P);
%    hold on
%    plot(Z);
%
% Other m-files required: vertices, polytope
% Subfunctions: none
% MAT-files required: none
%
% See also: intervalhull,  vertices

% Author: Matthias Althoff
% Written: 11-October-2008
% Last update: 26-November-2008
% Last revision: ---

%------------- BEGIN CODE --------------

dim=length(direction);
orient=eye(dim);

%check if direction is not the origin
if norm(direction)>0
    newGen=direction/norm(direction);
else
    newGen=[1;zeros(dim-1,1)];
end

%retrieve most aligned generator from orient
for iGen=1:length(orient(1,:))
     h(iGen)=abs(newGen'*orient(:,iGen)/norm(orient(:,iGen)));
end

[val,ind]=sort(h);
pickedIndices=ind(1:(end-1));

newP=[newGen,orient(:,pickedIndices)];

%compute oriented rectangular hull
V=get(V,'V');
Vtrans=inv(newP)*V;
Min=min(Vtrans,[],2);
Max=max(Vtrans,[],2);
IHtrans=interval(Min,Max);
Ztrans=zonotope(IHtrans);
Zenclose=newP*Ztrans;
 
%------------- END OF CODE --------------