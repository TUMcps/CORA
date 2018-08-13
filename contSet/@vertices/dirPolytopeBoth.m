function [Zenclose] = dirPolytopeBoth(V,direction1,direction2)
% polytope - Converts a zonotope to a polytope representation
%
% Syntax:  
%    [Zenclose] = dirPolytopeBoth(V,direction1,direction2)
%
% Inputs:
%    V - vertices object
%    direction1 - one of the directions the enclosure should face
%    direction2 - one of the directions the enclosure should face
%
% Outputs:
%    Zenclose - enclosing zonotope
%
% Example: 
%       ---
%
% Other m-files required: vertices, polytope
% Subfunctions: none
% MAT-files required: none
%
% See also: dirPolytopeNew

% Author:       Matthias Althoff
% Written:      14-February-2011
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

dim=length(direction1);
orient=eye(dim);

%check if direction is not the origin for direction1
if norm(direction1)>0
    newGen1=direction1/norm(direction1);
else
    newGen1=[1;zeros(dim-1,1)];
end

%check if direction is not the origin for direction2
if norm(direction2)>0
    newGen2=direction2/norm(direction2);
else
    newGen2=[1;zeros(dim-1,1)];
end

%retrieve most aligned generator from orient
for iGen=1:length(orient(1,:))
     h1(iGen)=abs(newGen1'*orient(:,iGen)/norm(orient(:,iGen)));
     h2(iGen)=abs(newGen2'*orient(:,iGen)/norm(orient(:,iGen)));
end

[val1,ind1]=sort(h1);
[val2,ind2]=sort(h2);

%intersect indices
pickedIndices = intersect(ind1(1:(end-1)),ind2(1:(end-1)));

if length(pickedIndices) <= (dim-2)
    newP=[newGen1, newGen2, orient(:,pickedIndices)];
else
    newP=[newGen1, orient(:,pickedIndices)];
end

%compute oriented rectangular hull
V=get(V,'V');
Vtrans=pinv(newP)*V;
Min=min(Vtrans,[],2);
Max=max(Vtrans,[],2);
IHtrans=intervalhull([Min,Max]);
Ztrans=zonotope(IHtrans);
Zenclose=newP*Ztrans;
 
%------------- END OF CODE --------------