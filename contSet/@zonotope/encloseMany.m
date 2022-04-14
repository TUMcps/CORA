function [Zenclose,rotMatrixInv] = encloseMany(Z,direction)
% encloseMany - function for the enclosure of many zonotopes
%
% Syntax:  
%    [Zenclose,rotMatrixInv] = encloseMany(Z,direction)
%
% Inputs:
%    Z - cell array of zonotopes to be enclosed
%    direction - mean direction, in which the zonotopes to be enclosed are
%                heading to
%
% Outputs:
%    Zenclose - enclosing zonotope (which is an oriented rectangular hull)
%    rotMatrix - rotation matrix of the oriented rectangular hull
%
% Example: 
%
% Other m-files required: 
% Subfunctions: none
% MAT-files required: none
%
% See also: dirPolytope

% Author:       Matthias Althoff
% Written:      15-January-2008
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%get dimension and original axis aligned orientation
dim=length(direction);
orient=eye(dim);

%replace one of the axis aligned vectors
newGen=direction/norm(direction);

%retrieve most aligned generator from orient
for iGen=1:length(orient(1,:))
     h(iGen)=abs(newGen'*orient(:,iGen)/norm(orient(:,iGen)));
end

[~,ind]=sort(h);
pickedIndices=ind(1:(end-1));

rotMatrix=[newGen,orient(:,pickedIndices)];

%obtain and collect vertices
Vsum=[];
for i=1:length(Z)
    Zred=reduce(Z{i},'parallelpiped');
    Vnew=vertices(Zred);
    Vsum=[Vsum,Vnew];
end

%rotate vertices
rotMatrixInv=inv(rotMatrix);
Vtrans=rotMatrixInv*Vsum;
%compute enclosing interval
Itrans=interval.enclosePoints(vertices(Vtrans));
%generate zonotope
Ztrans=zonotope(Itrans);
%rotate zonotope back
Zenclose=rotMatrix*Ztrans;
    
%------------- END OF CODE --------------