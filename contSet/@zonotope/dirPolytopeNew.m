function [Zenclose] = dirPolytopeNew(Z,direction)
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
% Written: 14-October-2008
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

dim=length(direction);
orient=eye(dim);

newGen=direction/norm(direction);

%retrieve most aligned generator from orient
for iGen=1:length(orient(1,:))
     h(iGen)=abs(newGen'*orient(:,iGen)/norm(orient(:,iGen)));
end

[val,ind]=sort(h);
pickedIndices=ind(1:(end-1));

newP=[newGen,orient(:,pickedIndices)];

%compute oriented rectangular hull
Ztrans=inv(newP)*Z;
IHtrans=intervalhull(Ztrans);
Ztrans=zonotope(IHtrans);
Zenclose=newP*Ztrans;
 
%------------- END OF CODE --------------