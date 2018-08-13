function r = radius(Z)
% radius - Computes the radius of a hypersphere enclosing a zonotope
%
% Syntax:  
%    r = radius(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    r - radius
%
% Example: 
%    Z=zonotope(rand(2,5));
%    r=radius(Z);
%
% Other m-files required: ---
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Matthias Althoff
% Written:      19-April-2010
% Last update:  27-July-2016
% Last revision:---

%------------- BEGIN CODE --------------

%extract generators
G=Z.Z(:,2:end);

%method 1
nrOfGenerators=length(G(1,:));
r=0;
%add length of generators
for i=1:nrOfGenerators
    rPart=norm(G(:,i));
    r=r+rPart;
end

%method 2
%compute enclosing interval hull
IH=intervalhull(Z);
%compute edge length
l=rad(IH);
%compute enclosing radius
rAlt=sqrt(l'*l);

%chose minimum
r=min(r,rAlt);

%------------- END OF CODE --------------