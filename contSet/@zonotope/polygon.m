function p = polygon(Z)
% polygon - converts a two-dimensional zonotope into a polygon and returns
%    its vertices
%
% Syntax:  
%    p = polygon(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    p - ordered set of points of a polygon
%
% Example: 
%    Z=zonotope([1 1 0; 0 0 1]);
%    p=polygon(Z);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Daniel Heß, Matthias Althoff
% Written:      28-June-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% delete zero generators
Z = deleteZeros(Z);

% obtain center and generator
c = center(Z); %center
G = generators(Z); %matrix of generators

% obtain number of generators
n = size(G,2);

% obtain size of enclosing intervalhull of first two dimensions
xmax = sum(abs(G(1,:)));
ymax = sum(abs(G(2,:)));

 
% Z with normalized direction: All generators pointing "up"
Gnorm = G;
Gnorm(:,G(2,:)<0)=G(:,G(2,:)<0)*-1;

%compute angles
angles = atan2(Gnorm(2,:),Gnorm(1,:));
angles(angles<0) = angles(angles<0) +2*pi;%handle numerical imprecision/deficiency in atan2, wraparound is not at pi?!?

% assert(~any(angles>pi));

%sort all generators by their angle
[~,IX] = sort(angles,'ascend');

%cumsum the generators in order of angle
p = zeros(2,n+1);
for i = 1:n
    p(:,i+1) = p(:,i) + 2*Gnorm(:,IX(i));
end

p(1,:) = p(1,:) + xmax - max(p(1,:));
p(2,:) = p(2,:) - ymax;

%flip/mirror upper half to get lower half of zonotope (point symmetry)            
p = [p(1,:),p(1,end)+p(1,1)-p(1,2:end);...
    p(2,:),p(2,end)+p(2,1)-p(2,2:end)];

%consider center
p(1,:) = c(1) + p(1,:);
p(2,:) = c(2) + p(2,:);

%------------- END OF CODE --------------