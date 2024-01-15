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
%    Z = zonotope([1 1 0; 0 0 1]);
%    p = polygon(Z);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Daniel Heß, Matthias Althoff
% Written:       28-June-2016
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty case
if representsa_(Z,'emptySet',eps)
    p = zeros(dim(Z),0); return
end

% only supported for 2D cases
if dim(Z) ~= 2
    throw(CORAerror('CORA:notSupported','Only supports for 2D zonotopes.'));
end

% delete zero generators
Z = compact_(Z,'zeros',eps);

% obtain center and generator matrix
c = Z.c;
G = Z.G;

% obtain number of generators
nrGens = size(G,2);

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
p = zeros(2,nrGens+1);
for i = 1:nrGens
    p(:,i+1) = p(:,i) + 2*Gnorm(:,IX(i));
end

p(1,:) = p(1,:) + xmax - max(p(1,:));
p(2,:) = p(2,:) - ymax;

%flip/mirror upper half to get lower half of zonotope (point symmetry)            
p = [p(1,:),p(1,end)+p(1,1)-p(1,2:end);...
    p(2,:),p(2,end)+p(2,1)-p(2,2:end)];

%consider center
p = c + p;

% ------------------------------ END OF CODE ------------------------------
