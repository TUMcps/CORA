function r = radius(Z)
% radius - computes the radius of a hypersphere enclosing a zonotope
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
%    Z = zonotope([1;-1],[-1 2 1; 3 2 0]);
%    r = radius(Z)
%
% Other m-files required: ---
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Authors:       Matthias Althoff
% Written:       19-April-2010
% Last update:   27-July-2016
%                27-August-2019 (MW)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%extract generators
G = Z.G;

%method 1
%add length of generators
r = sum(vecnorm(G));

%method 2
%convert to interval (axis-aligned box around zonotope)
IH=interval(Z);
%compute half of edge length
l=rad(IH);
%compute enclosing radius
rAlt=norm(l);

%choose minimum
r=min(r,rAlt);

% ------------------------------ END OF CODE ------------------------------
