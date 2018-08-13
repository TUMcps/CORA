function [P] = multiPolytope(firstZ,Z)
% multiPolytope - Converts a set of zonotope from a G- to a 
% H-representation and intersects them
%
% Syntax:  
%    [P] = multiPolytope(firstZ,Z)
%
% Inputs:
%    firstZ - first zonotope objects
%    Z - cell array of zonotope objects
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

% Author:       Matthias Althoff
% Written:      30-September-2008
% Last update:  26-February-2009
%               30-July-2016
% Last revision:---

%------------- BEGIN CODE --------------

%compute first polytope
[P] = polytope(firstZ);

%intersect with other polytopes
for i=2:length(Z)
    Pnew = polytope(Z{i});
    P=P&Pnew; %intersection
end

%------------- END OF CODE --------------