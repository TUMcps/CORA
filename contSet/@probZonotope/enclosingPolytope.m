function [P] = enclosingPolytope(pZ)
% enclosingPolytope - Converts the mean of a probabilistic zonotope to a polytope 
% representation
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
% Written: 18-September-2007
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%reduce probabilistic zonotope
[P]=reduce(zonotope(pZ.Z),'best');

%------------- END OF CODE --------------