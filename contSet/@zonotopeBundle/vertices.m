function [V] = vertices(obj)
% vertices - Returns potential vertices of a zonotope bundle
% WARNING: Do not use this function for high order zonotope bundles -
% computational complexity grows exponential!
%
% Syntax:  
%    [V] = vertices(Zbundle)
%
% Inputs:
%    Zbundle - zonotope bundle object
%
% Outputs:
%    V - vertices object
%
% Example: 
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: intervallhull,  polytope

% Author:       Matthias Althoff
% Written:      18-August-2016 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%obtain polytope
P = polytope(obj);

%obtain vertices
V = vertices(P);

%------------- END OF CODE --------------