function V = vertices_(zB,varargin)
% vertices_ - Returns potential vertices of a zonotope bundle
%    WARNING: Do not use this function for high order zonotope bundles as
%    the computational complexity grows exponentially!
%
% Syntax:  
%    V = vertices_(zB)
%
% Inputs:
%    zB - zonoBundle object
%
% Outputs:
%    V - matrix storing the vertices 
%
% Example: 
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polytope

% Author:       Matthias Althoff
% Written:      18-August-2016 
% Last update:  ---
% Last revision:27-March-2023 (MW, rename vertices_)

%------------- BEGIN CODE --------------

%obtain polytope
P = mptPolytope(zB);

%obtain vertices (input check in polytope-function)
V = vertices(P);

%------------- END OF CODE --------------