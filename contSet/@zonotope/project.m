function Z = project(Z,dims)
% project - projects a zonotope onto the specified dimensions
%
% Syntax:  
%    Z = project(Z,dims)
%
% Inputs:
%    Z - (zonotope) zonotope
%    dims - dimensions for projection
%
% Outputs:
%    Z - (zonotope) projected zonotope
%
% Example: 
%    Z = zonotope([1 -1 0; 0 0 -1; 1 0 1]);
%    Z = project(Z, [1 3])
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      15-September-2008
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

Z.Z = Z.Z(dims,:);

%------------- END OF CODE --------------