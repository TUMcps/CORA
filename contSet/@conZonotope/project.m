function obj = project(obj,dim)
% project - project a constrained zonotope object to a subspace
%
% Syntax:  
%    obj = project(obj,dim)
%
% Inputs:
%    obj - constrained zonotope object
%    dim - dimensions of the projection
%
% Outputs:
%
% Example: 
%    Z = [0 1 0 1 0;0 1 2 -1 0;0 0 0 0 1];
%    A = [-2 1 -1 0];
%    b = 2;
%    cZono = conZonotope(Z,A,b);
%    projZono = project(cZono,[1,2]);
%    
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Niklas Kochdumper
% Written:      11-May-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

obj.Z = obj.Z(dim,:);

%------------- END OF CODE --------------