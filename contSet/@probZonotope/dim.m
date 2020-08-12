function d = dim(probZ)
% dim - return dimension of probabilistic zonotope
%
% Syntax:  
%    d = dim(probZ)
%
% Inputs:
%    probZ - probabilistic zonotope object
%
% Outputs:
%    d - dimension of probZ
%
% Example: 
%    Z1 = [10 1 -2; 0 1 1];
%    Z2 = [0.6 1.2; 0.6 -1.2];
%    probZ = probZonotope(Z1,Z2,2);
%    d = dim(probZ);
%
% Other m-files required: center
% Subfunctions: none
% MAT-files required: none
%
% See also: rank.m
%
% Author:        Mark Wetzlinger
% Written:       10-June-2020
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

d = length(center(probZ));

%------------- END OF CODE --------------