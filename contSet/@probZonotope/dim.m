function n = dim(probZ)
% dim - returns the dimension of the ambient space of a probabilistic
%    zonotope
%
% Syntax:
%    n = dim(probZ)
%
% Inputs:
%    probZ - probZonotope object
%
% Outputs:
%    n - dimension of the ambient space
%
% Example: 
%    Z1 = [10 1 -2; 0 1 1];
%    Z2 = [0.6 1.2; 0.6 -1.2];
%    probZ = probZonotope(Z1,Z2,2);
%    n = dim(probZ);
%
% Other m-files required: center
% Subfunctions: none
% MAT-files required: none
%
% See also: rank.m

% Authors:       Mark Wetzlinger
% Written:       10-June-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

n = length(center(probZ));

% ------------------------------ END OF CODE ------------------------------
