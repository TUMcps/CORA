function n = dim(cZ)
% dim - returns the dimension of the ambient space of a constrained zonotope
%
% Syntax:
%    n = dim(cZ)
%
% Inputs:
%    cZ - conZonotope object
%
% Outputs:
%    n - dimension of the ambient space
%
% Example: 
%    Z = [0 1 0 1;0 1 2 -1];
%    A = [-2 1 -1]; b = 2;
%    cZ = conZonotope(Z,A,b);
% 
%    n = dim(cZ)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       15-September-2019 
% Last update:   14-March-2021 (MW, different approach)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

n = size(cZ.c,1);

% ------------------------------ END OF CODE ------------------------------
