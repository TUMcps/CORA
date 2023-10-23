function probZ = project(probZ,dims)
% project - projects a probabilistic zonotope onto the specified dimensions
%
% Syntax:
%    probZ = project(probZ,dims)
%
% Inputs:
%    probZ - (probZonotope) probabilistic zonotope
%    dims - dimensions for projection
%
% Outputs:
%    probZ - (probZonotope) projected probabilistic zonotope
%
% Example: 
%    Z1 = [10 1 -2; 0 1 1; -1 0 3];
%    Z2 = [0.6 1.2 1; 0.6 -1.2 0; 0.3 1 -0.7];
%    probZ = probZonotope(Z1,Z2,2);
%    probZ = project(probZ, [1 3])
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       17-July-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% project set-based part
probZ.Z = probZ.Z(dims,:);

% project probabilistic part
probZ.cov = probZ.cov(dims,dims);

% ------------------------------ END OF CODE ------------------------------
