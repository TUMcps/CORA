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

% Authors:       Matthias Althoff
% Written:       15-September-2008
% Last update:   20-October-2023 (TL, correct projection for G\in\R^{n x 0})
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

Z.c = Z.c(dims,:);
if size(Z.G,1) > 0
    Z.G = Z.G(dims,:);
end

% ------------------------------ END OF CODE ------------------------------
