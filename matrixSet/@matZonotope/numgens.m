function h = numgens(matZ)
% numgens - returns the number of generators of a matZonotope
%
% Syntax:
%    h = numgens(matZ)
%
% Inputs:
%    matZ - matZonotope object
%
% Outputs:
%    h - numeric, number of generators
%
% Example:
%    C = eye(2);
%    G(:,:,1) = eye(1);
%    G(:,:,2) = eye(2);
%    matZ = matZonotope(C,G);
%    h = numgens(matZ)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Tobias Ladner
% Written:       25-April-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% G has dimensions (n x m x h)
G = matZ.G;
h = size(G,3);

end

% ------------------------------ END OF CODE ------------------------------
