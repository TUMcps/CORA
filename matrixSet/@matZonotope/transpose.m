function matZ_transposed = transpose(matZ)
% transpose - transposes the center and the generators
%
% Syntax:
%    matZ_transposed = transpose(matZ)
%
% Inputs:
%    matZ - matZonotope object
%
% Outputs:
%    matZ_transposed - matZonotope object
%
% Example: 
%    C = [0 0; 1 0];
%    G = []; G(:,:,1) = [1 3; -1 2]; G(:,:,2) = [2 0; 1 -1];
%    matZ = matZonotope(C,G);
%    matZ_transposed = transpose(matZ);
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

matZ_transposed = matZ;
matZ_transposed.C = matZ.C';
matZ_transposed.G = pagetranspose(matZ.G);

end

% ------------------------------ END OF CODE ------------------------------
