function n = dim(matZ,rc)
% dim - returns the dimension of the matrix zonotope
%
% Syntax:
%    n = dim(matZ)
%    n = dim(matZ,rc)
%
% Inputs:
%    matZ - matZonotope object
%    rc - 1 for row dimension, 2 for column dimension
%
% Outputs:
%    n - array with row and column dimension
%
% Example: 
%    C = [0 0; 0 0];
%    G{1} = [1 3; -1 2]; G{2} = [2 0; 1 -1];
%    matZ = matZonotope(C,G);
%    n = dim(matZ);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Mark Wetzlinger
% Written:       03-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if nargin == 1
    % read dimension from center
    n = size(matZ.center);

else
    % read dimension from center
    n = size(matZ.center,rc);

end

% ------------------------------ END OF CODE ------------------------------
