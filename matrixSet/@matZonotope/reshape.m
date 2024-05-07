function matZ_reshaped = reshape(matZ,n,m)
% reshape - reshapes the center and the generators
%
% Syntax:
%    matZ_transposed = reshape(matZ,n,m)
%
% Inputs:
%    matZ - matZonotope object
%    n - new number of rows
%    m - new number of columns
%
% Outputs:
%    matZ_reshaped - matZonotope object
%
% Example: 
%    C = [0 0; 1 0];
%    G = []; G(:,:,1) = [1 3; -1 2]; G(:,:,2) = [2 0; 1 -1];
%    matZ = matZonotope(C,G);
%    matZ_transposed = reshape(matZ,1,4);
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

if nargin < 3
    throw(CORAerror('CORA:notEnoughInputArgs',3))
elseif nargin > 3
    throw(CORAerror('CORA:tooManyInputArgs',3))
end

matZ_reshaped = matZ;
matZ_reshaped.C = reshape(matZ.C,n,m);
[~,~,h] = size(matZ.G);
matZ_reshaped.G = reshape(matZ.G,n,m,h);

end

% ------------------------------ END OF CODE ------------------------------
