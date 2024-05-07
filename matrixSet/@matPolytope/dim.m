function n = dim(matP,rc)
% dim - returns the dimension of the matrix polytope
%
% Syntax:
%    n = dim(matP,rc)
%
% Inputs:
%    matP - matPolytope object
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

% Authors:       Mark Wetzlinger, Tobias Ladner
% Written:       03-April-2023
% Last update:   02-May-2024 (TL, simplified due to new structure of V)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
if nargin == 1
    rc = 1:2;
end

n = size(matP.V,rc);

% ------------------------------ END OF CODE ------------------------------
