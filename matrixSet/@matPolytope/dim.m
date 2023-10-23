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

% Authors:       Mark Wetzlinger
% Written:       03-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if nargin == 1
    if matP.verts == 0
        n = [0,0];
    else
        % read dimension from first vertex
        n = size(matP.vertices{1});
    end

else
    if matP.verts == 0
        n = 0;
    else
        % read row or column dimension
        n = size(matP.vertices{1},rc);
    end

end

% ------------------------------ END OF CODE ------------------------------
