function N = numverts(matP)
% numverts - returns the number of vertices of a matPolytope
%
% Syntax:
%    N = numverts(matP)
%
% Inputs:
%    matZ - matZonotope object
%
% Outputs:
%    N - numeric, number of vertices
%
% Example:
%    V(:,:,1) = [1 2; 0 1];
%    V(:,:,2) = [1 3; -1 2];
%    matP = matPolytope(V);
%    N = numverts(matP)
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

% V has dimensions (n x m x h)
N = size(matP.V,3);

end

% ------------------------------ END OF CODE ------------------------------
