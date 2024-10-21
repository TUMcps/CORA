function res = test_matZonotope_transpose
% test_matZonotope_transpose - unit test function for transpose
% 
% Syntax:
%    res = test_matZonotope_transpose
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       26-April-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% scalar
C = 0;
G = []; G(:,:,1) = 1; G(:,:,2) = -2;
matZscalar = matZonotope(C,G);
assert(all(dim(transpose(matZscalar)) == [1 1]));

% nx1 vector
C = [0; 1; 1];
G = []; G(:,:,1) = [1; -1; -2]; G(:,:,2) = [-2; 0; 1];
matZvector = matZonotope(C,G);
assert(all(dim(transpose(matZvector)) == [1 3]));

% matrix
C = [0 2; 1 -1; 1 -2];
G = []; G(:,:,1) = [1 1; -1 0; -2 1]; G(:,:,2) = [-2 0; 0 1; 1 -1];
matZ = matZonotope(C,G);
assert(compareMatrices(C,center(matZ)));
assert(all(dim(transpose(matZ)) == [2 3]));

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
