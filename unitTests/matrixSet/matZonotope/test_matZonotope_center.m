function res = test_matZonotope_center
% test_matZonotope_center - unit test function for dimension center
% 
% Syntax:
%    res = test_matZonotope_center
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

% empty matrix zonotope
matZempty = matZonotope();
assert(isempty(center(matZempty)));

% scalar
C = 0;
G = []; G(:,:,1) = 1; G(:,:,2) = -2;
matZscalar = matZonotope(C,G);
assert(compareMatrices(C,center(matZscalar)));

% nx1 vector
C = [0; 1; 1];
G = []; G(:,:,1) = [1; -1; -2]; G(:,:,2) = [-2; 0; 1];
matZvector = matZonotope(C,G);
assert(compareMatrices(C,center(matZvector)));

% matrix
C = [0 2; 1 -1; 1 -2];
G = []; G(:,:,1) = [1 1; -1 0; -2 1]; G(:,:,2) = [-2 0; 0 1; 1 -1];
matZ = matZonotope(C,G);
assert(compareMatrices(C,center(matZ)));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
