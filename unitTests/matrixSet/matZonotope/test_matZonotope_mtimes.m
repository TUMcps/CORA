function res = test_matZonotope_mtimes
% test_matZonotope_mtimes - unit test function for mtimes
% 
% Syntax:
%    res = test_matZonotope_mtimes
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
% Written:       25-April-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

% init some matrix zonotopes
% scalar
C = 0;
G = []; G(:,:,1) = 1; G(:,:,2) = -2;
matZscalar = matZonotope(C,G);
% nx1 vector
C = [0; 1; 1];
G = []; G(:,:,1) = [1; -1; -2]; G(:,:,2) = [-2; 0; 1];
matZvector = matZonotope(C,G);
% matrix
C = [0 2; 1 -1; 1 -2];
G = []; G(:,:,1) = [1 1; -1 0; -2 1]; G(:,:,2) = [-2 0; 0 1; 1 -1];
matZ = matZonotope(C,G);

% scalar case
factor = 2;
matZres = factor * matZscalar;
matZres = matZscalar * factor ;
resvec(end+1) = all(dim(matZres) == [1,1]);
matZres = factor * matZvector;
matZres = matZvector * factor ;
resvec(end+1) = all(dim(matZres) == [3,1]);
matZres = factor * matZ;
matZres = matZ * factor;
resvec(end+1) = all(dim(matZres) == [3,2]);

% matrx case
factor = magic(3);
matZres = factor * matZscalar;
matZres = matZscalar * factor ;
resvec(end+1) = all(dim(matZres) == [3,3]);
matZres = factor * matZvector;
resvec(end+1) = all(dim(matZres) == [3,1]);
matZres = factor * matZ;
resvec(end+1) = all(dim(matZres) == [3,2]);

% matZ * matZ
matZres = matZ.' * matZvector;
resvec(end+1) = all(dim(matZres) == [2,1]);

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
