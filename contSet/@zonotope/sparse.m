function Z_sparse = sparse(Z)
% sparse - Converts a zonotope (center and generator) into sparse
%          representation
%
% Syntax:
%    Z_sparse = sparse(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    Z_sparse - zonotope with sparse center vector and generator matrix
%
% Example: 
%    c = rand(10,1);
%    c(c < 0.8) = 0;
%    G = rand(10,10);
%    G(G < 0.8) = 0;
%    Z = zonotope(c,G);
%    Z_sparse = sparse(Z);
%    disp([Z.c Z.G]);
%    disp([Z_sparse.c Z_sparse.G]);
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Maximilian Perschl
% Written:       02-June-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

c_sparse = sparse(center(Z));
G_sparse = sparse(generators(Z));

Z_sparse = zonotope(c_sparse,G_sparse);

% ------------------------------ END OF CODE ------------------------------
