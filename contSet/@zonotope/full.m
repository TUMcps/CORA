function Z_full = full(Z)
% full - Converts a zonotope (center and generator) into full
%          representation (from sparse)
%
% Syntax:
%    Z_full = full(Z)
%
% Inputs:
%    Z - zonotope object (assumed to be sparse)
%
% Outputs:
%    Z_full - zonotope with full center vector and generator matrix
%
% Example: 
%    c = rand(10,1);
%    c(c < 0.8) = 0;
%    G = rand(10,10);
%    G(G < 0.8) = 0;
%    Z = zonotope(sparse(c),sparse(G));
%    Z_full = full(Z);
%    disp([Z.c Z.G]);
%    disp([Z_full.c Z_full.G]);
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

c_full = full(center(Z));
G_full = full(generators(Z));

Z_full = zonotope(c_full,G_full);

% ------------------------------ END OF CODE ------------------------------
