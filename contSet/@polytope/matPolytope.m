function matP = matPolytope(P)
% matPolytope - converts the given polytope into a matPolytope object
%
% Syntax:
%    matP = matPolytope(P)
%
% Inputs:
%    P - polytope object
%
% Outputs:
%    matP - matPolytope
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Matthias Althoff, Tobias Ladner, Mark Wetzlinger
% Written:       21-June-2010
% Last update:   02-May-2024 (TL, moved out of constructor)
%                12-July-2024 (MW, fix broken function)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% get vertices from polytope class
V = vertices(P);

% rewrite vertices in n x 1 x N format (n = dimension, N = number of
% vertices)
V = reshape(V,[dim(P),1,size(V,2)]);

% init matrix polytope
matP = matPolytope(V);

% ------------------------------ END OF CODE ------------------------------
