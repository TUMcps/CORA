function P = polytope(pgon)
% polytope - convert a convex polygon to a polytope object
%
% Syntax:
%    P = polytope(pgon)
%
% Inputs:
%    pgon - polygon
%
% Outputs:
%    P - polytope
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       13-March-2020
% Last update:   11-October-2024 (TL, rewrote)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% compute convex hull
pgon_conv = convHull(pgon);

% get vertices
V_conv = vertices(pgon_conv);

% init polytope
P = polytope(V_conv);

end

% ------------------------------ END OF CODE ------------------------------
