function res = isConvex(pgon)
% isConvex - check if the polygon is convex
%
% Syntax:
%    res = isConvex(pgon)
%
% Inputs:
%    pgon - polygon
%
% Outputs:
%    res - logical
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Niklas Kochdumper
% Written:       13-March-2020
% Last update:   11-October-2024 (TL, rewrote)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% polygon is convex if its convex hull has the same number of
% vertices as the polygon itself already has
res = size(vertices_(convHull(pgon)),2) == size(vertices_(pgon),2);

end

% ------------------------------ END OF CODE ------------------------------
