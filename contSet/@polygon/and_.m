function pgon = and_(pgon1, pgon2, varargin)
% and_ - computes the intersection of two polygons
%
% Syntax:
%    pgon = and_(pgon1, pgon2)
%
% Inputs:
%    pgon1 - polygon
%    pgon2 - polygon
%
% Outputs:
%    pgon - polygon
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/and

% Authors:       Niklas Kochdumper
% Written:       13-March-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% use polyshape/intersect
pgon = polygon(intersect(pgon1.set, pgon2.set));

end

% ------------------------------ END OF CODE ------------------------------
